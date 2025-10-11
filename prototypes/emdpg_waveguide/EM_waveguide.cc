/**
 * @brief Necessary includes for the implementation of the DPG method for the
 * time-harmonic maxwell equations.
 *
 * @author Oreste Marquis, Bruno Blais, Matthias Maier, 2025
 */

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec_sz.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_faces.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <complex>
#include <fstream>
#include <iostream>

/**
 * @namespace LA
 * @brief Alias namespace for dealii::LinearAlgebraTrilinos.
 *
 * This namespace provides a convenient shorthand for accessing the
 * Trilinos-based linear algebra functionality from the deal.II library.
 * All classes and functions from dealii::LinearAlgebraTrilinos can be
 * accessed via the LA namespace.
 */
namespace LA
{
  using namespace dealii::LinearAlgebraTrilinos;
} // namespace LA
namespace EM_DPG
{
  using namespace dealii;

  // This is a helper function that maps a 3D tensor to a 2D face in
  // H^-1/2(curl).
  template <int dim>
  DEAL_II_ALWAYS_INLINE inline Tensor<1, dim, std::complex<double>>
  map_H12(const Tensor<1, dim, std::complex<double>> &tensor,
          const Tensor<1, dim>                       &normal)
  {
    auto result = cross_product_3d(normal, cross_product_3d(tensor, normal));

    return result;
  }

  // We create a class to handle the input parameters.
  template <int dim>
  class Parameters : public ParameterAcceptor
  {
  public:
    // CONSTRUCTOR
    Parameters();

    // Constexpr values
    static constexpr double               speed_of_light = 299792458.;
    static constexpr std::complex<double> imag{0., 1.};

    // Material properties
    std::complex<double> epsilon_r;
    std::complex<double> mu_r;
    double               sigma_r;

    // Waveguide geometry
    double waveguide_a; // Width of the waveguide
    double waveguide_b; // Height of the waveguide
    double waveguide_length;

    // Frequency and excitation modes
    double frequency;
    int    mode_x;
    int    mode_y;

    // Mesh parameters
    unsigned int              n_refinements;
    std::vector<unsigned int> n_repetitions;
    unsigned int              initial_refinements;
    bool                      adaptive_refinement;
    int                       max_refinement_level;

    // Order of the DGQ elements
    unsigned int degree;
    unsigned int delta_degree;
  };

  template <int dim>
  Parameters<dim>::Parameters()
    : ParameterAcceptor("Parameters")
  {
    // Defaults properties are the ones of vacuum
    epsilon_r = {1., 0.};
    mu_r      = {1., 0.};
    sigma_r   = 0.;

    add_parameter("epsilon r", epsilon_r, "Relative permittivity of medium");

    add_parameter("mu r", mu_r, "Relative permeability of medium");

    add_parameter("sigma r", sigma_r, "Relative conductivity of medium");

    // Default waveguide geometry is 0.25 x 0.25 x 1
    waveguide_a      = 0.25; // Width of the waveguide
    waveguide_b      = 0.25; // Height of the waveguide
    waveguide_length = 1.0;  // Length of the waveguide

    add_parameter("waveguide a", waveguide_a, "Width of the waveguide");

    add_parameter("waveguide b", waveguide_b, "Height of the waveguide");

    add_parameter("waveguide length",
                  waveguide_length,
                  "Length of the waveguide");

    // Default frequency and excitation modes
    frequency = 670356315.23; // This corresponds to the TE10 mode with a
                              // wavelenght of 1 for a = 0.25
    mode_x = 1;               // Mode in x direction
    mode_y = 0;               // Mode in y direction

    add_parameter("frequency", frequency, "Frequency of the signal");

    add_parameter("mode x", mode_x, " TE mode in x direction");

    add_parameter("mode y", mode_y, " TE mode in y direction");

    // Default mesh parameters
    n_refinements        = 0;
    initial_refinements  = 0;
    n_repetitions        = {1, 1, 4};
    adaptive_refinement  = false;
    max_refinement_level = 8;

    add_parameter("n refinements",
                  n_refinements,
                  "Number of refinements cycle for mesh refinement");

    add_parameter("initial refinements",
                  initial_refinements,
                  "Number of refinements for the coarse mesh before starting");

    add_parameter("n repetitions",
                  n_repetitions,
                  "Number of repetitions for the subdivided hyper rectangle");

    add_parameter("adaptive refinement",
                  adaptive_refinement,
                  "If true, the code will use adaptive mesh refinement, "
                  "otherwise it will use uniform mesh refinement");

    add_parameter("max refinement level",
                  max_refinement_level,
                  "Maximum number of refinements from the coarse mesh");

    // Default finite element parameters
    degree = 0;
    add_parameter("degree", degree, "order of the DGQ finite element space");

    delta_degree = 1;
    add_parameter("delta degree",
                  delta_degree,
                  "order difference between the trial and test space");
  }

  /* Create analytical solution class for the electric field (E).
   * The corresponding solution for the electric field follows the physicist's
   * time convention with the imaginary unit in the exponent being negative
   * (i.e. exp(-i wt) ). The function manage both the real and complex
   * components of the solution. Since the frequency, the material properties,
   * and the geometry are arguments for the parameter file, they are used in the
   * constructor to initialize the different wavenumber. Note that everything is
   * addimensionalize using the following : \tilde{\mathbf{E}} =
   * \frac{\mathbf{E}}{E_0}, \tilde{\mathbf{H}} = \frac{Z_0}{E_0} \mathbf{H},
   *       \varepsilon_r = \frac{\varepsilon}{\varepsilon_0},
   *       \mu_r = \frac{\mu}{\mu_0},
   *       Z_r = \frac{Z}{Z_0},
   *       \sigma_r = \frac{\sigma}{\omega \varepsilon_0},
   *       \tilde{\omega} =  \frac{\omega}{k_0 c_0},
   *       \tilde{\mathbf{J}} = \frac{1}{\omega \varepsilon E_0 } \mathbf{J}
   *       \tilde{k} = L \mathbf{k},
   *       \tilde{\nabla} = L \nabla.
   */

  template <int dim>
  class AnalyticalSolution_E : public Function<dim>
  {
  public:
    // Constructor
    AnalyticalSolution_E(const Parameters<dim> &parameters)
      : Function<dim>(6)
      , parameters(parameters)
      , k(2. * M_PI * parameters.frequency / parameters.speed_of_light)
      , k_x(parameters.mode_x * M_PI / parameters.waveguide_a)
      , k_y(parameters.mode_y * M_PI / parameters.waveguide_b)
      , k_z(std::sqrt(k * k - k_x * k_x - k_y * k_y))
      , k_c(std::sqrt(k_x * k_x + k_y * k_y))
      , Z_r(std::sqrt(parameters.mu_r / parameters.epsilon_r))
    {}

    virtual double
    value(const Point<dim> &p, const unsigned int component) const override;

  private:
    const Parameters<dim>     &parameters;
    const double               k;
    const double               k_x;
    const double               k_y;
    const double               k_z;
    const double               k_c;
    const std::complex<double> Z_r;
  };

  template <int dim>
  double
  AnalyticalSolution_E<dim>::value(const Point<dim>  &p,
                                   const unsigned int component) const
  {
    const std::complex<double> &imag = parameters.imag;
    std::complex<double>        factor =
      imag * k * Z_r / (k_c * k_c) * std::exp(imag * k_z * p[2]);

    if (component == 0)
      return (-factor * k_y * std::cos(k_x * p[0]) * std::sin(k_y * p[1]))
        .real();
    else if (component == 1)
      return (factor * k_x * std::sin(k_x * p[0]) * std::cos(k_y * p[1]))
        .real();
    else if (component == 2)
      return 0.;
    else if (component == 3)
      return (-factor * k_y * std::cos(k_x * p[0]) * std::sin(k_y * p[1]))
        .imag();
    else if (component == 4)
      return (factor * k_x * std::sin(k_x * p[0]) * std::cos(k_y * p[1]))
        .imag();
    else if (component == 5)
      return 0.;
    else
      throw std::runtime_error(
        "Too many components for the analytical solution");
  }

  // Same here, we create the analytical solution class for the magnetic field
  // (H).
  template <int dim>
  class AnalyticalSolution_H : public Function<dim>
  {
  public:
    // Constructor
    AnalyticalSolution_H(const Parameters<dim> &parameters)
      : Function<dim>(6)
      , parameters(parameters)
      , k(2. * M_PI * parameters.frequency / parameters.speed_of_light)
      , k_x(parameters.mode_x * M_PI / parameters.waveguide_a)
      , k_y(parameters.mode_y * M_PI / parameters.waveguide_b)
      , k_z(std::sqrt(k * k - k_x * k_x - k_y * k_y))
      , k_c(std::sqrt(k_x * k_x + k_y * k_y))
    {}

    virtual double
    value(const Point<dim> &p, const unsigned int component) const override;

  private:
    const Parameters<dim> &parameters;
    const double           k;
    const double           k_x;
    const double           k_y;
    const double           k_z;
    const double           k_c;
  };

  template <int dim>
  double
  AnalyticalSolution_H<dim>::value(const Point<dim>  &p,
                                   const unsigned int component) const
  {
    const std::complex<double> &imag   = parameters.imag;
    std::complex<double>        factor = std::exp(imag * k_z * p[2]);

    if (component == 0)
      return (-factor * imag * k_z * k_x / (k_c * k_c) * std::sin(k_x * p[0]) *
              std::cos(k_y * p[1]))
        .real();
    else if (component == 1)
      return (-factor * imag * k_z * k_y / (k_c * k_c) * std::cos(k_x * p[0]) *
              std::sin(k_y * p[1]))
        .real();
    else if (component == 2)
      return (factor * std::cos(k_x * p[0]) * std::cos(k_y * p[1])).real();
    else if (component == 3)
      return (-factor * imag * k_z * k_x / (k_c * k_c) * std::sin(k_x * p[0]) *
              std::cos(k_y * p[1]))
        .imag();
    else if (component == 4)
      return (-factor * imag * k_z * k_y / (k_c * k_c) * std::cos(k_x * p[0]) *
              std::sin(k_y * p[1]))
        .imag();
    else if (component == 5)
      return (factor * std::cos(k_x * p[0]) * std::cos(k_y * p[1])).imag();
    else
      throw std::runtime_error(
        "Too many components for the analytical solution");
  }

  // Main class for the DPG solver. The method rely on multiple DoFHandler and
  // FESystem to manage the different functional spaces. The DoFHandlers that we
  // rely on are the following:
  // - The <code>dof_handler_trial_interior</code> is for the unknowns in the
  // interior of the cells;
  // - The <code>dof_handler_trial_skeleton</code> is for the unknowns in the
  // skeleton;
  // - The <code>dof_handler_test</code> is for the test functions. Although we
  // do not use the unknowns associated with this DoFHandler, it enables us to
  // evaluate the test function we will use in DPG.

  // The same applies for the three FESystem:
  // <code>fe_system_trial_interior</code>,
  // <code>fe_system_trial_skeleton</code> and <code>fe_system_test</code>. In
  // each one of these, we will store the relevant finite element space in the
  // same order to avoid confusion. The first component will therefore always be
  // related to the real part of the electric field, the second component to the
  // its imaginary part, the third component to the real part of the magnetic
  // field and the fourth component to its imaginary part.
  template <int dim>
  class Time_Harmonic_Maxwell_DPG : public ParameterAcceptor
  {
  public:
    Time_Harmonic_Maxwell_DPG(const Parameters<dim> &parameters);
    void
    run();

  private:
    // <code>make_grid</code> creates the initial mesh.
    void
    make_grid(const unsigned int initial_refinement);

    // <code>refine_grid</code> refines the mesh uniformly or adaptively using
    // the built-in error estimator depending on the input parameters. Note that
    // at the moment the adaptive refinement does not work in parallel for
    // fe_nedelecSZ elements. From quick tests, it seems to be related to
    // hanging nodes constraints. Use the standard Nedelec elements if you want
    // to use adaptive refinement in parallel.
    void
    refine_grid();

    // The <code>setup_system</code> function initializes the three DoFHandlers,
    // the system matrix and right-hand side and establishes the boundary
    // conditions that rely on AffineConstraints
    void
    setup_system();

    // The <code>assemble_system</code> assembles both the right-hand side and
    // the system matrix. This function is used twice per problem and it has
    // two functions.
    // - When <code>solve_interior = false</code> the system is assembled and is
    // locally condensed such that the resulting system only contains the
    // skeleton unknowns. This is achieved by local condensation.
    // - When <code>solve_interior = true</code> the system is assembled and the
    // skeleton degrees of freedom solution is used to
    // reconstruct the interior solution and the built-in error estimator right
    // after.
    void
    assemble_system(bool solve_interior = false);

    // <code>solve_linear_system_skeleton</code> solves the linear system of
    // equations for the skeleton degree of freedom.
    void
    solve_skeleton();

    // <code>output_results</code> write the skeleton and the interior unknowns
    // into two different VTU and PVTU files.
    void
    output_results(unsigned int cycle);

    // <code> calculate_L2_error</code> calculates the $L^2$ norm of the error
    // using the analytical solution to verify convergence order, but it also
    // calculates the error in the energy norm for each cell to obtain the
    // h-refinement indicator.
    void
    calculate_L2_error();

    /* run time parameters */
    const Parameters<dim>     &parameters;
    const double               k;               // Wavenumber
    const double               k_z;             // z-component of the wavenumber
    const std::complex<double> Z_r;             // Relative impedance
    std::complex<double>       imag = {0., 1.}; // Imaginary unit

    // MPI-related variables
    MPI_Comm           mpi_communicator;
    ConditionalOStream pcout;
    TimerOutput        computing_timer;

    // The parallel triangulation
    parallel::distributed::Triangulation<dim> triangulation;

    // The variables for the interior unknowns.
    const FESystem<dim> fe_trial_interior;
    DoFHandler<dim>     dof_handler_trial_interior;
    LA::MPI::Vector     locally_relevant_solution_interior;
    IndexSet            locally_owned_dofs_trial_interior;
    IndexSet            locally_relevant_dofs_trial_interior;

    // The variables for the skeleton and, consequently, the linear system.
    const FESystem<dim>       fe_trial_skeleton;
    DoFHandler<dim>           dof_handler_trial_skeleton;
    LA::MPI::Vector           locally_relevant_solution_skeleton;
    IndexSet                  locally_owned_dofs_trial_skeleton;
    IndexSet                  locally_relevant_dofs_trial_skeleton;
    LA::MPI::Vector           system_rhs;
    LA::MPI::SparseMatrix     system_matrix;
    AffineConstraints<double> constraints;

    // The variables for the test spaces.
    const FESystem<dim> fe_test;
    DoFHandler<dim>     dof_handler_test;
    LA::MPI::Vector     locally_relevant_residual_estimator;
    IndexSet            locally_owned_dofs_test;
    IndexSet            locally_relevant_dofs_test;

    // The vector container for the energy norm on each cell used for adaptive
    // refinement.
    Vector<double> estimated_error_per_cell;

    // A container for the $L^2$ error and other related quantities.
    ConvergenceTable error_table;

    // Extractors that will be used at multiple places in the class to
    // select the relevant components for the calculation. Those are created
    // here to avoid repetition throughout the class.
    const FEValuesExtractors::Vector extractor_E_real;
    const FEValuesExtractors::Vector extractor_E_imag;
    const FEValuesExtractors::Vector extractor_H_real;
    const FEValuesExtractors::Vector extractor_H_imag;
  };

  // DPG constructor
  // In the constructor, we assign the relevant finite element to each FESystem
  // following the nomenclature described above in the class description:
  // - <code>fe_system_trial_interior</code> contains $\Re(\mathbf{E})$,
  // $\Im(\mathbf{E})$, $\Re(\mathbf{H})$, $\Im(\mathbf{H})$ ;
  // - <code>fe_system_trial_skeleton</code> contains $\Re(\hat{E}_\parallel)$,
  // $\Im(\hat{E}_\parallel)$, $\Re(\hat{H}_\parallel)$,
  // $\Im(\hat{H}_\parallel)$;
  // - <code>fe_system_test</code> contains $\Re(\mathbf{F})$,
  // $\Im(\mathbf{F})$, $\Re(\mathbf{I})$, $\Im(\mathbf{I})$.
  // We also initialize probleme dependent parameters and parallel related
  // variables.

  template <int dim>
  Time_Harmonic_Maxwell_DPG<dim>::Time_Harmonic_Maxwell_DPG(
    const Parameters<dim> &parameters)
    : ParameterAcceptor("Time_Harmonic_Maxwell_DPG")
    , parameters(parameters)
    , k(2. * M_PI * parameters.frequency / parameters.speed_of_light)
    , k_z(std::sqrt(
        k * k - std::pow(parameters.mode_x * M_PI / parameters.waveguide_a, 2) -
        std::pow(parameters.mode_y * M_PI / parameters.waveguide_b, 2)))
    , Z_r(std::sqrt(parameters.mu_r / parameters.epsilon_r))
    , mpi_communicator(MPI_COMM_WORLD)
    , pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    , computing_timer(mpi_communicator,
                      pcout,
                      TimerOutput::summary,
                      TimerOutput::wall_times)
    , triangulation(mpi_communicator,
                    typename Triangulation<dim>::MeshSmoothing(
                      Triangulation<dim>::smoothing_on_refinement |
                      Triangulation<dim>::smoothing_on_coarsening))
    , fe_trial_interior(FE_DGQ<dim>(parameters.degree) ^ dim,
                        FE_DGQ<dim>(parameters.degree) ^ dim,
                        FE_DGQ<dim>(parameters.degree) ^ dim,
                        FE_DGQ<dim>(parameters.degree) ^ dim)
    , // (E, E_imag, H, H_imag)
    dof_handler_trial_interior(triangulation)
    , fe_trial_skeleton(FE_NedelecSZ<dim>(parameters.degree),
                        FE_NedelecSZ<dim>(parameters.degree),
                        FE_NedelecSZ<dim>(parameters.degree),
                        FE_NedelecSZ<dim>(parameters.degree))
    , // (E_hat, E_hat_imag, H_hat, H_hat_imag)
    dof_handler_trial_skeleton(triangulation)
    , fe_test(FE_NedelecSZ<dim>(parameters.degree + parameters.delta_degree),
              FE_NedelecSZ<dim>(parameters.degree + parameters.delta_degree),
              FE_NedelecSZ<dim>(parameters.degree + parameters.delta_degree),
              FE_NedelecSZ<dim>(parameters.degree + parameters.delta_degree))
    , // (F, F_imag, I, I_imag)
    dof_handler_test(triangulation)
    , extractor_E_real(0)
    , extractor_E_imag(dim)
    , extractor_H_real(2 * dim)
    , extractor_H_imag(3 * dim)
  {
    AssertThrow(dim == 3,
                ExcMessage("The implementation only works for dim==3."));

    AssertThrow(parameters.delta_degree >= 1,
                ExcMessage("The delta_degree needs to be at least 1."));

    AssertThrow(k > 0,
                ExcMessage("The wavenumber must be positive, its the "
                           "magnitude of the wave vector."));
  }

  // Here the <code>make_grid</code> function creates the initial rectangular
  // mesh using a structured grid using
  // dealii::GridGenerator::subdivided_hyper_rectangle.
  template <int dim>
  void
  Time_Harmonic_Maxwell_DPG<dim>::make_grid(
    const unsigned int initial_refinement)
  {
    pcout << "*--- Making grid ---*" << std::endl;

    Point<dim> lower_left_waveguide(0., 0., 0.);
    Point<dim> upper_right_waveguide(parameters.waveguide_a,
                                     parameters.waveguide_b,
                                     parameters.waveguide_length);

    GridGenerator::subdivided_hyper_rectangle(triangulation,
                                              parameters.n_repetitions,
                                              lower_left_waveguide,
                                              upper_right_waveguide,
                                              true);
    triangulation.refine_global(initial_refinement);
  }

  // Here the <code>refine_grid</code> function refines the mesh uniformly or
  // adaptively using the built-in error estimator computed from the previous
  // solution.
  template <int dim>
  void
  Time_Harmonic_Maxwell_DPG<dim>::refine_grid()
  {
    pcout << "*--- Refining grid ---*" << std::endl;

    TimerOutput::Scope t(computing_timer, "refine_grid");

    if (parameters.adaptive_refinement)
      {
        parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
          triangulation, estimated_error_per_cell, 0.5, 0.3);

        // We added the possibility to limit the maximum refinement level.
        for (auto &cell : triangulation.active_cell_iterators())
          {
            if (cell->is_locally_owned())
              {
                if (cell->level() >= parameters.max_refinement_level &&
                    cell->refine_flag_set())
                  cell->clear_refine_flag();
              }
          }

        triangulation.execute_coarsening_and_refinement();
      }
    else
      {
        triangulation.refine_global();
      }
  }

  // In this function, we initialize the three different DoFHandlers and apply
  // the constraints related to the boundary conditions. We also create the
  // sparsity pattern and initialize the system matrix and right-hand side.
  template <int dim>
  void
  Time_Harmonic_Maxwell_DPG<dim>::setup_system()
  {
    TimerOutput::Scope t(computing_timer, "setup_system");

    dof_handler_trial_skeleton.distribute_dofs(fe_trial_skeleton);
    dof_handler_trial_interior.distribute_dofs(fe_trial_interior);
    dof_handler_test.distribute_dofs(fe_test);

    unsigned int n_dofs_interior = dof_handler_trial_interior.n_dofs();
    unsigned int n_dofs_skeleton = dof_handler_trial_skeleton.n_dofs();
    unsigned int n_dofs_test     = dof_handler_test.n_dofs();

    pcout << "*--- Setting up system ---*" << std::endl;
    pcout << std::endl
          << "Number of dofs for the interior: " << n_dofs_interior
          << std::endl;

    pcout << "Number of dofs for the skeleton: " << n_dofs_skeleton
          << std::endl;

    pcout << "Number of dofs for the test space: " << n_dofs_test << std::endl;
    pcout << std::endl;

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        error_table.add_value("dofs_interior", n_dofs_interior);
        error_table.add_value("dofs_skeleton", n_dofs_skeleton);
        error_table.add_value("dofs_test", n_dofs_test);
      }

    // Initialize the index sets for the different DoFHandlers
    locally_owned_dofs_trial_interior =
      dof_handler_trial_interior.locally_owned_dofs();
    locally_owned_dofs_trial_skeleton =
      dof_handler_trial_skeleton.locally_owned_dofs();
    locally_owned_dofs_test = dof_handler_test.locally_owned_dofs();

    locally_relevant_dofs_trial_interior =
      DoFTools::extract_locally_relevant_dofs(dof_handler_trial_interior);
    locally_relevant_dofs_trial_skeleton =
      DoFTools::extract_locally_relevant_dofs(dof_handler_trial_skeleton);
    locally_relevant_dofs_test =
      DoFTools::extract_locally_relevant_dofs(dof_handler_test);

    // We apply the boundary conditions using AffineConstraints. The boundary
    // conditions are n x E = 0 on the PEC boundaries.
    constraints.clear();
    constraints.reinit(locally_owned_dofs_trial_skeleton,
                       locally_relevant_dofs_trial_skeleton);
    DoFTools::make_hanging_node_constraints(dof_handler_trial_skeleton,
                                            constraints);

    // -x boundary n x E = 0
    VectorTools::project_boundary_values_curl_conforming_l2(
      dof_handler_trial_skeleton,
      0,
      Functions::ZeroFunction<dim>(4 * dim),
      0,
      constraints);
    VectorTools::project_boundary_values_curl_conforming_l2(
      dof_handler_trial_skeleton,
      dim,
      Functions::ZeroFunction<dim>(4 * dim),
      0,
      constraints);

    // x boundary n x E = 0
    VectorTools::project_boundary_values_curl_conforming_l2(
      dof_handler_trial_skeleton,
      0,
      Functions::ZeroFunction<dim>(4 * dim),
      1,
      constraints);
    VectorTools::project_boundary_values_curl_conforming_l2(
      dof_handler_trial_skeleton,
      dim,
      Functions::ZeroFunction<dim>(4 * dim),
      1,
      constraints);

    // -y boundary n x E = 0
    VectorTools::project_boundary_values_curl_conforming_l2(
      dof_handler_trial_skeleton,
      0,
      Functions::ZeroFunction<dim>(4 * dim),
      2,
      constraints);
    VectorTools::project_boundary_values_curl_conforming_l2(
      dof_handler_trial_skeleton,
      dim,
      Functions::ZeroFunction<dim>(4 * dim),
      2,
      constraints);

    // y boundary n x E = 0
    VectorTools::project_boundary_values_curl_conforming_l2(
      dof_handler_trial_skeleton,
      0,
      Functions::ZeroFunction<dim>(4 * dim),
      3,
      constraints);
    VectorTools::project_boundary_values_curl_conforming_l2(
      dof_handler_trial_skeleton,
      dim,
      Functions::ZeroFunction<dim>(4 * dim),
      3,
      constraints);

    // Because we want to use high-order Nedelec elements for the face, but
    // dealii does not support a FE_FaceNedelec, we hack our way around by using
    // full cell elements, but freezing the interior dofs using
    // constrain_dofs_to_zero. To do so, we create a container for all dof
    // indices and a container for the face dof indices and we loop on all the
    // faces of each cell to find which dofs are on the face and flag them. The
    // remaining dofs are then constrained to zero.

    std::vector<types::global_dof_index> cell_dof_indices(
      fe_trial_skeleton.n_dofs_per_cell());
    std::vector<types::global_dof_index> face_dof_indices(
      fe_trial_skeleton.n_dofs_per_face());

    std::vector<bool> is_dof_on_face(fe_trial_skeleton.n_dofs_per_cell(),
                                     false);

    // Loop on all skeleton dofs and set the interior constraints to zero
    for (const auto &cell : dof_handler_trial_skeleton.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            // Get all dof indices on the cell
            cell->get_dof_indices(cell_dof_indices);

            // Loop on all the faces of the cell
            for (const auto &face : cell->face_iterators())
              {
                face->get_dof_indices(face_dof_indices);

                // Loop on all dofs on the face
                for (const auto &face_dof : face_dof_indices)
                  {
                    // Find the first iterator in the cell dof indices that
                    // matches the face dof (find where the face dof is in the
                    // cell dof indices)
                    const auto it = std::find(cell_dof_indices.begin(),
                                              cell_dof_indices.end(),
                                              face_dof);

                    // If the dof is on the face (find returns the second
                    // iterator if no match is found), set the corresponding
                    // flag to true
                    if (it != cell_dof_indices.end())
                      {
                        is_dof_on_face[std::distance(cell_dof_indices.begin(),
                                                     it)] = true;
                      }
                  }
              }

            // Loop on all dofs on the cell and constrain the interior ones to
            // zero
            for (unsigned int index = 0; index < cell_dof_indices.size();
                 ++index)
              {
                // If the dof is not on a face, then it is an interior dof
                if (!is_dof_on_face[index])
                  {
                    constraints.constrain_dof_to_zero(cell_dof_indices[index]);
                  }
              }
          }
      }

    constraints.close();

    // Here we initialize the relevant vectors to solve the linear system
    // associated with the skeleton unknowns. We also initialize all the
    // parallel vector that will contain the solution (interior and on the
    // skeleton) and the residual estimator.
    DynamicSparsityPattern dsp(locally_relevant_dofs_trial_skeleton);
    DoFTools::make_sparsity_pattern(dof_handler_trial_skeleton,
                                    dsp,
                                    constraints,
                                    false);
    SparsityTools::distribute_sparsity_pattern(
      dsp,
      locally_owned_dofs_trial_skeleton,
      mpi_communicator,
      locally_relevant_dofs_trial_skeleton);

    system_matrix.reinit(locally_owned_dofs_trial_skeleton,
                         locally_owned_dofs_trial_skeleton,
                         dsp,
                         mpi_communicator);

    locally_relevant_solution_skeleton.reinit(
      locally_owned_dofs_trial_skeleton,
      locally_relevant_dofs_trial_skeleton,
      mpi_communicator);
    system_rhs.reinit(locally_owned_dofs_trial_skeleton, mpi_communicator);
    locally_relevant_solution_interior.reinit(
      locally_owned_dofs_trial_interior,
      locally_relevant_dofs_trial_interior,
      mpi_communicator);
    locally_relevant_residual_estimator.reinit(locally_owned_dofs_test,
                                               locally_relevant_dofs_test,
                                               mpi_communicator);
  };

  // Here we assemble the system matrix and right-hand side and apply the Robin
  // boundary conditions during the assembly.
  template <int dim>
  void
  Time_Harmonic_Maxwell_DPG<dim>::assemble_system(const bool solve_interior)
  {
    pcout << "*--- Assembling system ---*" << std::endl;

    TimerOutput::Scope t(computing_timer, "assemble_system");

    // We initialize vectors to store the locally owned solution
    LA::MPI::Vector locally_owned_solution_interior(
      locally_owned_dofs_trial_interior, mpi_communicator);
    LA::MPI::Vector locally_owned_residual_estimator(locally_owned_dofs_test,
                                                     mpi_communicator);

    // We first define quadrature rules and related variables. Since the
    // quadrature formula should be the same for both trial and test FE and that
    // the test space have higher polynomial degree than the others by
    // construction, we use it to define the quadrature formula.
    const QGauss<dim>     quadrature_formula(fe_test.degree + 1);
    const QGauss<dim - 1> face_quadrature_formula(fe_test.degree + 1);
    const unsigned int    n_q_points      = quadrature_formula.size();
    const unsigned int    n_face_q_points = face_quadrature_formula.size();

    // We then create the corresponding FEValues and FEFaceValues objects. Note
    // that only the test space needs gradients because of the ultraweak
    // formulation. Similarly, because everything is on the same triangulation,
    // we only need to update the quadrature points and JxW values in one of the
    // spaces. Here we choose the trial space.
    FEValues<dim> fe_values_trial_interior(fe_trial_interior,
                                           quadrature_formula,
                                           update_values |
                                             update_quadrature_points |
                                             update_JxW_values);

    FEValues<dim> fe_values_test(fe_test,
                                 quadrature_formula,
                                 update_values | update_gradients);

    FEFaceValues<dim> fe_face_values_trial_skeleton(fe_trial_skeleton,
                                                    face_quadrature_formula,
                                                    update_values |
                                                      update_quadrature_points |
                                                      update_normal_vectors |
                                                      update_JxW_values);

    FEFaceValues<dim> fe_face_values_test(fe_test,
                                          face_quadrature_formula,
                                          update_values);

    // We also create all the relevant matrices and vector to build the DPG
    // system. To do so we first need the number of dofs per cell for each of
    // the finite element spaces.
    const unsigned int dofs_per_cell_test = fe_test.n_dofs_per_cell();
    const unsigned int dofs_per_cell_trial_interior =
      fe_trial_interior.n_dofs_per_cell();
    const unsigned int dofs_per_cell_trial_skeleton =
      fe_trial_skeleton.n_dofs_per_cell();

    // To avoid unecessary computations, we will precompute the shape functions
    // of all our spaces once and then use those precomputed values to assemble
    // the local matrices. Consequently, we need to create containers to store
    // these values. Since those are complex-valued functions, we use two
    // different containers for the real and imaginary parts because the
    // conjugate of a complex tensor is not implemented in deal.II.
    std::vector<Tensor<1, dim, std::complex<double>>> F(dofs_per_cell_test);
    std::vector<Tensor<1, dim, std::complex<double>>> F_conj(
      dofs_per_cell_test);
    std::vector<Tensor<1, dim, std::complex<double>>> I(dofs_per_cell_test);
    std::vector<Tensor<1, dim, std::complex<double>>> I_conj(
      dofs_per_cell_test);
    std::vector<Tensor<1, dim, std::complex<double>>> curl_F(
      dofs_per_cell_test);
    std::vector<Tensor<1, dim, std::complex<double>>> curl_F_conj(
      dofs_per_cell_test);
    std::vector<Tensor<1, dim, std::complex<double>>> curl_I(
      dofs_per_cell_test);
    std::vector<Tensor<1, dim, std::complex<double>>> curl_I_conj(
      dofs_per_cell_test);
    std::vector<Tensor<1, dim, std::complex<double>>> F_face(
      dofs_per_cell_test);
    std::vector<Tensor<1, dim, std::complex<double>>> F_face_conj(
      dofs_per_cell_test);
    std::vector<Tensor<1, dim, std::complex<double>>> I_face_conj(
      dofs_per_cell_test);
    std::vector<Tensor<1, dim, std::complex<double>>> n_cross_I_face(
      dofs_per_cell_test);
    std::vector<Tensor<1, dim, std::complex<double>>> n_cross_I_face_conj(
      dofs_per_cell_test);

    std::vector<Tensor<1, dim, std::complex<double>>> E(
      dofs_per_cell_trial_interior);
    std::vector<Tensor<1, dim, std::complex<double>>> H(
      dofs_per_cell_trial_interior);
    std::vector<Tensor<1, dim, std::complex<double>>> E_hat(
      dofs_per_cell_trial_skeleton);
    std::vector<Tensor<1, dim, std::complex<double>>> n_cross_E_hat(
      dofs_per_cell_trial_skeleton);
    std::vector<Tensor<1, dim, std::complex<double>>> n_cross_H_hat(
      dofs_per_cell_trial_skeleton);

    // Also, to avoid multiple if calls to understand where assemble each term
    // in the matrix, we will use containers to store dofs relationships for
    // each of the local matrices.

    // G matrix for the Riesz map and relationships between test functions
    std::vector<std::pair<unsigned int, unsigned int>> G_FF;
    std::vector<std::pair<unsigned int, unsigned int>> G_FI;
    std::vector<std::pair<unsigned int, unsigned int>> G_IF;
    std::vector<std::pair<unsigned int, unsigned int>> G_II;

    // B matrix for the bilinear form and relationships between interior trial
    // and test
    std::vector<std::pair<unsigned int, unsigned int>> B_FE;
    std::vector<std::pair<unsigned int, unsigned int>> B_IE;
    std::vector<std::pair<unsigned int, unsigned int>> B_FH;
    std::vector<std::pair<unsigned int, unsigned int>> B_IH;

    // B_hat matrix for the bilinear form and relationships between skeleton
    // trial and test
    std::vector<std::pair<unsigned int, unsigned int>> B_hat_FE;
    std::vector<std::pair<unsigned int, unsigned int>> B_hat_IE;
    std::vector<std::pair<unsigned int, unsigned int>> B_hat_FH;

    // l vector for the linear form and relationships between test functions
    std::vector<unsigned int> l_F;

    // Reserve memory to avoid reallocations of each matrix
    G_FF.reserve(dofs_per_cell_test * dofs_per_cell_test);
    G_FI.reserve(dofs_per_cell_test * dofs_per_cell_test);
    G_IF.reserve(dofs_per_cell_test * dofs_per_cell_test);
    G_II.reserve(dofs_per_cell_test * dofs_per_cell_test);

    B_FE.reserve(dofs_per_cell_trial_interior * dofs_per_cell_test);
    B_IE.reserve(dofs_per_cell_trial_interior * dofs_per_cell_test);
    B_FH.reserve(dofs_per_cell_trial_interior * dofs_per_cell_test);
    B_IH.reserve(dofs_per_cell_trial_interior * dofs_per_cell_test);

    B_hat_FE.reserve(dofs_per_cell_trial_skeleton * dofs_per_cell_test);
    B_hat_IE.reserve(dofs_per_cell_trial_skeleton * dofs_per_cell_test);
    B_hat_FH.reserve(dofs_per_cell_trial_skeleton * dofs_per_cell_test);

    l_F.reserve(dofs_per_cell_test);

    // Here we create the DPG local matrices and vector used for the assembly
    // before condensation.
    LAPACKFullMatrix<double> G_matrix(dofs_per_cell_test, dofs_per_cell_test);

    LAPACKFullMatrix<double> B_matrix(dofs_per_cell_test,
                                      dofs_per_cell_trial_interior);

    LAPACKFullMatrix<double> B_hat_matrix(dofs_per_cell_test,
                                          dofs_per_cell_trial_skeleton);

    Vector<double> l_vector(dofs_per_cell_test);

    // We create the condensation matrices.
    LAPACKFullMatrix<double> M1_matrix(dofs_per_cell_trial_interior,
                                       dofs_per_cell_trial_interior);
    LAPACKFullMatrix<double> M2_matrix(dofs_per_cell_trial_interior,
                                       dofs_per_cell_trial_skeleton);
    LAPACKFullMatrix<double> M3_matrix(dofs_per_cell_trial_skeleton,
                                       dofs_per_cell_trial_skeleton);
    LAPACKFullMatrix<double> M4_matrix(dofs_per_cell_trial_interior,
                                       dofs_per_cell_test);
    LAPACKFullMatrix<double> M5_matrix(dofs_per_cell_trial_skeleton,
                                       dofs_per_cell_test);

    // During the calculation of matrix vector product, we require intermediary
    // matrices and vector that we also allocate here.
    LAPACKFullMatrix<double> tmp_matrix(dofs_per_cell_trial_skeleton,
                                        dofs_per_cell_trial_interior);

    LAPACKFullMatrix<double> tmp_matrix2(dofs_per_cell_trial_skeleton,
                                         dofs_per_cell_trial_skeleton);

    LAPACKFullMatrix<double> tmp_matrix3(dofs_per_cell_trial_skeleton,
                                         dofs_per_cell_test);

    Vector<double> tmp_vector(dofs_per_cell_trial_interior);
    Vector<double> tmp_vector2(dofs_per_cell_test);

    // We create the cell matrix and the RHS that will be distributed in the
    // full system after the assembly along with the index’s mapping.
    FullMatrix<double> cell_matrix(dofs_per_cell_trial_skeleton,
                                   dofs_per_cell_trial_skeleton);
    Vector<double>     cell_skeleton_rhs(dofs_per_cell_trial_skeleton);

    std::vector<types::global_dof_index> local_dof_indices(
      dofs_per_cell_trial_skeleton);

    // Finally, when reconstructing the interior solution from the skeleton, we
    // require additional vectors that we allocate here.
    Vector<double> cell_interior_rhs(dofs_per_cell_trial_interior);
    Vector<double> cell_interior_solution(dofs_per_cell_trial_interior);
    Vector<double> cell_skeleton_solution(dofs_per_cell_trial_skeleton);
    Vector<double> cell_residual(dofs_per_cell_test);

    // We define some constants that will be used during the assembly. Those
    // would change according to the material parameters, but here we only have
    // one material.
    const std::complex<double> ikZr      = imag * k * Z_r;
    const std::complex<double> conj_ikZr = conj(ikZr);
    const std::complex<double> ikepsilon_Zr =
      imag * k / Z_r * (1. - imag * parameters.sigma_r / parameters.epsilon_r);
    const std::complex<double> conj_ikepsilon_Zr = conj(ikepsilon_Zr);
    const std::complex<double> kz_kZr            = k_z / (k * Z_r);
    const std::complex<double> conj_kz_kZr       = conj(kz_kZr);

    // We also create objects to store the incident fields that will be used to
    // compute the excitation for the robin boundary condition.
    Tensor<1, dim, std::complex<double>> H_inc;
    Tensor<1, dim, std::complex<double>> E_inc;
    Tensor<1, dim, std::complex<double>> g_inc;
    AnalyticalSolution_E<dim>            analytical_solution_E(parameters);
    AnalyticalSolution_H<dim>            analytical_solution_H(parameters);

    // As it is standard we first loop over the cells of the triangulation. Here
    // we have the choice of the DoFHandler to perform this loop. We use the
    // DoFHandler associated with the trial space.
    for (const auto &cell : dof_handler_trial_interior.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            // We reinitialize the FEValues objects to the current cell.
            fe_values_trial_interior.reinit(cell);

            // We will also need to reinitialize the FEValues for the test
            // space and make sure that is the same cell as the one used for the
            // trial space.
            const typename DoFHandler<dim>::active_cell_iterator cell_test =
              cell->as_dof_handler_iterator(dof_handler_test);
            fe_values_test.reinit(cell_test);

            // Similarly, we reinitialize the FEValues for the trial space on
            // the skeleton, but this will not be used before we also loop on
            // the cells faces.
            const typename DoFHandler<dim>::active_cell_iterator cell_skeleton =
              cell->as_dof_handler_iterator(dof_handler_trial_skeleton);

            // We then reinitialize all the matrices that we are aggregating
            // information for each cell.
            G_matrix     = 0;
            B_matrix     = 0;
            B_hat_matrix = 0;
            l_vector     = 0;

            // We also need to reinitialize the $M_1$ condensation matrix
            // between each iteration on cell to get rid of its inverse status.
            M1_matrix = 0;

            // Reset the dofs relationships containers
            G_FF.clear();
            G_FI.clear();
            G_IF.clear();
            G_II.clear();

            B_FE.clear();
            B_IE.clear();
            B_FH.clear();
            B_IH.clear();

            l_F.clear();

            // We fill the dofs relationship containers at the cell level. To do
            // so, we first loop on the test space dofs.
            for (unsigned int i : fe_values_test.dof_indices())
              {
                // Get the information on which element the dof is
                const unsigned int current_element_test_i =
                  fe_test.system_to_base_index(i).first.first;

                // Fill the l vector relationship
                if ((current_element_test_i == 0) ||
                    (current_element_test_i == 1))
                  {
                    l_F.emplace_back(i);
                  }

                // Loop over the dofs test a second time to fill the dofs
                // relationship for the G matrix (Riesz map)
                for (unsigned int j : fe_values_test.dof_indices())
                  {
                    const unsigned int current_element_test_j =
                      fe_test.system_to_base_index(j).first.first;
                    if (((current_element_test_i == 0) ||
                         (current_element_test_i == 1)) &&
                        ((current_element_test_j == 0) ||
                         (current_element_test_j == 1)))
                      {
                        G_FF.emplace_back(i, j);
                      }
                    if (((current_element_test_i == 0) ||
                         (current_element_test_i == 1)) &&
                        ((current_element_test_j == 2) ||
                         (current_element_test_j == 3)))
                      {
                        G_FI.emplace_back(i, j);
                      }
                    if (((current_element_test_i == 2) ||
                         (current_element_test_i == 3)) &&
                        ((current_element_test_j == 0) ||
                         (current_element_test_j == 1)))
                      {
                        G_IF.emplace_back(i, j);
                      }
                    if (((current_element_test_i == 2) ||
                         (current_element_test_i == 3)) &&
                        ((current_element_test_j == 2) ||
                         (current_element_test_j == 3)))
                      {
                        G_II.emplace_back(i, j);
                      }
                  }

                // Then we loop over the dofs trial space to fill the dofs
                // relationship for the B matrix (bilinear form)
                for (unsigned int j : fe_values_trial_interior.dof_indices())
                  {
                    const unsigned int current_element_trial_j =
                      fe_trial_interior.system_to_base_index(j).first.first;

                    if (((current_element_test_i == 0) ||
                         (current_element_test_i == 1)) &&
                        ((current_element_trial_j == 0) ||
                         (current_element_trial_j == 1)))
                      {
                        B_FE.emplace_back(i, j);
                      }
                    if (((current_element_test_i == 0) ||
                         (current_element_test_i == 1)) &&
                        ((current_element_trial_j == 2) ||
                         (current_element_trial_j == 3)))
                      {
                        B_FH.emplace_back(i, j);
                      }
                    if (((current_element_test_i == 2) ||
                         (current_element_test_i == 3)) &&
                        ((current_element_trial_j == 0) ||
                         (current_element_trial_j == 1)))
                      {
                        B_IE.emplace_back(i, j);
                      }
                    if (((current_element_test_i == 2) ||
                         (current_element_test_i == 3)) &&
                        ((current_element_trial_j == 2) ||
                         (current_element_trial_j == 3)))
                      {
                        B_IH.emplace_back(i, j);
                      }
                  }
              }

            // Now we loop over all quadrature points of the cell
            for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
              {
                // To avoid unnecessary computation, we fill the shape values
                // containers for the real and imaginary parts of the electric
                // and magnetic fields and the Dof relationship at the current
                // quadrature point.
                const double &JxW = fe_values_trial_interior.JxW(q_point);

                for (unsigned int i : fe_values_test.dof_indices())
                  {
                    F[i] =
                      fe_values_test[extractor_E_real].value(i, q_point) +
                      imag * fe_values_test[extractor_E_imag].value(i, q_point);
                    F_conj[i] =
                      fe_values_test[extractor_E_real].value(i, q_point) -
                      imag * fe_values_test[extractor_E_imag].value(i, q_point);

                    curl_F[i] =
                      fe_values_test[extractor_E_real].curl(i, q_point) +
                      imag * fe_values_test[extractor_E_imag].curl(i, q_point);
                    curl_F_conj[i] =
                      fe_values_test[extractor_E_real].curl(i, q_point) -
                      imag * fe_values_test[extractor_E_imag].curl(i, q_point);

                    I[i] =
                      fe_values_test[extractor_H_real].value(i, q_point) +
                      imag * fe_values_test[extractor_H_imag].value(i, q_point);
                    I_conj[i] =
                      fe_values_test[extractor_H_real].value(i, q_point) -
                      imag * fe_values_test[extractor_H_imag].value(i, q_point);

                    curl_I[i] =
                      fe_values_test[extractor_H_real].curl(i, q_point) +
                      imag * fe_values_test[extractor_H_imag].curl(i, q_point);
                    curl_I_conj[i] =
                      fe_values_test[extractor_H_real].curl(i, q_point) -
                      imag * fe_values_test[extractor_H_imag].curl(i, q_point);
                  }

                for (unsigned int i : fe_values_trial_interior.dof_indices())
                  {
                    E[i] =
                      fe_values_trial_interior[extractor_E_real].value(
                        i, q_point) +
                      imag * fe_values_trial_interior[extractor_E_imag].value(
                               i, q_point);
                    H[i] =
                      fe_values_trial_interior[extractor_H_real].value(
                        i, q_point) +
                      imag * fe_values_trial_interior[extractor_H_imag].value(
                               i, q_point);
                  }

                // Now we loop on each relationship container to assemble the
                // relevant matrices
                for (const auto &[i, j] : G_FF)
                  {
                    G_matrix(i, j) +=
                      (((F[j] * F_conj[i]) + (curl_F[j] * curl_F_conj[i]) +
                        (conj_ikepsilon_Zr * F[j] * ikepsilon_Zr * F_conj[i])) *
                       JxW)
                        .real();
                  }

                for (const auto &[i, j] : G_FI)
                  {
                    G_matrix(i, j) += (((curl_I[j] * ikepsilon_Zr * F_conj[i]) -
                                        (conj_ikZr * I[j] * curl_F_conj[i])) *
                                       JxW)
                                        .real();
                  }

                for (const auto &[i, j] : G_IF)
                  {
                    G_matrix(i, j) +=
                      (((conj_ikepsilon_Zr * F[j] * curl_I_conj[i]) -
                        (curl_F[j] * ikZr * I_conj[i])) *
                       JxW)
                        .real();
                  }

                for (const auto &[i, j] : G_II)
                  {
                    G_matrix(i, j) +=
                      (((I[j] * I_conj[i]) + (curl_I[j] * curl_I_conj[i]) +
                        (conj_ikZr * I[j] * ikZr * I_conj[i])) *
                       JxW)
                        .real();
                  }

                for (const auto &[i, j] : B_FE)
                  {
                    B_matrix(i, j) +=
                      (ikepsilon_Zr * E[j] * F_conj[i] * JxW).real();
                  }

                for (const auto &[i, j] : B_FH)
                  {
                    B_matrix(i, j) += (H[j] * curl_F_conj[i] * JxW).real();
                  }

                for (const auto &[i, j] : B_IE)
                  {
                    B_matrix(i, j) += (E[j] * curl_I_conj[i] * JxW).real();
                  }

                for (const auto &[i, j] : B_IH)
                  {
                    B_matrix(i, j) -= (ikZr * H[j] * I_conj[i] * JxW).real();
                  }

                for (const auto &i : l_F)
                  {
                    l_vector[i] += 0.0;
                  }
              }

            // We now build the skeleton terms. Similarly, we choose to loop on
            // the skeleton trial space faces.
            for (const auto &face : cell_skeleton->face_iterators())
              {
                // We reinitialize the FEFaceValues objects to the current
                // faces.
                fe_face_values_test.reinit(cell_test, face);
                fe_face_values_trial_skeleton.reinit(cell_skeleton, face);

                // Reset the face dofs relationships
                G_FF.clear();
                G_FI.clear();
                G_IF.clear();
                G_II.clear();

                B_hat_FH.clear();
                B_hat_IE.clear();
                B_hat_FE.clear();

                l_F.clear();

                // We fill the dofs relationship containers at the face level.
                // To do so, we first loop on the test space dofs.
                for (unsigned int i : fe_face_values_test.dof_indices())
                  {
                    // Get the information on which element the dof is
                    const unsigned int current_element_test_i =
                      fe_test.system_to_base_index(i).first.first;

                    // If we are at either boundary 4 or 5, we are on the Robin
                    // B.C. and the load vector and the B_hat matrix will have a
                    // contribution in addition to a modification of the riesz
                    // map (G matrix) because of the energy norm that we want to
                    // minimize there.
                    if ((face->boundary_id() == 4) ||
                        (face->boundary_id() == 5))
                      {
                        if ((current_element_test_i == 0) ||
                            (current_element_test_i == 1))
                          {
                            l_F.emplace_back(i);
                          }

                        // Loop over the dofs test to fill the G_matrix dofs
                        // relationship
                        for (unsigned int j : fe_face_values_test.dof_indices())
                          {
                            const unsigned int current_element_test_j =
                              fe_test.system_to_base_index(j).first.first;

                            if (((current_element_test_i == 0) ||
                                 (current_element_test_i == 1)) &&
                                ((current_element_test_j == 0) ||
                                 (current_element_test_j == 1)))
                              {
                                G_FF.emplace_back(i, j);
                              }
                            if (((current_element_test_i == 0) ||
                                 (current_element_test_i == 1)) &&
                                ((current_element_test_j == 2) ||
                                 (current_element_test_j == 3)))
                              {
                                G_FI.emplace_back(i, j);
                              }
                            if (((current_element_test_i == 2) ||
                                 (current_element_test_i == 3)) &&
                                ((current_element_test_j == 0) ||
                                 (current_element_test_j == 1)))
                              {
                                G_IF.emplace_back(i, j);
                              }
                            if (((current_element_test_i == 2) ||
                                 (current_element_test_i == 3)) &&
                                ((current_element_test_j == 2) ||
                                 (current_element_test_j == 3)))
                              {
                                G_II.emplace_back(i, j);
                              }
                          }
                        // Loop over the dofs trial space to fill the B_hat
                        // matrix dofs relationship for the Robin boundary
                        // condition
                        for (unsigned int j :
                             fe_face_values_trial_skeleton.dof_indices())
                          {
                            const unsigned int current_element_trial_j =
                              fe_trial_skeleton.system_to_base_index(j)
                                .first.first;

                            if (((current_element_test_i == 0) ||
                                 (current_element_test_i == 1)) &&
                                ((current_element_trial_j == 0) ||
                                 (current_element_trial_j == 1)))
                              {
                                B_hat_FE.emplace_back(i, j);
                              }
                            if (((current_element_test_i == 2) ||
                                 (current_element_test_i == 3)) &&
                                ((current_element_trial_j == 0) ||
                                 (current_element_trial_j == 1)))
                              {
                                B_hat_IE.emplace_back(i, j);
                              }
                          }
                      }
                    else
                      {
                        // If not on Robin B.C., assemble all the other
                        // relevantskeleton terms
                        for (unsigned int j :
                             fe_face_values_trial_skeleton.dof_indices())
                          {
                            const unsigned int current_element_trial_j =
                              fe_trial_skeleton.system_to_base_index(j)
                                .first.first;

                            if (((current_element_test_i == 0) ||
                                 (current_element_test_i == 1)) &&
                                ((current_element_trial_j == 2) ||
                                 (current_element_trial_j == 3)))
                              {
                                B_hat_FH.emplace_back(i, j);
                              }
                            if (((current_element_test_i == 2) ||
                                 (current_element_test_i == 3)) &&
                                ((current_element_trial_j == 0) ||
                                 (current_element_trial_j == 1)))
                              {
                                B_hat_IE.emplace_back(i, j);
                              }
                          }
                      }
                  }

                // Loop over all face quadrature points
                for (unsigned int q_point = 0; q_point < n_face_q_points;
                     ++q_point)
                  {
                    // Initialize reusable variables
                    const auto &position =
                      fe_face_values_trial_skeleton.quadrature_point(q_point);
                    const auto &normal =
                      fe_face_values_trial_skeleton.normal_vector(q_point);
                    const double JxW_face =
                      fe_face_values_trial_skeleton.JxW(q_point);

                    // As for the cell, we first loop over the test dofs to fill
                    // the face values containers
                    for (unsigned int i : fe_face_values_test.dof_indices())
                      {
                        F_face[i] =
                          fe_face_values_test[extractor_E_real].value(i,
                                                                      q_point) +
                          imag * fe_face_values_test[extractor_E_imag].value(
                                   i, q_point);
                        F_face_conj[i] =
                          fe_face_values_test[extractor_E_real].value(i,
                                                                      q_point) -
                          imag * fe_face_values_test[extractor_E_imag].value(
                                   i, q_point);

                        I_face_conj[i] =
                          fe_face_values_test[extractor_H_real].value(i,
                                                                      q_point) -
                          imag * fe_face_values_test[extractor_H_imag].value(
                                   i, q_point);

                        n_cross_I_face[i] = cross_product_3d(
                          normal,
                          fe_face_values_test[extractor_H_real].value(i,
                                                                      q_point) +
                            imag * fe_face_values_test[extractor_H_imag].value(
                                     i, q_point));
                        n_cross_I_face_conj[i] = cross_product_3d(
                          normal,
                          fe_face_values_test[extractor_H_real].value(i,
                                                                      q_point) -
                            imag * fe_face_values_test[extractor_H_imag].value(
                                     i, q_point));
                      }

                    // Then, similarly we loop over the trial dofs to fill the
                    // face values containers. Note that the to be in
                    // H^-1/2(curl), the fields needs to have the tangential
                    // property mapping (n x (E x n)) which effectively extract
                    // the tangential component of the field at the face. So
                    // here we apply this operation using the map_H12 function
                    // that we defined earlier. Stricly speeking, nx(E_parallel)
                    // = n x E, and we would not need to use the map_H12
                    // function, but we keep it for consistency.
                    for (unsigned int i :
                         fe_face_values_trial_skeleton.dof_indices())
                      {
                        E_hat[i] = map_H12(
                          fe_face_values_trial_skeleton[extractor_E_real].value(
                            i, q_point) +
                            imag *
                              fe_face_values_trial_skeleton[extractor_E_imag]
                                .value(i, q_point),
                          normal);

                        n_cross_E_hat[i] = cross_product_3d(
                          normal,
                          map_H12(
                            fe_face_values_trial_skeleton[extractor_E_real]
                                .value(i, q_point) +
                              imag *
                                fe_face_values_trial_skeleton[extractor_E_imag]
                                  .value(i, q_point),
                            normal));

                        n_cross_H_hat[i] = cross_product_3d(
                          normal,
                          map_H12(
                            fe_face_values_trial_skeleton[extractor_H_real]
                                .value(i, q_point) +
                              imag *
                                fe_face_values_trial_skeleton[extractor_H_imag]
                                  .value(i, q_point),
                            normal));
                      }

                    // Here we get the excitation at the inlet port. To do so we
                    // force the electromagnetic fields with our analytical
                    // solution at the boundary.
                    if (face->boundary_id() == 4)
                      {
                        H_inc[0] =
                          analytical_solution_H.value(position, 0) +
                          imag * analytical_solution_H.value(position, 3);
                        H_inc[1] =
                          analytical_solution_H.value(position, 1) +
                          imag * analytical_solution_H.value(position, 4);
                        H_inc[2] =
                          analytical_solution_H.value(position, 2) +
                          imag * analytical_solution_H.value(position, 5);

                        E_inc[0] =
                          analytical_solution_E.value(position, 0) +
                          imag * analytical_solution_E.value(position, 3);
                        E_inc[1] =
                          analytical_solution_E.value(position, 1) +
                          imag * analytical_solution_E.value(position, 4);
                        E_inc[2] =
                          analytical_solution_E.value(position, 2) +
                          imag * analytical_solution_E.value(position, 5);

                        g_inc = cross_product_3d(normal, H_inc) +
                                map_H12(kz_kZr * E_inc, normal);
                      }
                    // If we are at the outlet, we want to have an absorbing
                    // boundary condition so we set the excitation to zero.
                    if (face->boundary_id() == 5)
                      {
                        g_inc = 0.0;
                      }

                    // Now we loop on each relationship container to assemble
                    // the relevant matrices.
                    for (const auto &[i, j] : G_FF)
                      {
                        G_matrix(i, j) += (conj_kz_kZr * F_face[j] * kz_kZr *
                                           F_face_conj[i] * JxW_face)
                                            .real();
                      }

                    for (const auto &[i, j] : G_FI)
                      {
                        G_matrix(i, j) += (n_cross_I_face[j] * kz_kZr *
                                           F_face_conj[i] * JxW_face)
                                            .real();
                      }

                    for (const auto &[i, j] : G_IF)
                      {
                        G_matrix(i, j) += (conj_kz_kZr * F_face[j] *
                                           n_cross_I_face_conj[i] * JxW_face)
                                            .real();
                      }

                    for (const auto &[i, j] : G_II)
                      {
                        G_matrix(i, j) += (n_cross_I_face[j] *
                                           n_cross_I_face_conj[i] * JxW_face)
                                            .real();
                      }

                    for (const auto &[i, j] : B_hat_FH)
                      {
                        B_hat_matrix(i, j) +=
                          (n_cross_H_hat[j] * F_face_conj[i] * JxW_face).real();
                      }

                    for (const auto &[i, j] : B_hat_IE)
                      {
                        B_hat_matrix(i, j) +=
                          (n_cross_E_hat[j] * I_face_conj[i] * JxW_face).real();
                      }

                    for (const auto &[i, j] : B_hat_FE)
                      {
                        B_hat_matrix(i, j) -=
                          (kz_kZr * E_hat[j] * F_face_conj[i] * JxW_face)
                            .real();
                      }

                    for (const auto &i : l_F)
                      {
                        l_vector[i] -=
                          (g_inc * F_face_conj[i] * JxW_face).real();
                      }
                  }
              } // End of face loop

            // Finally, after having assembled all the matrices and vectors, we
            // build the condensed version of the system.

            // We only need the inverse of the Gram matrix $G$, so we
            // invert it.
            G_matrix.invert();

            // We construct $M_4 = B^\dagger G^{-1}$ and $M_5 = \hat{B}^\dagger
            // G^{-1}$ with it:
            B_matrix.Tmmult(M4_matrix, G_matrix);
            B_hat_matrix.Tmmult(M5_matrix, G_matrix);

            // Then using $M_4$ we compute the condensed matrix $M_1 = B^\dagger
            // G^{-1} B$ and $M_2 = B^\dagger G^{-1} \hat{B}$:
            M4_matrix.mmult(M1_matrix, B_matrix);
            M4_matrix.mmult(M2_matrix, B_hat_matrix);

            // We also compute the matrix $M_3 = \hat{B}^\dagger G^{-1} \hat{B}$
            M5_matrix.mmult(M3_matrix, B_hat_matrix);

            // Finally, as for the $G$ matrix, we invert the $M_1$
            // matrix:
            M1_matrix.invert();

            // If the flag solves interior is set to true, we have already the
            // solution on the skeleton and only need to perform $u_h = M_1^{-1}
            // (M_4 l - M_2 \hat{u}_h)$ on each cell. When this is obtained, we
            // can perform at the same time the residual estimator (\Psi =
            // G^{-1}(l-B u_h
            // - \hat{B}\hat{u}_h)).
            if (solve_interior)
              {
                // We first get the solution vector for this cell.
                cell_skeleton->get_dof_values(
                  locally_relevant_solution_skeleton, cell_skeleton_solution);

                // Then we do the matrix-vector product to obtain the interior
                // unknowns.
                M2_matrix.vmult(tmp_vector, cell_skeleton_solution);
                M4_matrix.vmult(cell_interior_rhs, l_vector);
                cell_interior_rhs -= tmp_vector;
                M1_matrix.vmult(cell_interior_solution, cell_interior_rhs);

                // Finally, we map the cell interior solution to the global
                // interior solution.
                cell->distribute_local_to_global(
                  cell_interior_solution, locally_owned_solution_interior);

                // Then we compute the residual estimator and also map it to the
                // global residual estimator vector.
                B_matrix.vmult(tmp_vector2, cell_interior_solution);
                B_hat_matrix.vmult_add(tmp_vector2, cell_skeleton_solution);
                l_vector -= tmp_vector2;
                G_matrix.vmult(cell_residual, l_vector);
                cell_test->distribute_local_to_global(
                  cell_residual, locally_owned_residual_estimator);
              }
            // If the flag solves interior is set to false, we have to compute
            // the local matrix and the local RHS for the condensed system.
            else
              {
                // So the cell matrix is obtained with the formula $(M_3 -
                // M_2^\dagger M_1^{-1} M_2)$:
                M2_matrix.Tmmult(tmp_matrix, M1_matrix);
                tmp_matrix.mmult(tmp_matrix2, M2_matrix);
                tmp_matrix2.add(-1.0, M3_matrix);
                tmp_matrix2 *= -1.0;
                // This line is used to convert the LAPACK matrix to a full
                // matrix so we can perform the distribution to the global
                // system below.
                cell_matrix = tmp_matrix2;

                // Then we compute the cell RHS using $(M_5 -
                // M_2^\dagger M_1^{-1} M_4)l -
                // G$.
                tmp_matrix.mmult(tmp_matrix3, M4_matrix);
                M5_matrix.add(-1.0, tmp_matrix3);
                M5_matrix.vmult(cell_skeleton_rhs, l_vector);

                // Map to global matrix
                cell_skeleton->get_dof_indices(local_dof_indices);
                constraints.distribute_local_to_global(cell_matrix,
                                                       cell_skeleton_rhs,
                                                       local_dof_indices,
                                                       system_matrix,
                                                       system_rhs);
              }
          }
      }

    // After the loop over the cells, we finalize the assembly by compressing
    // the vectors because of the MPI parallelization according to the flag
    // solve_interior.
    if (solve_interior)
      {
        locally_owned_solution_interior.compress(VectorOperation::insert);
        locally_owned_residual_estimator.compress(VectorOperation::insert);

        locally_relevant_solution_interior  = locally_owned_solution_interior;
        locally_relevant_residual_estimator = locally_owned_residual_estimator;
      }
    else
      {
        system_matrix.compress(VectorOperation::add);
        system_rhs.compress(VectorOperation::add);
      }
  }

  // The solve function is in charge of solving the linear system assembled and
  // has nothing specific to DPG per se. Nonetheless, the method allows us to
  // use the Conjugate Gradient iterative solver. Note that because we do not
  // have any preconditioner, the number of iterations can be quite high. For
  // simplicity, we put a high upper limit on the number of iterations, but in
  // practice one would want to change this function to have a more robust
  // solver. The tolerance for the convergence here is defined proportional to
  // the $L^2$ norm of the RHS vector so the stopping criterion is scaled aware.
  template <int dim>
  void
  Time_Harmonic_Maxwell_DPG<dim>::solve_skeleton()
  {
    pcout << "*--- Solving system ---*" << std::endl;

    TimerOutput::Scope t(computing_timer, "solve");

    LA::MPI::Vector completely_distributed_solution(
      locally_owned_dofs_trial_skeleton, mpi_communicator);
    SolverControl solver_control(dof_handler_trial_skeleton.n_dofs(),
                                 1e-10 * system_rhs.l2_norm());
    LA::SolverCG  solver(solver_control);
    TrilinosWrappers::PreconditionIdentity preconditioner;
    solver.solve(system_matrix,
                 completely_distributed_solution,
                 system_rhs,
                 preconditioner);

    constraints.distribute(completely_distributed_solution);
    completely_distributed_solution.update_ghost_values();
    locally_relevant_solution_skeleton = completely_distributed_solution;
    locally_relevant_solution_skeleton.update_ghost_values();

    unsigned int n_CG_iter = solver_control.last_step();
    pcout << "   " << n_CG_iter
          << " CG iterations needed to obtain convergence." << std::endl;
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        error_table.add_value("n_iter", n_CG_iter);
      }
  }

  // This function is standard and is in charge of outputting the
  // solution to a file that can be visualized with Paraview or Visit (VTU and
  // PVTU format). However, because the skeleton solution lives only on the mesh
  // faces, we made use of the DataOutFaces class to output it properly. Note
  // that the subdomains are only outputed in the interior solution file.
  template <int dim>
  void
  Time_Harmonic_Maxwell_DPG<dim>::output_results(const unsigned int cycle)
  {
    pcout << "*--- Outputting results ---*" << std::endl;

    TimerOutput::Scope t(computing_timer, "output_results");

    // Organize the solution output
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler_trial_interior);

    std::vector<std::string> solution_interior_names(dim, "E_real");
    for (unsigned int i = 0; i < dim; ++i)
      {
        solution_interior_names.emplace_back("E_imag");
      }
    for (unsigned int i = 0; i < dim; ++i)
      {
        solution_interior_names.emplace_back("H_real");
      }
    for (unsigned int i = 0; i < dim; ++i)
      {
        solution_interior_names.emplace_back("H_imag");
      }

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
        4 * dim, DataComponentInterpretation::component_is_part_of_vector);

    // Output the solution
    data_out.add_data_vector(locally_relevant_solution_interior,
                             solution_interior_names,
                             DataOut<dim>::type_automatic,
                             data_component_interpretation);

    // Out the subdomains
    Vector<float> subdomain(triangulation.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
    data_out.add_data_vector(subdomain, "subdomain");

    data_out.build_patches(fe_trial_interior.degree);
    data_out.write_vtu_with_pvtu_record(
      "./output/", "solution_waveguide", cycle, mpi_communicator, 3);

    // Organize the error estimator output
    DataOut<dim> data_out_test;
    data_out_test.attach_dof_handler(dof_handler_test);
    data_out_test.attach_triangulation(triangulation);

    std::vector<std::string> residuals_names(dim, "F_real");
    for (unsigned int i = 0; i < dim; ++i)
      {
        residuals_names.emplace_back("F_imag");
      }
    for (unsigned int i = 0; i < dim; ++i)
      {
        residuals_names.emplace_back("I_real");
      }
    for (unsigned int i = 0; i < dim; ++i)
      {
        residuals_names.emplace_back("I_imag");
      }

    // Output the residual estimator and the estimated error per cell
    data_out_test.add_data_vector(locally_relevant_residual_estimator,
                                  residuals_names,
                                  DataOut<dim>::type_automatic,
                                  data_component_interpretation);

    data_out_test.add_data_vector(
      estimated_error_per_cell,
      std::string("error_norm"),
      DataOut<dim>::type_automatic,
      std::vector<DataComponentInterpretation::DataComponentInterpretation>(
        1, DataComponentInterpretation::component_is_scalar));

    data_out_test.build_patches(fe_test.degree);
    data_out_test.write_vtu_with_pvtu_record(
      "./output/", "residual_waveguide", cycle, mpi_communicator, 3);

    // Organize the skeleton solution output
    DataOutFaces<dim> data_out_faces(false);
    data_out_faces.attach_dof_handler(dof_handler_trial_skeleton);

    std::vector<std::string> solution_skeleton_names(dim, "E_hat_real");
    for (unsigned int i = 0; i < dim; ++i)
      {
        solution_skeleton_names.emplace_back("E_hat_imag");
      }
    for (unsigned int i = 0; i < dim; ++i)
      {
        solution_skeleton_names.emplace_back("H_hat_real");
      }
    for (unsigned int i = 0; i < dim; ++i)
      {
        solution_skeleton_names.emplace_back("H_hat_imag");
      }

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation_skeleton(
        4 * dim, DataComponentInterpretation::component_is_part_of_vector);

    // Output the skeleton solution
    data_out_faces.add_data_vector(locally_relevant_solution_skeleton,
                                   solution_skeleton_names,
                                   DataOutFaces<dim>::type_automatic,
                                   data_component_interpretation_skeleton);

    data_out_faces.build_patches(fe_trial_skeleton.degree);
    data_out_faces.write_vtu_with_pvtu_record(
      "./output/", "solution-face_waveguide", cycle, mpi_communicator, 3);
  }

  // In the function calculate_L2_error, we compute the $L^2$ error of each
  // component of our solution, i.e., for the real and imaginary parts of the
  // electric and magnetic fields for both the interior and skeleton solutions.
  // Because we want to have all the error for all the different components of
  // our solution vectors separately, we cannot use the
  // VectorTools::integrate_difference function directly. Instead, we will
  // perform the computation "by hand" by looping over all the cells and faces,
  // interpolating the solution at the quadrature points, and computing the
  // error with respect to the analytical solution at those points. Since we are
  // performing this loop, at the same time we can also compute the error on
  // each cell using the energy norm that we defined for the problem. This will
  // be useful for the adaptive mesh refinement later on.
  template <int dim>
  void
  Time_Harmonic_Maxwell_DPG<dim>::calculate_L2_error()
  {
    pcout << "*--- Calculating L2 error ---*" << std::endl;

    // We first need to instantiate the FEValues and FEFaceValues objects as it
    // has been done during the assembly to manage the solution when we loop
    // over the cells and faces.
    QGauss<dim>           quadrature_formula(fe_test.degree + 1);
    FEValues<dim>         fe_values_trial_interior(fe_trial_interior,
                                           quadrature_formula,
                                           update_values |
                                             update_quadrature_points |
                                             update_JxW_values);
    const QGauss<dim - 1> face_quadrature_formula(fe_test.degree + 1);
    FEFaceValues<dim>     fe_face_values_trial_skeleton(fe_trial_skeleton,
                                                    face_quadrature_formula,
                                                    update_values |
                                                      update_quadrature_points |
                                                      update_normal_vectors |
                                                      update_JxW_values);

    FEValues<dim>     fe_values_test(fe_test,
                                 quadrature_formula,
                                 update_values | update_gradients |
                                   update_quadrature_points);
    FEFaceValues<dim> fe_face_values_test(fe_test,
                                          face_quadrature_formula,
                                          update_values |
                                            update_quadrature_points);

    const unsigned int n_q_points      = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    // We create variables that will store all the integration result we are
    // interested in.
    double L2_error_E_real     = 0;
    double L2_error_E_imag     = 0;
    double L2_error_E_hat_real = 0;
    double L2_error_E_hat_imag = 0;
    double L2_error_H_real     = 0;
    double L2_error_H_imag     = 0;
    double L2_error_H_hat_real = 0;
    double L2_error_H_hat_imag = 0;

    // We create variables to store the analytical solutions
    Tensor<1, dim, std::complex<double>> analytical_E;
    Tensor<1, dim, std::complex<double>> analytical_H;

    // When looping on each cell or face, we will extract the different field
    // solution obtain numerically. The containers used to store the
    // interpolated solution at the quadrature points are declared below.
    std::vector<Tensor<1, dim>> local_E_values_real(n_q_points);
    std::vector<Tensor<1, dim>> local_E_values_imag(n_q_points);
    std::vector<Tensor<1, dim>> local_H_values_real(n_q_points);
    std::vector<Tensor<1, dim>> local_H_values_imag(n_q_points);
    std::vector<Tensor<1, dim>> local_E_hat_values_real(n_face_q_points);
    std::vector<Tensor<1, dim>> local_E_hat_values_imag(n_face_q_points);
    std::vector<Tensor<1, dim>> local_H_hat_values_real(n_face_q_points);
    std::vector<Tensor<1, dim>> local_H_hat_values_imag(n_face_q_points);

    // We create similar containers for the residual estimator fields (which are
    // related to the shape fonctions F and I). These will be used to compute
    // the error in the energy norm on each cell.
    std::vector<Tensor<1, dim>> local_F_values_real(n_q_points);
    std::vector<Tensor<1, dim>> local_F_values_imag(n_q_points);
    std::vector<Tensor<1, dim>> local_I_values_real(n_q_points);
    std::vector<Tensor<1, dim>> local_I_values_imag(n_q_points);
    std::vector<Tensor<1, dim>> local_curl_F_values_real(n_q_points);
    std::vector<Tensor<1, dim>> local_curl_F_values_imag(n_q_points);
    std::vector<Tensor<1, dim>> local_curl_I_values_real(n_q_points);
    std::vector<Tensor<1, dim>> local_curl_I_values_imag(n_q_points);
    std::vector<Tensor<1, dim>> local_F_hat_values_real(n_face_q_points);
    std::vector<Tensor<1, dim>> local_F_hat_values_imag(n_face_q_points);
    std::vector<Tensor<1, dim>> local_I_hat_values_real(n_face_q_points);
    std::vector<Tensor<1, dim>> local_I_hat_values_imag(n_face_q_points);
    Tensor<1, dim>              local_n_I_hat_values_real;
    Tensor<1, dim>              local_n_I_hat_values_imag;

    // We reinitialize the estimated error per cell vector and update the ghost
    // values of the solution and residual estimator vectors to make sure that
    // they don't use any outdated values.
    estimated_error_per_cell.reinit(triangulation.n_active_cells());
    locally_relevant_solution_interior.update_ghost_values();
    locally_relevant_residual_estimator.update_ghost_values();

    // Here we initialize analytical solution object from the previously defined
    // classes.
    AnalyticalSolution_E<dim> analytical_solution_E(parameters);
    AnalyticalSolution_H<dim> analytical_solution_H(parameters);

    // Initialize variables that are used in the energy norm definition.
    const auto ikZr      = imag * k * Z_r;
    const auto conj_ikZr = conj(ikZr);
    const auto ikepsilon_Zr =
      imag * k / Z_r * (1. - imag * parameters.sigma_r / parameters.epsilon_r);
    const auto conj_ikepsilon_Zr = conj(ikepsilon_Zr);
    const auto kz_kZr            = k_z / (k * Z_r);
    const auto conj_kz_kZr       = conj(kz_kZr);

    // We create reusable variables to store the fields in complex form to ease
    // the calculations.
    std::complex<double> complex_shape_function;

    // We first loop over all the cells of our mesh and extract the solution in
    // the interior to compute its L2 error and the residual of each field on
    // each cell.
    for (const auto &cell : dof_handler_trial_interior.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            const unsigned int cell_index = cell->active_cell_index();
            fe_values_trial_interior.reinit(cell);
            const typename DoFHandler<dim>::active_cell_iterator cell_skeleton =
              cell->as_dof_handler_iterator(dof_handler_trial_skeleton);
            const typename DoFHandler<dim>::active_cell_iterator cell_test =
              cell->as_dof_handler_iterator(dof_handler_test);
            fe_values_test.reinit(cell_test);

            fe_values_trial_interior[extractor_E_real].get_function_values(
              locally_relevant_solution_interior, local_E_values_real);
            fe_values_trial_interior[extractor_E_imag].get_function_values(
              locally_relevant_solution_interior, local_E_values_imag);
            fe_values_trial_interior[extractor_H_real].get_function_values(
              locally_relevant_solution_interior, local_H_values_real);
            fe_values_trial_interior[extractor_H_imag].get_function_values(
              locally_relevant_solution_interior, local_H_values_imag);

            fe_values_test[extractor_E_real].get_function_values(
              locally_relevant_residual_estimator, local_F_values_real);
            fe_values_test[extractor_E_imag].get_function_values(
              locally_relevant_residual_estimator, local_F_values_imag);
            fe_values_test[extractor_H_real].get_function_values(
              locally_relevant_residual_estimator, local_I_values_real);
            fe_values_test[extractor_H_imag].get_function_values(
              locally_relevant_residual_estimator, local_I_values_imag);

            fe_values_test[extractor_E_real].get_function_curls(
              locally_relevant_residual_estimator, local_curl_F_values_real);
            fe_values_test[extractor_E_imag].get_function_curls(
              locally_relevant_residual_estimator, local_curl_F_values_imag);
            fe_values_test[extractor_H_real].get_function_curls(
              locally_relevant_residual_estimator, local_curl_I_values_real);
            fe_values_test[extractor_H_imag].get_function_curls(
              locally_relevant_residual_estimator, local_curl_I_values_imag);

            // Then we loop on the cell quadrature points because it is where we
            // will compute the error.
            const auto &quadrature_points =
              fe_values_trial_interior.get_quadrature_points();

            for (const unsigned int q_index :
                 fe_values_trial_interior.quadrature_point_indices())
              {
                const double JxW      = fe_values_trial_interior.JxW(q_index);
                const auto  &position = quadrature_points[q_index];

                // Loop over dimensions to compute each component of the error
                // because for a vector field is simply the sum of the squared
                // scalar L² errors of each component.
                for (unsigned int i = 0; i < dim; ++i)
                  {
                    // Calculate the L2 error for E
                    L2_error_E_real +=
                      pow((local_E_values_real[q_index][i] -
                           analytical_solution_E.value(position, i)),
                          2) *
                      JxW;
                    L2_error_E_imag +=
                      pow((local_E_values_imag[q_index][i] -
                           analytical_solution_E.value(position, i + dim)),
                          2) *
                      JxW;

                    // Calculate the L2 error for H
                    L2_error_H_real +=
                      pow((local_H_values_real[q_index][i] -
                           analytical_solution_H.value(position, i)),
                          2) *
                      JxW;
                    L2_error_H_imag +=
                      pow((local_H_values_imag[q_index][i] -
                           analytical_solution_H.value(position, i + dim)),
                          2) *
                      JxW;

                    // Calculate the Energy norm of Psi (the residual estimator)
                    // ||F||^2
                    complex_shape_function =
                      local_F_values_real[q_index][i] +
                      imag * local_F_values_imag[q_index][i];
                    estimated_error_per_cell(cell_index) +=
                      (pow(complex_shape_function.real(), 2) +
                       pow(complex_shape_function.imag(), 2)) *
                      JxW;

                    // ||I||^2
                    complex_shape_function =
                      local_I_values_real[q_index][i] +
                      imag * local_I_values_imag[q_index][i];
                    estimated_error_per_cell(cell_index) +=
                      (pow(complex_shape_function.real(), 2) +
                       pow(complex_shape_function.imag(), 2)) *
                      JxW;

                    // ||curl F - ikZr* I||^2
                    complex_shape_function =
                      local_curl_F_values_real[q_index][i] +
                      imag * local_curl_F_values_imag[q_index][i] -
                      conj_ikZr * (local_I_values_real[q_index][i] +
                                   imag * local_I_values_imag[q_index][i]);
                    estimated_error_per_cell(cell_index) +=
                      (pow(complex_shape_function.real(), 2) +
                       pow(complex_shape_function.imag(), 2)) *
                      JxW;

                    // ||curl I + ikepsilon_Zr* F||^2
                    complex_shape_function =
                      local_curl_I_values_real[q_index][i] +
                      imag * local_curl_I_values_imag[q_index][i] +
                      conj_ikepsilon_Zr *
                        (local_F_values_real[q_index][i] +
                         imag * local_F_values_imag[q_index][i]);
                    estimated_error_per_cell(cell_index) +=
                      (pow(complex_shape_function.real(), 2) +
                       pow(complex_shape_function.imag(), 2)) *
                      JxW;
                  }
              }

            // Then we loop on the face for the skeleton solution error and use
            // the same approach.
            for (const auto &face : cell->face_iterators())
              {
                // Reinitialization
                fe_face_values_trial_skeleton.reinit(cell_skeleton, face);
                fe_face_values_test.reinit(cell_test, face);

                // Extract local solution
                fe_face_values_trial_skeleton[extractor_E_real]
                  .get_function_values(locally_relevant_solution_skeleton,
                                       local_E_hat_values_real);
                fe_face_values_trial_skeleton[extractor_E_imag]
                  .get_function_values(locally_relevant_solution_skeleton,
                                       local_E_hat_values_imag);

                fe_face_values_trial_skeleton[extractor_H_real]
                  .get_function_values(locally_relevant_solution_skeleton,
                                       local_H_hat_values_real);
                fe_face_values_trial_skeleton[extractor_H_imag]
                  .get_function_values(locally_relevant_solution_skeleton,
                                       local_H_hat_values_imag);

                // The energy norm only has a face contribution the face where
                // we applied a Robin Boundary condition, i.e., the inlet port
                // and outlet port faces. So we only extract the residual
                // estimator solution on these faces to save some computational
                // time.
                if ((face->boundary_id() == 4) || (face->boundary_id() == 5))
                  {
                    // Extract local solution
                    fe_face_values_test[extractor_E_real].get_function_values(
                      locally_relevant_residual_estimator,
                      local_F_hat_values_real);
                    fe_face_values_test[extractor_E_imag].get_function_values(
                      locally_relevant_residual_estimator,
                      local_F_hat_values_imag);

                    fe_face_values_test[extractor_H_real].get_function_values(
                      locally_relevant_residual_estimator,
                      local_I_hat_values_real);
                    fe_face_values_test[extractor_H_imag].get_function_values(
                      locally_relevant_residual_estimator,
                      local_I_hat_values_imag);
                  }

                // Compute the L2 error
                const auto &quadrature_points =
                  fe_face_values_trial_skeleton.get_quadrature_points();

                for (const unsigned int &q_index :
                     fe_face_values_trial_skeleton.quadrature_point_indices())
                  {
                    const double JxW =
                      fe_face_values_trial_skeleton.JxW(q_index);
                    const auto &normal =
                      fe_face_values_trial_skeleton.normal_vector(q_index);
                    const auto &position = quadrature_points[q_index];

                    // This term will be needed for the energy norm calculation.
                    if ((face->boundary_id() == 4) ||
                        (face->boundary_id() == 5))
                      {
                        local_n_I_hat_values_real =
                          cross_product_3d(normal,
                                           local_I_hat_values_real[q_index]);
                        local_n_I_hat_values_imag =
                          cross_product_3d(normal,
                                           local_I_hat_values_imag[q_index]);
                      }

                    // Calculate the L2 error for E
                    for (unsigned int i = 0; i < dim; ++i)
                      {
                        L2_error_E_hat_real +=
                          pow((local_E_hat_values_real[q_index][i] -
                               analytical_solution_E.value(position, i)),
                              2) *
                          JxW;
                        L2_error_E_hat_imag +=
                          pow((local_E_hat_values_imag[q_index][i] -
                               analytical_solution_E.value(position, i + dim)),
                              2) *
                          JxW;

                        L2_error_H_hat_real +=
                          pow((local_H_hat_values_real[q_index][i] -
                               analytical_solution_H.value(position, i)),
                              2) *
                          JxW;
                        L2_error_H_hat_imag +=
                          pow((local_H_hat_values_imag[q_index][i] -
                               analytical_solution_H.value(position, i + dim)),
                              2) *
                          JxW;

                        if ((face->boundary_id() == 4) ||
                            (face->boundary_id() == 5))
                          {
                            // || n x I + kz_kZr * F||^2
                            complex_shape_function =
                              conj_kz_kZr *
                                (local_F_hat_values_real[q_index][i] +
                                 imag * local_F_hat_values_imag[q_index][i]) +
                              (local_n_I_hat_values_real[i] +
                               imag * local_n_I_hat_values_imag[i]);

                            estimated_error_per_cell(cell_index) +=
                              (pow(complex_shape_function.real(), 2) +
                               pow(complex_shape_function.imag(), 2)) *
                              JxW;
                          }
                      }
                  }
              }
          }
      }

    // Finally, we perform a global sum to obtain the final L2 errors and output
    // the results.
    const double L2_error_E_real_global =
      std::sqrt(Utilities::MPI::sum(L2_error_E_real, mpi_communicator));
    const double L2_error_E_imag_global =
      std::sqrt(Utilities::MPI::sum(L2_error_E_imag, mpi_communicator));
    const double L2_error_H_real_global =
      std::sqrt(Utilities::MPI::sum(L2_error_H_real, mpi_communicator));
    const double L2_error_H_imag_global =
      std::sqrt(Utilities::MPI::sum(L2_error_H_imag, mpi_communicator));
    const double L2_error_E_hat_real_global =
      std::sqrt(Utilities::MPI::sum(L2_error_E_hat_real, mpi_communicator));
    const double L2_error_E_hat_imag_global =
      std::sqrt(Utilities::MPI::sum(L2_error_E_hat_imag, mpi_communicator));
    const double L2_error_H_hat_real_global =
      std::sqrt(Utilities::MPI::sum(L2_error_H_hat_real, mpi_communicator));
    const double L2_error_H_hat_imag_global =
      std::sqrt(Utilities::MPI::sum(L2_error_H_hat_imag, mpi_communicator));
    pcout << "L2 E real part error is : " << L2_error_E_real_global
          << std::endl;
    pcout << "L2 E imag part error is : " << L2_error_E_imag_global
          << std::endl;
    pcout << "L2 H real part error is : " << L2_error_H_real_global
          << std::endl;
    pcout << "L2 H imag part error is : " << L2_error_H_imag_global
          << std::endl;
    pcout << "L2 E skeleton real part error is : " << L2_error_E_hat_real_global
          << std::endl;
    pcout << "L2 E skeleton imag part error is : " << L2_error_E_hat_imag_global
          << std::endl;
    pcout << "L2 H skeleton real part error is : " << L2_error_H_hat_real_global
          << std::endl;
    pcout << "L2 H skeleton imag part error is : " << L2_error_H_hat_imag_global
          << std::endl;

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        error_table.add_value("eL2_E_r", L2_error_E_real_global);
        error_table.add_value("eL2_E_i", L2_error_E_imag_global);
        error_table.add_value("eL2_H_r", L2_error_H_real_global);
        error_table.add_value("eL2_H_i", L2_error_H_imag_global);
        error_table.add_value("eL2_E_hat_r", L2_error_E_hat_real_global);
        error_table.add_value("eL2_E_hat_i", L2_error_E_hat_imag_global);
        error_table.add_value("eL2_H_hat_r", L2_error_H_hat_real_global);
        error_table.add_value("eL2_H_hat_i", L2_error_H_hat_imag_global);
      }
  }

  // The run function is the main loop of the program using all the previously
  // defined functions. It is also where the convergence rates are obtained
  // after all the refinement cycles. Note again, after solving the skeleton
  // system, we call the assembly function another time to solve for the
  // interior.
  template <int dim>
  void
  Time_Harmonic_Maxwell_DPG<dim>::run()
  {
    // Guarding against unsupported cases
    if (dim != 3)
      {
        throw std::runtime_error("Problems only works for 3D");
      }

    unsigned int cycle = 0;

    // We use the residual indicator as a stopping criteria for the refinement
    // cycles if the adaptive refinement flag is set to true.
    double residual_indicator = 1.0;
    // Stopping criteria for adaptive refinement which is the maximum residual
    // across all cells
    double stopping_criteria = 1e-5;

    pcout << "Running with Trilinos on "
          << Utilities::MPI::n_mpi_processes(mpi_communicator)
          << " MPI rank(s)..." << std::endl;

    while (residual_indicator > stopping_criteria)
      {
        pcout << "===========================================" << std::endl
              << "Cycle " << cycle << ':' << std::endl;

        if (cycle == 0)
          {
            make_grid(parameters.initial_refinements);
          }
        else
          {
            refine_grid();
          }

        const unsigned int n_cells = triangulation.n_global_active_cells();
        const double       cell_size =
          Utilities::MPI::max(GridTools::maximal_cell_diameter<dim>(
                                triangulation),
                              mpi_communicator);

        pcout << "Number of active cells: " << n_cells << std::endl;

        // Only rank 0 outputs / writes to table
        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
          {
            error_table.add_value("cycle", cycle);
            error_table.add_value("n_cells", n_cells);
            error_table.add_value("cell_size", cell_size);
          }

        setup_system();
        assemble_system();
        solve_skeleton();
        assemble_system(true); // Solve the interior problem
        calculate_L2_error();

        output_results(cycle);

        pcout << "===========================================" << std::endl;
        pcout << "End of cycle " << cycle << std::endl;
        pcout << "===========================================" << std::endl;

        computing_timer.print_summary();
        computing_timer.reset();
        pcout << std::endl;

        cycle++;
        const double max_residual =
          Utilities::MPI::max(estimated_error_per_cell.linfty_norm(),
                              mpi_communicator);
        residual_indicator = max_residual;

        if (cycle > parameters.n_refinements &&
            parameters.adaptive_refinement == false)
          break;
      }

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        error_table.evaluate_convergence_rates(
          "eL2_E_r", "n_cells", ConvergenceTable::reduction_rate_log2, 3);
        error_table.evaluate_convergence_rates(
          "eL2_E_i", "n_cells", ConvergenceTable::reduction_rate_log2, 3);
        error_table.evaluate_convergence_rates(
          "eL2_H_r", "n_cells", ConvergenceTable::reduction_rate_log2, 3);
        error_table.evaluate_convergence_rates(
          "eL2_H_i", "n_cells", ConvergenceTable::reduction_rate_log2, 3);
        error_table.evaluate_convergence_rates(
          "eL2_E_hat_r", "n_cells", ConvergenceTable::reduction_rate_log2, 3);
        error_table.evaluate_convergence_rates(
          "eL2_E_hat_i", "n_cells", ConvergenceTable::reduction_rate_log2, 3);
        error_table.evaluate_convergence_rates(
          "eL2_H_hat_r", "n_cells", ConvergenceTable::reduction_rate_log2, 3);
        error_table.evaluate_convergence_rates(
          "eL2_H_hat_i", "n_cells", ConvergenceTable::reduction_rate_log2, 3);
      }

    pcout << "===========================================" << std::endl;
    pcout << "Convergence table:" << std::endl;

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        error_table.write_text(std::cout);
      }
  }
} // end of namespace EM_DPG

// @sect3{The <code>main</code> function}

// This is the main function of the program.
int
main(int argc, char *argv[])
{
  const unsigned int dim = 3;

  try
    {
      using namespace dealii;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      // Initialize the parameters
      EM_DPG::Parameters<dim> parameters;
      ParameterAcceptor::initialize("parameters.prm");

      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        {
          std::cout << "==========================================="
                    << std::endl
                    << "Trial order: " << parameters.degree << std::endl
                    << "Test order: "
                    << parameters.delta_degree + parameters.degree << std::endl
                    << "==========================================="
                    << std::endl
                    << std::endl;
        }

      EM_DPG::Time_Harmonic_Maxwell_DPG<dim> maxwell_dpg(parameters);
      maxwell_dpg.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
