#ifndef lethe_vof_algebraic_reinitialization_scratch_data_h
#define lethe_vof_algebraic_reinitialization_scratch_data_h

#include <core/bdf.h>
#include <core/time_integration_utilities.h>

#include <solvers/vof.h>

#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/numerics/vector_tools.h>


/**
 * @brief Stores the information required by the assembly procedure for the VOF
 * algebraic reinitialization of the interface.
 * Computes the reinitialized phase and shape function values and gradients at
 * all Gauss quadrature points and stores them into arrays.
 * Serves as link between the evaluation at a Gauss point of a variables of
 * interest and their use in the assembly, which is carried out by the
 * @p assemblers (VOFAlgebraicReinitializationAssemblerBase).
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 */

template <int dim>
class VOFAlgebraicReinitializationScratchData
{
public:
  // TODO AMISHGA complete here
  /**
   * @brief Constructor of the scratch data object.
   * Creates the fe_values that will be used to fill the member variables.
   * Allocates the the memory necessary for all member variables, but does not
   * compute anything since that must be done at the cell level.
   *
   * // TODO AMISHGA check if this needed and same for the properties_manager
   * @param[in] simulation_control SimulationControl object that holds
   * information related to the control of the steady-state or transient
   * simulation.
   *
   * @param[in] properties_manager PhysicalPropertiesManager object that stores
   * physical property models.
   *
   * @param[in] fe_vof_algebraic_reinit FiniteElement object used for VOF
   * algebraic interface reinitialization.
   *
   * @param[in] quadrature Quadrature rule used for the assembly of the matrix
   * and the right-hand side.
   *
   * @param[in] mapping Mapping of the domain used when solving the VOF
   * equation.
   *
   * @param[in] fe_vof FiniteElement object used for VOF.
   */
  VOFAlgebraicReinitializationScratchData(
    const std::shared_ptr<SimulationControl> &simulation_control,
    const PhysicalPropertiesManager          &properties_manager,
    const FiniteElement<dim>                 &fe_vof_algebraic_reinit,
    const Quadrature<dim>                    &quadrature,
    const Mapping<dim>                       &mapping,
    const FiniteElement<dim>                 &fe_vof)
    : simulation_control(simulation_control)
    , properties_manager(properties_manager)
    , fe_values_vof_algebraic_reinit(mapping,
                                     fe_vof_algebraic_reinit,
                                     quadrature,
                                     update_values | update_gradients |
                                       update_quadrature_points |
                                       update_JxW_values)
    , fe_values_vof(mapping,
                    fe_vof,
                    quadrature,
                    update_values | update_gradients)
  {
    allocate();
  }

  VOFAlgebraicReinitializationScratchData(
    VOFAlgebraicReinitializationScratchData<dim> &sd)
    : simulation_control(sd.simulation_control)
    , properties_manager(sd.properties_manager)
    , fe_values_vof_algebraic_reinit(
        sd.fe_values_vof_algebraic_reinit.get_mapping(),
        sd.fe_values_vof_algebraic_reinit.get_fe(),
        sd.fe_values_vof_algebraic_reinit.get_quadrature(),
        update_values | update_gradients | update_quadrature_points |
          update_JxW_values)
    , fe_values_vof(sd.fe_values_vof.get_mapping(),
                    sd.fe_values_vof.get_fe(),
                    sd.fe_values_vof.get_quadrature(),
                    update_values | update_gradients)
  {
    allocate();
  }
  /**
   * @brief Allocates the memory necessary memory for all members of the
   * scratch.
   */
  void
  allocate();


  /**
   * @brief Reinitializes the content of the
   * VOFAlgebraicReinitializationScratchData object for a given cell using
   * FeValues and contents of solution and previous solution arrays.
   *
   * @tparam VectorType Vector type used for the solvers.
   *
   * @param[in] cell Cell over which the assembly is carried.
   *
   * @param[in] current_solution Present value of the solution for the
   * reinitialized phase fraction.
   *
   * @param[in] previous_solutions Vector of previous time step solution values
   * for the reinitialized phase fraction.
   */
  template <typename VectorType>
  void
  reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
         const VectorType                                     &current_solution,
         const std::vector<VectorType> &previous_solutions)
  {
    fe_values_vof_algebraic_reinit.reinit(cell);
    this->quadrature_points =
      fe_values_vof_algebraic_reinit
        .get_quadrature_points(); // TODO AMISHGA check if needed
    auto &fe_vof_algebraic_reinit = fe_values_vof_algebraic_reinit.get_fe();

    // TODO AMISHGA refactor to utilities (inline function)
    if (dim == 2)
      this->cell_size =
        std::sqrt(4. * cell->measure() / M_PI) / fe_vof_algebraic_reinit.degree;
    else if (dim == 3)
      this->cell_size = pow(6 * cell->measure() / M_PI, 1. / 3.) /
                        fe_vof_algebraic_reinit.degree;

    // Gather present solutions
    fe_values_vof_algebraic_reinit.get_function_values(
      current_solution, this->present_phase_algebraic_reinit_values);
    fe_values_vof_algebraic_reinit.get_function_gradients(
      current_solution, this->present_phase_algebraic_reinit_gradients);

    // Gather previous solutions for time stepping
    for (unsigned int p = 0; p < previous_solutions.size(); ++p)
      {
        fe_values_vof_algebraic_reinit.get_present_fe_values(
          previous_solutions[p],
          this->previous_phase_algebraic_reinit_values[p]);
      }
  }

  /**
   * @brief Reinitialize the cell with the VOF solution calculated with the VOF
   * advection equation.
   *
   * @tparam VectorType Vector type used for the solvers.
   *
   * @param[in] cell Cell over which the assembly is carried.
   *
   * @param[in] current_solution Present value of the solution for the
   * VOF phase fraction.
   *
   * @param[in] previous_solutions Vector of previous time step solution values
   * for the VOF phase fraction.
   */
  template <typename VectorType>
  void
  reinit_vof(const typename DoFHandler<dim>::active_cell_iterator &cell,
             const VectorType              &current_solution,
             const std::vector<VectorType> &previous_solutions)
  {
    fe_values_vof.reinit(cell);
    fe_values_vof.get_function_values(current_solution,
                                      present_vof_phase_values);
    fe_values_vof.get_function_gradients(current_solution,
                                         present_vof_phase_gradients);

    // TODO AMISHGA check if necessary

    //    // Gather previous values
    //    for (unsigned int p = 0; p < previous_solutions.size(); ++p)
    //      {
    //        fe_values_vof.get_function_gradients(
    //          previous_solutions[p], this->previous_vof_phase_gradients);
    //      }
  }

  // TODO AMISHGA check if variables are necessary and do a cleanup
  const std::shared_ptr<SimulationControl> simulation_control;

  // Physical properties
  const PhysicalPropertiesManager      properties_manager;
  std::map<field, std::vector<double>> fields;

  // FEValues for the VOF algebraic reinitialization problem
  FEValues<dim> fe_values_vof_algebraic_reinit;
  unsigned int  n_dofs;
  unsigned int  n_q_points;
  double        cell_size;

  // Quadrature
  std::vector<double>     JxW;
  std::vector<Point<dim>> quadrature_points;

  // VOF algebraic reinitialization values
  std::vector<double>              present_phase_algebraic_reinit_values;
  std::vector<Tensor<1, dim>>      present_phase_algebraic_reinit_gradients;
  std::vector<std::vector<double>> previous_phase_algebraic_reinit_values;
  //  std::vector<Tensor<1, dim>> previous_phase_algebraic_reinit_gradients;

  // Shape functions
  std::vector<std::vector<double>>         ksi;
  std::vector<std::vector<Tensor<1, dim>>> grad_ksi;
  //  std::vector<std::vector<Tensor<2, dim>>> hess_ksi;
  //  std::vector<std::vector<double>>         laplacian_ksi;


  /**
   * Scratch component for VOF
   */
  std::vector<double>         present_vof_phase_values;
  std::vector<Tensor<1, dim>> present_vof_phase_gradients;
  //  std::vector<std::vector<Tensor<1, dim>>> previous_vof_phase_gradients
  //  //TODO need to extrapolate ?
  //    std::vector<std::vector<double>> previous_vof_phase_values;
  FEValues<dim> fe_values_vof;
};



#endif
