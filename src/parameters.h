#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/conditional_ostream.h>

#ifndef LETHE_GLS_PARAMETERS_H
#define LETHE_GLS_PARAMETERS_H

using namespace dealii;

namespace Parameters
{

  void declareAllParameters(ParameterHandler &prm);


  struct SimulationControl
  {
    // Method used for time progression (steady, unsteady)
    enum TimeSteppingMethod { steady, backward, bdf2};
    TimeSteppingMethod method;

    // Initial time step
    double dt;

    // End time
    double timeEnd;

    // Adaptative time stepping
    bool   adapt;

    // Max CFL
    double maxCFL;

    // Number of mesh adaptation (steady simulations)
    unsigned int nbMeshAdapt;

    // Folder for simulation output
    std::string output_folder;

    // Prefix for simulation output
    std::string output_name;

    // Frequency of the output
    unsigned int outputFrequency;

    // Subdivsions of the results in the output
    unsigned int subdivision;

    static void declare_parameters (ParameterHandler &prm);
    void parse_parameters (ParameterHandler &prm);
  };


  struct PhysicalProperties
  {
    // Kinematic viscosity (mu/rho)
    double viscosity;

    static void declare_parameters (ParameterHandler &prm);
    void parse_parameters (ParameterHandler &prm);
  };


  struct Timer
  {
    // Time measurement in the simulation. None, at each iteration, only at the end
    enum Type { none, iteration, end};
    Type type;
    static void declare_parameters (ParameterHandler &prm);
    void parse_parameters (ParameterHandler &prm);
  };

  struct AnalyticalSolution
  {
    // Residual precision
    unsigned int errorPrecision;
    static void declare_parameters (ParameterHandler &prm);
    void parse_parameters (ParameterHandler &prm);
  };

  struct Forces
  {
    // Type of verbosity for the iterative solver
    enum  Verbosity { quiet, verbose };
    Verbosity verbosity;

    // Enable force post-processing
    bool calculate_force;

    // Enable torque post-processing
    bool calculate_torque;

    // Frequency of the output
    unsigned int calculation_frequency;

    // Frequency of the output
    unsigned int output_frequency;

    // Output precision
    unsigned int output_precision;

    // Display precision
    unsigned int display_precision;

    // Prefix for simulation output
    std::string force_output_name;

    // Prefix for the torque output
    std::string torque_output_name;

    static void declare_parameters (ParameterHandler &prm);
    void parse_parameters (ParameterHandler &prm);
  };

  struct FEM
  {
    // Interpolation order velocity
    unsigned int velocityOrder;

    // Interpolation order pressure
    unsigned int pressureOrder;

    // Apply high order mapping everywhere
    bool qmapping_all;

    static void declare_parameters (ParameterHandler &prm);
    void parse_parameters (ParameterHandler &prm);
  };

  struct NonLinearSolver
  {
    // Type of verbosity for the iterative solver
    enum  Verbosity { quiet, verbose };
    Verbosity verbosity;

    // Tolerance
    double tolerance;

    // Maximal number of iterations for the Newton solver
    unsigned int maxIterations;

    // Residual precision
    unsigned int display_precision;

    static void declare_parameters (ParameterHandler &prm);
    void parse_parameters (ParameterHandler &prm);
  };


  struct LinearSolver
  {
    // Type of linear solver
    enum SolverType { gmres, bicgstab, amg};
    SolverType solver;

    // Type of verbosity for the iterative solver
    enum  Verbosity { quiet, verbose };
    Verbosity verbosity;

    // Residual precision
    unsigned int residual_precision;

    //Right hand side correction for constant pressure, does not work well right now
    bool   rhsCorr;

    // Relative residual of the iterative solver
    double relative_residual;

    // Minimum residual of the iterative solver
    double minimum_residual;

    // Maximum number of iterations
    int max_iterations;

    // ILU or ILUT fill
    double ilu_fill;

    // ILU or ILUT absolute tolerance (1e-3 works well)
    double ilu_atol;

    // ILU or ILUT relative tolerance (1.00 works well)
    double ilu_rtol;

    // AMG aggregation threshold
    double amg_aggregation_threshold;

    // AMG number of cycles
    unsigned int amg_n_cycles;

    // AMG W_cycle
    bool amg_w_cycles;

    // AMG Smoother sweeps
    unsigned int amg_smoother_sweeps;

    // AMG Smoother overalp
    unsigned int amg_smoother_overlap;

    static void declare_parameters (ParameterHandler &prm);
    void parse_parameters (ParameterHandler &prm);
  };

  struct Mesh
  {
    // GMSH or dealii primitive
    enum Type {gmsh, primitive};
    Type type;

    // File name of the mesh
    std::string fileName;

    // Initial refinement level of primitive mesh
    unsigned int initialRefinement;

    static void declare_parameters (ParameterHandler &prm);
    void parse_parameters (ParameterHandler &prm);
  };



  struct MeshAdaptation
  {

    // Type of mesh adaptation
    enum Type {none, uniform,kelly};
    Type type;

    // Decision factor for KELLY refinement (number or fraction)
    enum FractionType {number,fraction};
    FractionType fractionType;

    // Maximum number of elements
    unsigned int maxNbElements;

    // Maximum refinement level
    unsigned int maxRefLevel;

    // Maximum refinement level
    unsigned int minRefLevel;

    // Refinement after frequency iter
    unsigned int frequency;

    // Refinement fractioni havent used ILUT much)
    double fractionRefinement;

    // Coarsening fraction
    double fractionCoarsening;

    static void declare_parameters (ParameterHandler &prm);
    void parse_parameters (ParameterHandler &prm);
  };

  
  enum BoundaryType {noslip, slip, function};

  template <int dim>
  class BoundaryFunction
  {
  public:
    // Velocity components
    Functions::ParsedFunction<dim> u;
    Functions::ParsedFunction<dim> v;
    Functions::ParsedFunction<dim> w;

    // Point for the center of rotation
    Point<dim>                     cor;
  };

  template <int dim>
  class BoundaryConditions
  {
  public:
    // List of boundary type for each number
    std::vector<BoundaryType> type;

    // Functions for (u,v,w) for all boundaries
    BoundaryFunction<dim> *bcFunctions;

    // Number of boundary conditions
    unsigned int size;
    unsigned int max_size;

    void parse_boundary (ParameterHandler &prm, unsigned int i_bc);
    void declareDefaultEntry (ParameterHandler &prm, unsigned int i_bc);
    void declare_parameters (ParameterHandler &prm);
    void parse_parameters (ParameterHandler &prm);
  };
  template <int dim>
  void BoundaryConditions<dim>::declareDefaultEntry (ParameterHandler &prm, unsigned int i_bc)
  {
    prm.declare_entry("type", "noslip",
                      Patterns::Selection("noslip|slip|function"),
                      "Type of boundary conditoin"
                      "Choices are <noslip|slip|function>.");

    prm.enter_subsection("u");
    bcFunctions[i_bc].u.declare_parameters(prm,1);
    prm.set("Function expression","0");
    prm.leave_subsection();

    prm.enter_subsection("v");
    bcFunctions[i_bc].v.declare_parameters(prm,1);
    prm.set("Function expression","0");
    prm.leave_subsection();

    prm.enter_subsection("w");
    bcFunctions[i_bc].w.declare_parameters(prm,1);
    prm.set("Function expression","0");
    prm.leave_subsection();

    prm.enter_subsection("cor");
    prm.declare_entry("x","0",Patterns::Double(),"X COR");
    prm.declare_entry("y","0",Patterns::Double(),"Y COR");
    prm.declare_entry("z","0",Patterns::Double(),"Z COR");
    prm.leave_subsection();

  }

  template <int dim>
  void BoundaryConditions<dim>::parse_boundary (ParameterHandler &prm, unsigned int i_bc)
  {
    const std::string op = prm.get("type");
    if (op == "noslip")
      type[i_bc] = noslip;
    if (op == "slip")
      type[i_bc] = slip;
    if (op == "function")
    {
      type[i_bc] = function;
      prm.enter_subsection("u");
      bcFunctions[i_bc].u.parse_parameters(prm);
      prm.leave_subsection();

      prm.enter_subsection("v");
      bcFunctions[i_bc].v.parse_parameters(prm);
      prm.leave_subsection();

      prm.enter_subsection("w");
      bcFunctions[i_bc].w.parse_parameters(prm);
      prm.leave_subsection();

      prm.enter_subsection("cor");
      bcFunctions[i_bc].cor[0]=prm.get_double("x");
      bcFunctions[i_bc].cor[1]=prm.get_double("y");
      if (dim==3) bcFunctions[i_bc].cor[2]=prm.get_double("z");
      prm.leave_subsection();
    }
  }

  template <int dim>
  void BoundaryConditions<dim>::declare_parameters (ParameterHandler &prm)
  {
    max_size=4;

    prm.enter_subsection("boundary conditions");
    {
      prm.declare_entry("number", "0",
                        Patterns::Integer(),
                        "Number of boundary conditions");
      type.resize(max_size);
      bcFunctions = new BoundaryFunction<dim>[max_size];

      prm.enter_subsection("bc 0");
      {
        declareDefaultEntry(prm,0);
      }
      prm.leave_subsection();

      prm.enter_subsection("bc 1");
      {
        declareDefaultEntry(prm,1);
      }
      prm.leave_subsection();

      prm.enter_subsection("bc 2");
      {
        declareDefaultEntry(prm,2);
      }
      prm.leave_subsection();

      prm.enter_subsection("bc 3");
      {
        declareDefaultEntry(prm,3);
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }

  template <int dim>
  void BoundaryConditions<dim>::parse_parameters (ParameterHandler &prm)
  {
    prm.enter_subsection("boundary conditions");
    {
      size = prm.get_integer("number");
      type.resize(size);


      if (size>=1)
      {
        prm.enter_subsection("bc 0");
        {
          parse_boundary(prm,0);
        }
        prm.leave_subsection();
      }
      if (size>=2)
      {
        prm.enter_subsection("bc 1");
        {
          parse_boundary(prm,1);
        }
        prm.leave_subsection();
      }
      if (size>=3)
      {
        prm.enter_subsection("bc 2");
        {
          parse_boundary(prm,2);
        }
        prm.leave_subsection();
      }
      if (size>=4)
      {
        prm.enter_subsection("bc 3");
        {
          parse_boundary(prm,3);
        }
        prm.leave_subsection();
      }
    }
    prm.leave_subsection();
  }

  // Type of initial conditions
  enum InitialConditionType {none, L2projection, viscous, nodal};

  template <int dim>
  class InitialConditions
  {
  public:
    InitialConditions():
      uvwp(dim+1)
    {}

    InitialConditionType type;

    // Artificial viscosity
    double viscosity;

    // Velocity components
    Functions::ParsedFunction<dim> uvwp;

    void declare_parameters (ParameterHandler &prm);
    void parse_parameters (ParameterHandler &prm);
  };




  template <int dim>
  void InitialConditions<dim>::declare_parameters (ParameterHandler &prm)
  {
    prm.enter_subsection("initial conditions");
    {
      prm.declare_entry("type", "none",
                        Patterns::Selection("none|L2projection|viscous|nodal"),
                        "Type of initial condition"
                        "Choices are <none|L2projection|viscous|nodal>.");

      prm.enter_subsection("uvwp");
      uvwp.declare_parameters(prm,dim);
      if (dim==2) prm.set("Function expression","0; 0; 0");
      if (dim==3) prm.set("Function expression","0; 0; 0; 0");
      prm.leave_subsection();

      prm.declare_entry("viscosity", "1",Patterns::Double(),"viscosity for viscous initial conditions");
    }
    prm.leave_subsection();
  }

  template <int dim>
  void InitialConditions<dim>::parse_parameters (ParameterHandler &prm)
  {
    prm.enter_subsection("initial conditions");
    {
      const std::string op = prm.get("type");
      if (op == "none")
        type = none;
      if (op == "L2projection")
        type = L2projection;
      if (op == "viscous")
        type = viscous;
      if (op== "nodal")
        type = nodal;

      viscosity = prm.get_double("viscosity");
      prm.enter_subsection("uvwp");
      uvwp.parse_parameters(prm);
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }

  FEM getFEMParameters2D(std::string);
  FEM getFEMParameters3D(std::string);
}
#endif
