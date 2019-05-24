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
    // Dimension of the simulation
    unsigned int dimension;

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

    enum Variable {velocity, pressure};
    Variable variable;

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

  struct Testing
  {
    // Time measurement in the simulation. None, at each iteration, only at the end
    bool enabled;
    static void declare_parameters (ParameterHandler &prm);
    void parse_parameters (ParameterHandler &prm);
  };

  struct Restart
  {
    // Time measurement in the simulation. None, at each iteration, only at the end
    std::string filename;
    bool restart;
    bool checkpoint;
    unsigned int frequency;
    static void declare_parameters (ParameterHandler &prm);
    void parse_parameters (ParameterHandler &prm);
  };




}
#endif
