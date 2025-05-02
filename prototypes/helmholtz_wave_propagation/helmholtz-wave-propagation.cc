// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/symmetric_tensor.h>
 
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
 
#include <iostream>
#include <fstream>
#include <complex>

using namespace dealii;

struct Settings
{
  bool
  try_parse(const std::string &prm_filename);

  enum ProblemType
  {
    mms,
    channel//,
    // channel_without_walls
  };

  enum MMSProblem
  {
    sine_wave,
    plane_wave
  };

  enum MMSGeometry
  {
    cube,
    ball
  };

  enum MMSBoundaryType
  {
    dirichlet,
    neumann,
    robin
  };

  ProblemType   problem_type;
  unsigned int  dimensions;
  unsigned int  fe_order;
  unsigned int  initial_refinement;

  MMSProblem      mms_problem;
  MMSGeometry     mms_geometry;
  MMSBoundaryType mms_boundary_type;  

  unsigned int x_cell_ratio;
  unsigned int y_cell_ratio;

  double  domain_length_ratio;
  double  transducer_length_ratio;
  double  wall_thickness_ratio;

  double  pml_thickness_ratio;
  double  pml_stretching;
  double  pml_attenuation;
  double  pml_exponent;

  double  transducer_amplitude;
  double  transducer_phase;
  double  transducer_frequency;

  double  wall_density;
  double  wall_speed_of_sound;
  double  wall_impedance;

  double  dynamic_viscosity;
  double  bulk_viscosity;
  double  density;
  double  speed_of_sound;
  double  impedance;

  double particle_radius;
  double particle_density;
  double particle_speed_of_sound;

  bool          convergence_study;
  unsigned int  n_cycles;

  std::string   output_name;
  std::string   output_path;
  unsigned int  n_timesteps;

  double wavelength;
  double wall_wavelength;

  double angular_frequency;

  std::vector<double>               look_up_density;
  std::vector<double>               look_up_speed_of_sound;
  std::vector<std::complex<double>> look_up_xi;
};

bool
Settings::try_parse(const std::string &prm_filename)
{
  ParameterHandler prm;

  // General Case parameters
  prm.declare_entry("problem type",
                    "channel",
                    Patterns::Selection("mms|channel|channel_without_walls"),
                    "Problem type <mms|channel|channel_without_walls>");
  prm.declare_entry("dim",
                    "2",
                    Patterns::Integer(),
                    "The problem dimension <1|2|3>");
  prm.declare_entry("element order",
                    "2",
                    Patterns::Integer(),
                    "Order of FE element <1|2|3|4>");
  prm.declare_entry("initial refinement",
                    "3",
                    Patterns::Integer(),
                    "Initial refinement level");

  // MMS Parameters
  prm.declare_entry("mms problem type",
                    "sine_wave",
                    Patterns::Selection("sine_wave|plane_wave"),
                    "MMS problem type <sine_wave|plane_wave>");
  prm.declare_entry("mms geometry type",
                    "cube",
                    Patterns::Selection("cube|ball"),
                    "MMS geometry type <cube|ball>");
  prm.declare_entry("mms boundary type",
                    "dirichlet",
                    Patterns::Selection("dirichlet|neumann|robin"),
                    "MMS boundary type <dirichlet|neumann|robin>");
  prm.declare_entry("mms n wavelengths",
                    "1.",
                    Patterns::Double(),
                    "Number of wavelengths in the domain");
  prm.declare_entry("mms frequency",
                    "1.",
                    Patterns::Double(),
                    "Wave frequency");

  // Channel Resolution parameters
  prm.declare_entry("x cell ratio",
                    "2",
                    Patterns::Integer(),
                    "Numbers of cells per wavelength in x (along channel)");
  prm.declare_entry("y cell ratio",
                    "2",
                    Patterns::Integer(),
                    "Numbers of cells per wavelength in y or z (across channel)");
  
  // Channel Domain parameters
  prm.declare_entry("domain length ratio",
                    "1",
                    Patterns::Double(),
                    "Domain length ratio");
  prm.declare_entry("transducer length ratio",
                    "0.2",
                    Patterns::Double(),
                    "Transducer length ratio");
  prm.declare_entry("wall thickness ratio",
                    "0.25",
                    Patterns::Double(),
                    "Wall thickness ratio");
  
  // PML parameters
  prm.declare_entry("pml thickness ratio",
                    "1.0",
                    Patterns::Double(),
                    "PML thickness ratio");
  prm.declare_entry("pml stretching",
                    "1.0",
                    Patterns::Double(),
                    "PML stretching");
  prm.declare_entry("pml attenuation",
                    "1.0",
                    Patterns::Double(),
                    "PML attenuation");
  prm.declare_entry("pml exponent",
                    "2.0",
                    Patterns::Double(),
                    "PML exponent <1|2|3>");

  // Transducer parameters
  prm.declare_entry("transducer amplitude",
                    "1e-4",
                    Patterns::Double(),
                    "Transducer amplitude");
  prm.declare_entry("transducer phase",
                    "0.0",
                    Patterns::Double(),
                    "Transducer phase");
  prm.declare_entry("transducer frequency",
                    "30e3",
                    Patterns::Double(),
                    "Transducer frequency");

  prm.declare_entry("wall boundary type",
                    "sound_hard",
                    Patterns::Selection("sound_hard|impedance"),
                    "Wall boundary type <sound_hard|impedance>");
  prm.declare_entry("wall density",
                    "2500.",
                    Patterns::Double(),
                    "Wall density");
  prm.declare_entry("wall speed of sound",
                    "4500.",
                    Patterns::Double(),
                    "Wall speed of sound");

  prm.declare_entry("particle radius",
                    "5e-3",
                    Patterns::Double(),
                    "Particle radius in m");
  prm.declare_entry("particle density",
                    "1000.",
                    Patterns::Double(),
                    "Particle density in kg/m^3");
  prm.declare_entry("particle speed of sound",
                    "1500.",
                    Patterns::Double(),
                    "Particle speed of sound in m/s");

  prm.declare_entry("dynamic viscosity",
                    "1e-3",
                    Patterns::Double(),
                    "Dynamic viscosity");
  prm.declare_entry("bulk viscosity",
                    "3e-3",
                    Patterns::Double(),
                    "Bulk viscosity");
  prm.declare_entry("density",
                    "1000.",
                    Patterns::Double(),
                    "Density");
  prm.declare_entry("speed of sound",
                    "1500.",
                    Patterns::Double(),
                    "Speed of sound");

  prm.declare_entry("convergence study",
                    "false",
                    Patterns::Bool(),
                    "Compute convergence rates <true|false>");
  prm.declare_entry("number of cycles",
                    "1",
                    Patterns::Integer(),
                    "Number of refinement cycles");

  prm.declare_entry("output name",
                    "solution",
                    Patterns::FileName(),
                    "Name for vtu files");
  prm.declare_entry("output path",
                    "./",
                    Patterns::FileName(),
                    "Path for vtu output files");
  prm.declare_entry("n timesteps",
                    "100",
                    Patterns::Integer(),
                    "Number of time steps");

  if (prm_filename.size() == 0)
    {
      std::cout
        << "****  Error: No input file provided!\n"
        << "****  Error: Call this program as './matrix_based_non_linear_poisson input.prm\n"
        << '\n'
        << "****  You may want to use one of the input files in this\n"
        << "****  directory, or use the following default values\n"
        << "****  to create an input file:\n";
#if DEAL_II_VERSION_GTE(9, 7, 0)
        prm.print_parameters(std::cout, ParameterHandler::DefaultStyle);
#else
        prm.print_parameters(std::cout, ParameterHandler::Text);
#endif
      return false;
    }

  try
    {
      prm.parse_input(prm_filename);
    }
  catch (std::exception &e)
    {
      std::cerr << e.what() << std::endl;
      return false;
    }

  if (prm.get("problem type") == "mms")
    this->problem_type = mms;
  else if (prm.get("problem type") == "channel")
    this->problem_type = channel;
  // else if (prm.get("problem type") == "channel_without_walls")
  //   this->problem_type = channel_without_walls;
  else
    AssertThrow(false, ExcNotImplemented());

  if (prm.get("mms problem type") == "sine_wave")
    this->mms_problem = sine_wave;
  else if (prm.get("mms problem type") == "plane_wave")
    this->mms_problem = plane_wave;
  else
    AssertThrow(false, ExcNotImplemented());

  if (prm.get("mms geometry type") == "cube")
    this->mms_geometry = cube;
  else if (prm.get("mms geometry type") == "ball")
    this->mms_geometry = ball;
  else
    AssertThrow(false, ExcNotImplemented());

  if (prm.get("mms boundary type") == "dirichlet")
    this->mms_boundary_type = dirichlet;
  else if (prm.get("mms boundary type") == "neumann")
    this->mms_boundary_type = neumann;
  else if (prm.get("mms boundary type") == "robin")
    this->mms_boundary_type = robin;
  else
    AssertThrow(false, ExcNotImplemented());

  this->dimensions         = prm.get_integer("dim");
  this->fe_order           = prm.get_integer("element order");
  this->initial_refinement = prm.get_integer("initial refinement");

  this->x_cell_ratio = prm.get_integer("x cell ratio");
  this->y_cell_ratio = prm.get_integer("y cell ratio");

  this->domain_length_ratio     = prm.get_double("domain length ratio");
  this->transducer_length_ratio = prm.get_double("transducer length ratio");
  this->wall_thickness_ratio    = prm.get_double("wall thickness ratio");
  
  this->pml_thickness_ratio = prm.get_double("pml thickness ratio");
  this->pml_stretching      = prm.get_double("pml stretching");
  this->pml_attenuation     = prm.get_double("pml attenuation");
  this->pml_exponent        = prm.get_double("pml exponent");

  this->transducer_amplitude  = prm.get_double("transducer amplitude");
  this->transducer_phase      = M_PI / 180. * prm.get_double("transducer phase");
  double transducer_frequency = prm.get_double("transducer frequency");

  this->wall_density        = prm.get_double("wall density");
  this->wall_speed_of_sound = prm.get_double("wall speed of sound");
  this->wall_impedance      = this->wall_density * this->wall_speed_of_sound;

  this->dynamic_viscosity = prm.get_double("dynamic viscosity");
  this->bulk_viscosity    = prm.get_double("bulk viscosity");
  this->density           = prm.get_double("density");
  this->speed_of_sound    = prm.get_double("speed of sound"); 
  this->impedance         = this->density * this->speed_of_sound;

  this->particle_radius        = prm.get_double("particle radius");
  this->particle_density       = prm.get_double("particle density");
  this->particle_speed_of_sound = prm.get_double("particle speed of sound");
  
  this->convergence_study = prm.get_bool("convergence study");
  this->n_cycles          = prm.get_integer("number of cycles");
  
  this->output_name = prm.get("output name");
  this->output_path = prm.get("output path");
  this->n_timesteps = prm.get_integer("n timesteps");

  this->wavelength = this->speed_of_sound / transducer_frequency;
  this->wall_wavelength = this->wall_speed_of_sound / transducer_frequency;

  std::complex<double> i(0., 1.);

  // Setup wave number based on the problem type
  if (this->problem_type == mms)
  {
    double wavelength = 1.;
    double wave_speed = 1.;

    this->look_up_density = {1.0};
    this->look_up_speed_of_sound = {wave_speed};

    double wavenumber = 2. * M_PI / wavelength;

    this->angular_frequency = wavenumber / wave_speed;
    this->look_up_xi = {wavenumber * wavenumber};
  }
  else
  {
    this->angular_frequency = 2. * M_PI * transducer_frequency;

    double beta = this->bulk_viscosity / (this->dynamic_viscosity + 1e-6) + 1./3.;
    double kinematic_viscosity = this->dynamic_viscosity / this->density; 
    double viscous_damping = this->angular_frequency * (1. + beta) 
          * kinematic_viscosity / (this->speed_of_sound * this->speed_of_sound);

    this->look_up_density = {this->density, this->wall_density};
    this->look_up_speed_of_sound = {this->speed_of_sound, this->wall_speed_of_sound};

    this->look_up_xi = {
            // fluid wavenumber^2
            this->angular_frequency * this->angular_frequency / 
            (this->speed_of_sound * this->speed_of_sound * 
            (1. + viscous_damping * i)),
            // wall wavenumber^2
            this->angular_frequency * this->angular_frequency /
            (this->wall_speed_of_sound * this->wall_speed_of_sound)
            };
  }

  return true;
}

template <int dim>
class ComputePressureIntensity : public DataPostprocessorScalar<dim>
{
public:
  ComputePressureIntensity(const Settings &parameters)
    : DataPostprocessorScalar<dim>("pressure_intensity", update_values)
    , parameters(parameters)
  {}

  virtual void
  evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &scalar_potential,
                         std::vector<Vector<double>> &computed_quantities) const override
  {
    AssertDimension (scalar_potential.solution_values.size(), computed_quantities.size());

    // Get the material properties for the current cell
    const typename DoFHandler<dim>::cell_iterator current_cell =
      scalar_potential.template get_cell<dim>();
    const double density = parameters.look_up_density[current_cell->material_id()];
    const double speed_of_sound = parameters.look_up_speed_of_sound[current_cell->material_id()];
    const std::complex<double>  xi = parameters.look_up_xi[current_cell->material_id()];

    for (unsigned int p = 0; p < scalar_potential.solution_values.size(); ++p)
      {
        AssertDimension (scalar_potential.solution_values[p].size(), 2);
        AssertDimension (computed_quantities[p].size(), 1);
        
        // Get complex pressure amplitude from complex scalar potential amplitude
        std::complex<double> i(0., 1.);
        std::complex<double> phi(scalar_potential.solution_values[p](0), 
                                  scalar_potential.solution_values[p](1));
        std::complex<double> pressure = density * speed_of_sound * speed_of_sound
                                        * xi * phi / (parameters.angular_frequency * i);

        // Compute the intensity of the pressure (magnitude of the complex pressure)
        computed_quantities[p](0) = std::abs(pressure);
      }
  }

private:
  const Settings &parameters;
};

template <int dim>
class ComputePressurePhase : public DataPostprocessorScalar<dim>
{
public:
  ComputePressurePhase(const Settings &parameters)
    : DataPostprocessorScalar<dim>("pressure_phase", update_values)
    , parameters(parameters)
  {}

  virtual void
  evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &scalar_potential,
                         std::vector<Vector<double>> &computed_quantities) const override
  {
    AssertDimension (scalar_potential.solution_values.size(), computed_quantities.size());

    // Get the material properties for the current cell
    const typename DoFHandler<dim>::cell_iterator current_cell =
      scalar_potential.template get_cell<dim>();
    const double density = parameters.look_up_density[current_cell->material_id()];
    const double speed_of_sound = parameters.look_up_speed_of_sound[current_cell->material_id()];
    const std::complex<double> xi = parameters.look_up_xi[current_cell->material_id()];

    for (unsigned int p = 0; p < scalar_potential.solution_values.size(); ++p)
      {
        AssertDimension (scalar_potential.solution_values[p].size(), 2);
        AssertDimension (computed_quantities[p].size(), 1);
        
        // Get complex pressure amplitude from complex scalar potential amplitude 
        std::complex<double> i(0., 1.);
        std::complex<double> phi(scalar_potential.solution_values[p](0), 
                                  scalar_potential.solution_values[p](1));
        std::complex<double> pressure = density * speed_of_sound * speed_of_sound
                                        * xi * phi / (parameters.angular_frequency * i);

        // Compute the phase of the pressure (argument of the complex pressure)
        computed_quantities[p](0) = std::arg(pressure);
      }
  }

private:
  const Settings &parameters;
};  

template <int dim>
class ComputeTimeHarmonicPressure : public DataPostprocessorScalar<dim>
{
public :
  ComputeTimeHarmonicPressure(const Settings &parameters)
    : DataPostprocessorScalar<dim>("time_harmonic_pressure", update_values)
    , parameters(parameters)
  {}

  virtual void
  evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &scalar_potential,
                         std::vector<Vector<double>> &computed_quantities) const override
  {
    AssertDimension (scalar_potential.solution_values.size(), computed_quantities.size());

    // Get the material properties for the current cell
    const typename DoFHandler<dim>::cell_iterator current_cell =
      scalar_potential.template get_cell<dim>();
    const double density = parameters.look_up_density[current_cell->material_id()];
    const double speed_of_sound = parameters.look_up_speed_of_sound[current_cell->material_id()];
    const std::complex<double> xi = parameters.look_up_xi[current_cell->material_id()];

    for (unsigned int p = 0; p < scalar_potential.solution_values.size(); ++p)
      {
        AssertDimension (scalar_potential.solution_values[p].size(), 2);
        AssertDimension (computed_quantities[p].size(), 1);

        // Get complex pressure amplitude from complex scalar potential amplitude
        std::complex<double> i(0., 1.);
        std::complex<double> phi(scalar_potential.solution_values[p](0), 
                                  scalar_potential.solution_values[p](1));
        std::complex<double> pressure = density * speed_of_sound * speed_of_sound
                                        * xi * phi / (parameters.angular_frequency * i);

        // Compute time-harmonic pressure p(t) = P e^{-iwt}
        std::complex<double> eiwt = std::polar(1., - parameters.angular_frequency * current_time);
        std::complex<double> time_harmonic_pressure = pressure * eiwt;
        computed_quantities[p](0) = time_harmonic_pressure.real(); // only store real valued (physical) pressure
      }
  }

  void set_time(const double time)
  {
    current_time = time;
  }

private:
  const Settings &parameters;
  double current_time;
};

template <int dim>
class ComputeVelocityIntensity : public DataPostprocessorVector<dim>
{
public:
  ComputeVelocityIntensity(const Settings &parameters)
    : DataPostprocessorVector<dim>("velocity_intensity", 
                                    update_gradients)
    , parameters(parameters)
  {}

  virtual void
  evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &scalar_potential,
                        std::vector<Vector<double>> &computed_quantities) const override
  {
    AssertDimension (scalar_potential.solution_gradients.size(), computed_quantities.size());
    for (unsigned int p = 0; p < scalar_potential.solution_gradients.size(); ++p)
      {
        AssertDimension (scalar_potential.solution_gradients[p].size(), 2); 
        AssertDimension (computed_quantities[p].size(), dim);
        for (unsigned int d = 0; d < dim; ++d)
          {
            // Get complex velocity amplitude from complex scalar potential amplitude
            std::complex<double> velocity_d(-scalar_potential.solution_gradients[p][0][d], 
                                            -scalar_potential.solution_gradients[p][1][d]);

            // Compute the intensity of the velocity (magnitude of the complex velocity)
            computed_quantities[p](d) = std::abs(velocity_d);
          }
      }
  }

private:
  const Settings &parameters;
};

template <int dim>
class ComputeVelocityPhase : public DataPostprocessorVector<dim>
{
public:
  ComputeVelocityPhase(const Settings &parameters)
    : DataPostprocessorVector<dim>("velocity_phase", 
                                    update_gradients)
    , parameters(parameters)
  {}

  virtual void
  evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &scalar_potential,
                        std::vector<Vector<double>> &computed_quantities) const override
  {
    AssertDimension (scalar_potential.solution_gradients.size(), computed_quantities.size());
    for (unsigned int p = 0; p < scalar_potential.solution_gradients.size(); ++p)
      {
        AssertDimension (scalar_potential.solution_gradients[p].size(), 2); 
        AssertDimension (computed_quantities[p].size(), dim);
        
        for (unsigned int d = 0; d < dim; ++d)
          {
            // Get complex velocity amplitude from complex scalar potential amplitude
            std::complex<double> velocity_d(-scalar_potential.solution_gradients[p][0][d], 
                                            -scalar_potential.solution_gradients[p][1][d]);

            // Compute the phase of the velocity (argument of the complex velocity)
            computed_quantities[p](d) = std::arg(velocity_d);
          }
      }
  }

private:
  const Settings &parameters;
};

template <int dim>
class ComputeTimeHarmonicVelocity : public DataPostprocessorVector<dim>
{
public:
  ComputeTimeHarmonicVelocity(const Settings &parameters)
    : DataPostprocessorVector<dim>("time_harmonic_velocity", 
                                    update_gradients)
    , parameters(parameters)
  {}

  virtual void
  evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &scalar_potential,
                        std::vector<Vector<double>> &computed_quantities) const override
  {
    AssertDimension (scalar_potential.solution_gradients.size(), computed_quantities.size());
    for (unsigned int p = 0; p < scalar_potential.solution_gradients.size(); ++p)
      {
        AssertDimension (scalar_potential.solution_gradients[p].size(), 2); 
        AssertDimension (computed_quantities[p].size(), dim);
        
        for (unsigned int d = 0; d < dim; ++d)
          {
            // Get complex velocity amplitude from complex scalar potential amplitude 
            std::complex<double> velocity_d(-scalar_potential.solution_gradients[p][0][d],
                                            -scalar_potential.solution_gradients[p][1][d]);

            // Compute time-harmonic velocity v(t) = V e^{-iwt}
            std::complex<double> eiwt = std::polar(1., - parameters.angular_frequency * current_time);
            std::complex<double> time_harmonic_velocity = velocity_d * eiwt;
            computed_quantities[p](d) = time_harmonic_velocity.real(); // only store real valued (physical) velocity

          }
      }
  }

  void set_time(const double time)
  {
    current_time = time;
  }

private:
  const Settings &parameters;
  double current_time;
};

template <int dim>
class ComputeFluidDisplacement : public DataPostprocessorVector<dim>
{
public:
  ComputeFluidDisplacement(const Settings &parameters)
    : DataPostprocessorVector<dim>("fluid_displacement",
                                    update_gradients)
    , parameters(parameters)
  {}

  virtual void
  evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &scalar_potential,
                        std::vector<Vector<double>> &computed_quantities) const override
  {
    AssertDimension (scalar_potential.solution_gradients.size(), computed_quantities.size());
    for (unsigned int p = 0; p < scalar_potential.solution_gradients.size(); ++p)
      {
        AssertDimension (scalar_potential.solution_gradients[p].size(), 2); 
        AssertDimension (computed_quantities[p].size(), dim);
        
        for (unsigned int d = 0; d < dim; ++d)
          {
            // Get complex fluid displacement amplitude from complex scalar potential amplitude
            std::complex<double> velocity_d(-scalar_potential.solution_gradients[p][0][d],
                                            -scalar_potential.solution_gradients[p][1][d]);
            std::complex<double> i = std::complex<double>(0., 1.);
            std::complex<double> eiwt = std::polar(1., - parameters.angular_frequency * current_time);
            // std::complex<double> eiwt_t0 = std::polar(1., - parameters.angular_frequency * 0.);

            // Compute time-harmonic displacement s(t) = int_t0^t u(t') dt = U int_t0^t e^{-iwt} dt 
            std::complex<double> time_harmonic_displacement = - velocity_d / (i * parameters.angular_frequency) * (eiwt); // - eiwt_t0);
            computed_quantities[p](d) = time_harmonic_displacement.real(); // only store real valued (physical) displacement
          }
      }
  }

  void set_time(const double time)
  {
    current_time = time;
  }

private:
  const Settings &parameters;
  double current_time;
};

template <int dim>
class ComputeAcousticRadiationForce : public DataPostprocessorVector<dim>
{
public:
  ComputeAcousticRadiationForce()
    : DataPostprocessorVector<dim>("acoustic_radiation_force", 
                                    update_gradients)
  {}

  virtual void
  evaluate_scalar_field(const DataPostprocessorInputs::Scalar<dim> &acoustic_energy,
                        std::vector<Vector<double>> &computed_quantities) const override
  {
    AssertDimension (acoustic_energy.solution_gradients.size(), computed_quantities.size());
    for (unsigned int p = 0; p < acoustic_energy.solution_gradients.size(); ++p)
    {
      AssertDimension (computed_quantities[p].size(), dim);
      for (unsigned int d = 0; d < dim; ++d)
      {
        // get force as derivative of potencial energy field (Only valid in the standing wave region)
        computed_quantities[p](d) = - acoustic_energy.solution_gradients[p][d];
      }
    }
  }
};

template<int dim>
class ComputeSteadyVelocityField : public DataPostprocessorVector<dim>
{
public:
  ComputeSteadyVelocityField(const Settings &parameters)
    : DataPostprocessorVector<dim>("steady_velocity_field", 
                                    update_gradients)
    , parameters(parameters)
  {}

  virtual void
  evaluate_scalar_field(const DataPostprocessorInputs::Scalar<dim> &acoustic_energy,
                        std::vector<Vector<double>> &computed_quantities) const override
  {
    AssertDimension (acoustic_energy.solution_gradients.size(), computed_quantities.size());
    for (unsigned int p = 0; p < acoustic_energy.solution_gradients.size(); ++p)
    {
      AssertDimension (computed_quantities[p].size(), dim);

      Tensor<1, dim> velocity, initial_velocity;
      Tensor<1, dim> acoustic_force;
      for (unsigned int d = 0; d < dim; ++d)
      {
        acoustic_force[d] = - acoustic_energy.solution_gradients[p][d];
      }

      // Initial guess for velocity using linear stockes drag 
      initial_velocity = acoustic_force / (6. * M_PI * parameters.dynamic_viscosity * parameters.particle_radius);
      velocity = newton_solve(initial_velocity, acoustic_force);

      AssertDimension (computed_quantities[p].size(), dim);
      for (unsigned int d = 0; d < dim; ++d)
      {
        // output the resulting velocity
        computed_quantities[p](d) = velocity[d];
      }

    }
  }

private:

  double
  compute_drag_coefficient(const double reynolds) const
  {
    // Compute drag coefficient based on Reynolds number and clift et al. formula
    double drag_coefficient = 24. / reynolds * (1. + 0.15 * std::pow(reynolds, 0.687))
                            + 0.42 / (1. + 42500. * std::pow(reynolds, - 1.16));
    return drag_coefficient;
  }

  Tensor<1, dim>
  compute_drag_force(const Tensor<1, dim> &velocity) const
  {
    double projected_area = 0.5 * M_PI * parameters.particle_radius * parameters.particle_radius;
    double reynolds = (parameters.density * velocity.norm() * 2 * parameters.particle_radius) 
                    / parameters.dynamic_viscosity;
    double drag_coefficient = compute_drag_coefficient(reynolds);

    Tensor<1, dim> drag_force = - 0.5 * parameters.density * drag_coefficient * projected_area 
                                * velocity.norm() * velocity;
    return drag_force;
  }
  
  Vector<double>
  compute_residual(const Tensor<1, dim> &velocity, 
                    const Tensor<1, dim> &acoustic_force) const
  {
    Vector<double> residual(dim);
    Tensor<1, dim> drag_force = compute_drag_force(velocity);
    for (unsigned int d = 0; d < dim; ++d)
      residual[d] = drag_force[d] + acoustic_force[d];
    return residual;
  }

  FullMatrix<double>
  compute_numerical_jacobian(const Tensor<1, dim> &velocity,
                              const Tensor<1, dim> &acoustic_force,
                              const double epsilon = 1e-6) const
  {
    FullMatrix<double> jacobian(dim, dim);
    Vector<double> residual = compute_residual(velocity, acoustic_force);
    for (unsigned int d = 0; d < dim; ++d)
    {
      Tensor<1, dim> perturbed_velocity = velocity;
      perturbed_velocity[d] += epsilon;
      Vector<double> perturbed_residual = compute_residual(perturbed_velocity, acoustic_force);
      for (unsigned int d2 = 0; d2 < dim; ++d2)
        jacobian[d][d2] = (perturbed_residual[d2] - residual[d2]) / epsilon;
    }
    return jacobian;
  }

  Tensor<1, dim>
  newton_solve(const Tensor<1, dim> &initial_velocity, 
               const Tensor<1, dim> &acoustic_force,
               const double tolerance = 1e-6,
               const unsigned int max_iter = 1000) const
  {
    Tensor<1, dim> velocity = initial_velocity;
    Vector<double> delta_velocity(dim);

    for (unsigned int iter = 0; iter < max_iter; ++iter)
    {
      // Compute the residual
      Vector<double> residual = compute_residual(velocity, acoustic_force);
      if (residual.l2_norm() < tolerance)
        break;

      // Compute the Jacobian matrix using finite differences
      FullMatrix<double> jacobian = compute_numerical_jacobian(velocity, acoustic_force);
      
      // Solve the linear system J dx = -r
      residual *= -1.;
      jacobian.gauss_jordan();
      jacobian.vmult(delta_velocity, residual);

      for (unsigned int d = 0; d < dim; ++d)
        velocity[d] += delta_velocity[d];

      // If max iteration reached, print warning
      if (iter == max_iter - 1)
        std::cout << "Warning: Newton solver did not converge after " 
                  << max_iter << " iterations." << std::endl;
    }

    return velocity;
  }

  const Settings parameters;
};



template <int dim>
class MMSAnalyticalSolution : public Function<dim>
{
public:
  MMSAnalyticalSolution(const Settings &parameters)
    : Function<dim>(2) 
    , parameters(parameters)
  {}
  
  virtual double
  value(const Point<dim> &p,
        const unsigned int component) const override
  {
    double return_value = 0.;
    if (parameters.problem_type == Settings::mms)
    {
      std::complex<double> xi = parameters.look_up_xi[0];

      switch (parameters.mms_problem)
      {
        case Settings::sine_wave:
        {
          double k = std::sqrt(xi.real() / dim);
          std::complex<double> sine_wave(1., 1.); // The real and imaginary parts will be the same
          for (unsigned int d1 = 0; d1 < dim; ++d1)
            sine_wave *= std::sin(k * p(d1));

          if (component == 0) {return_value = sine_wave.real();}
          else if (component == 1) {return_value = sine_wave.imag();}
          break;
        }

        case Settings::plane_wave:
        {
          Tensor<1, dim> k;
          for (unsigned int d1 = 0; d1 < dim; ++d1)
            k[d1] = std::sqrt(xi.real() / dim);

          double wave_phase = k * p;
          std::complex<double> plane_wave = std::polar(1., wave_phase);
          if (component == 0) {return_value = plane_wave.real();}
          else if (component == 1) {return_value = plane_wave.imag();}
          break;
        }
      }
    }

    return return_value;
  }

  virtual Tensor<1, dim>
  gradient(const Point<dim> &p,
            const unsigned int component) const override
  {
    Tensor<1, dim> return_value;
    if (parameters.problem_type == Settings::mms)
    {
      std::complex<double> xi = parameters.look_up_xi[0];

      switch (parameters.mms_problem)
      {
        case Settings::sine_wave:
        {
          double k = std::sqrt(xi.real() / dim);
          for (unsigned int d1 = 0; d1 < dim; ++d1)   
          {
            std::complex<double> sine_wave_grad_d1(1., 1.);
            sine_wave_grad_d1 *= k * std::cos(k * p(d1));
            for (unsigned int d2 = 0; d2 < dim; ++d2)
              if (d2 != d1) {sine_wave_grad_d1 *= std::sin(k * p(d2));}
              
            if (component == 0) {return_value[d1] = sine_wave_grad_d1.real();}
            else if (component == 1) {return_value[d1] = sine_wave_grad_d1.imag();}
          }
          break;
        }

        case Settings::plane_wave:
        {
          std::complex<double> i(0., 1.);
          std::complex<double> plane_wave(value(p, 0), value(p, 1));
          Tensor<1, dim> k;
          for (unsigned int d1 = 0; d1 < dim; ++d1)
            k[d1] = std::sqrt(xi.real() / dim);

          for (unsigned int d1 = 0; d1 < dim; ++d1)  
          {
            if (component == 0) {return_value[d1] = (i * k[d1] * plane_wave).real();}
            else if (component == 1) {return_value[d1] = (i * k[d1] * plane_wave).imag();}
          }
          break;
        }
      }
    }
    return return_value;
  }

private:
  const Settings &parameters;
};

template <int dim>
class PMLStretching : public TensorFunction<2, dim>
{
public :
  PMLStretching(const Settings &parameters, bool real_stretching = true)
    : TensorFunction<2, dim>()
    , parameters(parameters)
    , real_stretching(real_stretching)
  {}

  virtual Tensor<2, dim>
  value(const Point<dim> &/*p*/) const override
  {
    Tensor<2, dim> return_value;
    Tensor<1, dim> pml_abs_distance;

    // Compute the distance within the PML region
    switch (parameters.problem_type)
    {
      case Settings::mms: 
        break;

      case Settings::channel:
        break;
    }

    Tensor<2, dim> identity = unit_symmetric_tensor<dim>();
    Tensor<2, dim> activation_tensor;
    for (unsigned int d = 0; d < dim; ++d)
    {
      double pml_thickness = parameters.pml_thickness_ratio * parameters.wavelength;
      activation_tensor[d][d] = std::pow(pml_abs_distance[d] / pml_thickness, parameters.pml_exponent);
    }

    if (real_stretching) {return_value = identity + (parameters.pml_stretching - 1.) * activation_tensor;}
    else {return_value = parameters.pml_attenuation / parameters.angular_frequency * activation_tensor;}
    return return_value;
  }

private :
  const Settings    &parameters;
  const bool        real_stretching;
};  
 
template <int dim>
class HelmholtzProblem
{
public:
  HelmholtzProblem(const Settings &parameters);

  virtual void 
  run();

protected :

  void 
  make_grid();

  void 
  setup_system();

  void 
  assemble_system();

  void 
  solve();

  void
  post_process_acoustic_potential();

  void
  post_process_convergence(unsigned int cycle);

  void
  post_process_quantities();

  virtual void
  output_results() const;

private :
  Triangulation<dim>    triangulation;
  const FESystem<dim>   fe;
  DoFHandler<dim>       dof_handler;

  Vector<double>        solution;

  SparsityPattern       sparsity_pattern;
  SparseMatrix<double>  system_matrix;

  Vector<double>        system_rhs;

  Settings              parameters;

  double  iteration_end;

  ConvergenceTable      convergence_table;
};

template <int dim>
HelmholtzProblem<dim>::HelmholtzProblem(const Settings &parameters)
  : fe(FE_Q<dim>(parameters.fe_order)^2)
  , dof_handler(triangulation)
  , parameters(parameters)
{}

// @brief Make the grid for the system based on the problem type. 
// Boundary ids are as follows:
//
// Dirichlet boundaries
// 10: MMS Dirichlet boundary
// 11: Sound soft boundary
//
// Neumann boundaries
// 20: MMS Neumann boundary
// 21: Sound hard boundary
// 22: Transducer boundary
//
// Robin boundaries
// 30: MMS Robin boundary
// 31: Air Impedance boundary
// 32: Impedance matching boundary

template <int dim>
void 
HelmholtzProblem<dim>::make_grid()
{
  switch (parameters.problem_type)
  {
    case Settings::mms:
    {
      switch (parameters.mms_geometry)
      {
        case Settings::cube:
        {
          // Create a cube of sides 1.0 centered at the origin
          GridGenerator::hyper_cube(triangulation, -0.5, 0.5, true);
          break;
        }
      
        case Settings::ball:
        {
          // Create a ball of diameter 1. centered at the origin
          GridGenerator::hyper_ball(triangulation, Point<dim>(), 0.5);
          break;
        }
      } 

      // In the case of the cube geometry, only set the x+ and x- boundaries (ids 0 and 1)
      // to the desired boundary condition type. In the case of the ball geometry, all boundaries
      // are set to the desired boundary condition type.
      for (const auto &cell : triangulation.active_cell_iterators())
      for (const auto &face : cell->face_iterators())
      { 
        if (face->at_boundary() && (face->boundary_id() <= 1))
        {
          if (parameters.mms_boundary_type == Settings::dirichlet)
            face->set_boundary_id(10);
          else if (parameters.mms_boundary_type == Settings::neumann)
            face->set_boundary_id(20);
          else if (parameters.mms_boundary_type == Settings::robin)
            face->set_boundary_id(30);
        }
        // Set the rest of boundaries as dirichlet
        else if (face->at_boundary())
          face->set_boundary_id(10);
      }
      break;
    }
    break;
    
    case Settings::channel:
    {
      // Domain of size L x H x H from -L/2 to L/2 in x and -H/2 to H/2 in y and z
      // Transducer from -L_t/2 to L_t/2 in x and -H/2 to H/2 in z
      // PML layer for x < -L_PML and x > L_PML
      Point<dim> p1, p2;
      std::vector<unsigned int> repetitions;
      double x_cells_domain = parameters.x_cell_ratio * parameters.domain_length_ratio;
      double x_cells_transducer = parameters.x_cell_ratio * parameters.transducer_length_ratio;

      double y_cells_fluid = parameters.y_cell_ratio * 0.5;
      double y_cells_wall_thin = parameters.y_cell_ratio * parameters.wall_thickness_ratio;
      double y_cells_wall_thick = parameters.y_cell_ratio * (0.25 - parameters.wall_thickness_ratio);
      
      double fluid_length = parameters.domain_length_ratio * parameters.wavelength;
      double transducer_length = parameters.transducer_length_ratio * parameters.wavelength;

      double fluid_height = 0.5 * parameters.wavelength;
      double wall_thickness_thin = parameters.wall_thickness_ratio * parameters.wall_wavelength;
      double wall_thickness_thick = 0.25 * parameters.wall_wavelength - wall_thickness_thin;

      // Create the fluid domain
      p1 = Point<dim>(-0.5*fluid_length, -0.5*fluid_height);
      p2 = Point<dim>(0.5*fluid_length, 0.5*fluid_height);
      repetitions = {(unsigned int)x_cells_domain, (unsigned int)y_cells_fluid};
      Triangulation<dim> fluid_triangulation;
      GridGenerator::subdivided_hyper_rectangle(fluid_triangulation, repetitions, p1, p2, true);
      for (const auto &cell : fluid_triangulation.active_cell_iterators())
      {
        cell->set_material_id(0); // Set fluid domain to material id 0
        for (const auto &face : cell->face_iterators())
        {
          // Set front and back faces to impedance matching boundary and the rest to nothing
          if (face->at_boundary() && ((face->boundary_id() == 0) || (face->boundary_id() == 1)))
            face->set_boundary_id(32);
          else if (face->at_boundary())
            face->set_boundary_id(0);
        }
      }

      // Create the wall domain top
      p1 = Point<dim>(-0.5*fluid_length, 0.5*fluid_height);
      p2 = Point<dim>(0.5*fluid_length, 0.5*fluid_height + wall_thickness_thin);
      repetitions = {(unsigned int)x_cells_domain, (unsigned int)y_cells_wall_thin};
      Triangulation<dim> wall_top_triangulation;
      GridGenerator::subdivided_hyper_rectangle(wall_top_triangulation, repetitions, p1, p2, true);
      for (const auto &cell : wall_top_triangulation.active_cell_iterators())
      {
        cell->set_material_id(1); // Set wall domain to material id 1
        for (const auto &face : cell->face_iterators())
        {
          // Set front and back faces to impedance matching boundary 
          if (face->at_boundary() && ((face->boundary_id() == 0) || (face->boundary_id() == 1)))
            face->set_boundary_id(32);
          // Set top boundary to impedance if outside the transducer
          else if (face->at_boundary() && (face->boundary_id() == 3) && 
                   ((face->center()[0] < -0.5*transducer_length) || 
                    (face->center()[0] > 0.5*transducer_length)))
            face->set_boundary_id(31);
            // face->set_boundary_id(0); // Sound Hard
          else if (face->at_boundary())
            face->set_boundary_id(0);
        }
      }

      // Create the wall domain bottom
      p1 = Point<dim>(-0.5*fluid_length, -0.5*fluid_height - wall_thickness_thin);
      p2 = Point<dim>(0.5*fluid_length, -0.5*fluid_height);
      repetitions = {(unsigned int)x_cells_domain, (unsigned int)y_cells_wall_thin};
      Triangulation<dim> wall_bottom_triangulation;
      GridGenerator::subdivided_hyper_rectangle(wall_bottom_triangulation, repetitions, p1, p2, true);
      for (const auto &cell : wall_bottom_triangulation.active_cell_iterators())
      {
        cell->set_material_id(1); // Set wall domain to material id 1
        for (const auto &face : cell->face_iterators())
        {
          // Set front and back faces to impedance matching boundary 
          if (face->at_boundary() && ((face->boundary_id() == 0) || (face->boundary_id() == 1)))
            face->set_boundary_id(32);
          // Set bottom boundary to impedance if outside the transducer
          else if (face->at_boundary() && (face->boundary_id() == 2) && 
                   ((face->center()[0] < -0.5*transducer_length) || 
                    (face->center()[0] > 0.5*transducer_length)))
            face->set_boundary_id(31);
            // face->set_boundary_id(0); // Sound Hard
          else if (face->at_boundary())
            face->set_boundary_id(0);
        }
      }

      // Create the thick wall domain top
      p1 = Point<dim>(-0.5*transducer_length, 0.5*fluid_height + wall_thickness_thin);
      p2 = Point<dim>(0.5*transducer_length, 0.5*fluid_height + wall_thickness_thin + wall_thickness_thick);
      repetitions = {(unsigned int)x_cells_transducer, (unsigned int)y_cells_wall_thick};
      Triangulation<dim> wall_top_thick_triangulation;
      GridGenerator::subdivided_hyper_rectangle(wall_top_thick_triangulation, repetitions, p1, p2, true);
      for (const auto &cell : wall_top_thick_triangulation.active_cell_iterators())
      {
        cell->set_material_id(1); // Set wall domain to material id 1
        for (const auto &face : cell->face_iterators())
        {
          // Set front and back faces to impedance 
          if (face->at_boundary() && ((face->boundary_id() == 0) || (face->boundary_id() == 1)))
            face->set_boundary_id(31);
            // face->set_boundary_id(0); // Sound Hard
          // Set top boundary to transducer
          else if (face->at_boundary() && (face->boundary_id() == 3))
            face->set_boundary_id(22);
          else if (face->at_boundary())
            face->set_boundary_id(0);
        }
      }
      
      // Create the thick wall domain bottom
      p1 = Point<dim>(-0.5*transducer_length, -0.5*fluid_height - wall_thickness_thin - wall_thickness_thick);
      p2 = Point<dim>(0.5*transducer_length, -0.5*fluid_height - wall_thickness_thin);
      repetitions = {(unsigned int)x_cells_transducer, (unsigned int)y_cells_wall_thick};
      Triangulation<dim> wall_bottom_thick_triangulation;
      GridGenerator::subdivided_hyper_rectangle(wall_bottom_thick_triangulation, repetitions, p1, p2, true);
      for (const auto &cell : wall_bottom_thick_triangulation.active_cell_iterators())
      {
        cell->set_material_id(1); // Set wall domain to material id 0
        for (const auto &face : cell->face_iterators())
        {
          // Set front and back faces to impedance 
          if (face->at_boundary() && ((face->boundary_id() == 0) || (face->boundary_id() == 1)))
            face->set_boundary_id(31);
            // face->set_boundary_id(0); // Sound Hard
          // Set bottom boundary to transducer
          else if (face->at_boundary() && (face->boundary_id() == 2))
            face->set_boundary_id(23);
          else if (face->at_boundary())
            face->set_boundary_id(0);
        }
      }

      // Merge Triangulations    
      triangulation.clear();
      dealii::GridGenerator::merge_triangulations({&fluid_triangulation, 
                                                    &wall_top_triangulation, 
                                                    &wall_bottom_triangulation,
                                                    &wall_top_thick_triangulation,
                                                    &wall_bottom_thick_triangulation}, 
                                                  triangulation,
                                                  1e-12,
                                                  false,
                                                  true
                                                  );

      break;
    }

  //   case Settings::channel:
  //     {
  //     // Domain of size L x H x H from -L/2 to L/2 in x and -H/2 to H/2 in y and z
  //     // Transducer from -L_t/2 to L_t/2 in x and -H/2 to H/2 in z
  //     // PML layer for x < -L_PML and x > L_PML 
  //     Point<dim> p1, p2;
  //     std::vector<unsigned int> repetitions;
  //     double x_cells = parameters.x_cell_ratio * parameters.domain_length_ratio;
  //     double y_cells = parameters.y_cell_ratio * parameters.domain_height_ratio;

  //     double domain_length = parameters.domain_length_ratio * parameters.wavelength;
  //     double domain_height = parameters.domain_height_ratio * parameters.wavelength;
  //     double transducer_length = parameters.transducer_length_ratio * parameters.wavelength;

  //     switch (parameters.dimensions)
  //     {
  //       case 2: 
  //       {
  //         p1 = Point<dim>(-0.5*domain_length, -0.5*domain_height);
  //         p2 = Point<dim>(0.5*domain_length, 0.5*domain_height);
  //         repetitions = {(unsigned int)x_cells, (unsigned int)y_cells};
  //         break;
  //       }
  //       case 3:
  //       {
  //         p1 = Point<dim>(-0.5*domain_length, -0.5*domain_height, -0.5*domain_height);
  //         p2 = Point<dim>(0.5*domain_length, 0.5*domain_height, 0.5*domain_height);
  //         repetitions = {(unsigned int)x_cells, (unsigned int)y_cells, (unsigned int)y_cells};
  //         break;
  //       }
  //       default:
  //         AssertThrow(false, ExcNotImplemented());
  //     }

  //     GridGenerator::subdivided_hyper_rectangle(triangulation, repetitions, p1, p2, true);
  //     for (const auto &cell : triangulation.active_cell_iterators())
  //       {
  //         // Set boundary ids for the domain
  //         for (const auto &face : cell->face_iterators())
  //           {
  //             if (face->at_boundary())
  //             {
  //               if ((face->boundary_id() == 0) || (face->boundary_id() == 1))
  //               {
  //                 face->set_boundary_id(32);
  //               }

  //               // if face is within - transducer /2 <x< transducer /2
  //               else if ((face->center()[0] >= -0.5 * transducer_length)
  //                     && (face->center()[0] <= 0.5 * transducer_length))
  //               {
  //                 if ((face->boundary_id() == 2) )//|| (face->boundary_id() == 3))
  //                 {
  //                   face->set_boundary_id(22);
  //                 }
  //                 else if (face->boundary_id() == 3)
  //                 {
  //                   face->set_boundary_id(23); // sound hard
  //                 }
  //                 else
  //                 {
  //                   face->set_boundary_id(31);
  //                 }
  //               }
  //               else
  //               {
  //                 face->set_boundary_id(31);
  //               }
  //             }
  //           } 
  //       }
  //     break;
  //   }
  }

  // Refine the domain based on the user input
  triangulation.refine_global(parameters.initial_refinement);

  // Output domain to file for debugging
  GridOut grid_out;
  std::ofstream out(parameters.output_path + parameters.output_name + "_grid.vtk");
  grid_out.write_vtk(triangulation, out);
  out.close();
}

template <int dim>
void 
HelmholtzProblem<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);
  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}

template <int dim>
void 
HelmholtzProblem<dim>::assemble_system()
{
  const QGauss<dim>     quadrature_formula(fe.degree + 1);
  const QGauss<dim - 1> face_quadrature_formula(fe.degree + 1);

  const unsigned int n_q_points       = quadrature_formula.size();
  const unsigned int n_face_q_points  = face_quadrature_formula.size();
  const unsigned int dofs_per_cell    = fe.n_dofs_per_cell();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  FEValues<dim>     fe_values(fe, 
                      quadrature_formula, 
                      update_values | update_gradients | 
                      update_quadrature_points | update_JxW_values);
  FEFaceValues<dim> fe_face_values(fe,
                      face_quadrature_formula,
                      update_values | update_quadrature_points |
                      update_normal_vectors | update_JxW_values);

  // Use fe extractors to get the real and imaginary parts of the solution
  const FEValuesExtractors::Scalar real(0);
  const FEValuesExtractors::Scalar imag(1);

  // Create vectors to store the test function values and gradients for each dof
  std::vector<double>         psi_real(dofs_per_cell),      psi_imag(dofs_per_cell);
  std::vector<Tensor<1, dim>> grad_psi_real(dofs_per_cell), grad_psi_imag(dofs_per_cell);

  // Store the PML stretching tensor (real and imaginary parts)
  PMLStretching<dim> pml_stretching_real(parameters, true);
  PMLStretching<dim> pml_stretching_imag(parameters, false);
  std::vector<Tensor<2, dim>> Lambda_real(n_q_points), Lambda_imag(n_q_points);

  // constants for the Helmholtz equation
  std::complex<double> i_complex(0., 1.);

  MMSAnalyticalSolution<dim> mms_analytical_solution(parameters);

  // Iterate over all cells
  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    fe_values.reinit(cell);
    cell_matrix = 0.;
    cell_rhs    = 0.;

    double density = parameters.look_up_density[cell->material_id()];
    double speed_of_sound = parameters.look_up_speed_of_sound[cell->material_id()];
    std::complex<double> xi = parameters.look_up_xi[cell->material_id()];

    pml_stretching_real.value_list(fe_values.get_quadrature_points(),
                                      Lambda_real);
    pml_stretching_imag.value_list(fe_values.get_quadrature_points(),
                                      Lambda_imag);

    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
    {
      Tensor<2, dim, std::complex<double>> Lambda_q = Lambda_real[q_point] 
                                                      + i_complex * Lambda_imag[q_point];
      Tensor<2, dim, std::complex<double>> LL_ij = transpose(invert(Lambda_q)) * invert(Lambda_q);
      std::complex<double> LL_det = determinant(LL_ij);

      Tensor<2, dim> A_ij_real, A_ij_imag;
      for (unsigned int d1 = 0; d1 < dim; ++d1)
        for (unsigned int d2 = 0; d2 < dim; ++d2)
        {
          A_ij_real[d1][d2] = (LL_ij[d1][d2] * LL_det).real();
          A_ij_imag[d1][d2] = (LL_ij[d1][d2] * LL_det).imag();
        }
        
      // Store the test function for real and imaginary pats for ease of use
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        psi_real[i] = fe_values[real].value(i, q_point);
        psi_imag[i] = fe_values[imag].value(i, q_point);
        grad_psi_real[i] = fe_values[real].gradient(i, q_point);
        grad_psi_imag[i] = fe_values[imag].gradient(i, q_point);
      }

      // Assemble the cell matrix for the domain integral
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
        { 
          // Real equation on cell domain
          cell_matrix(i, j) += (
            (xi * LL_det).real() * psi_real[i] * psi_real[j] 
            // - grad_psi_real[i] * A_ij_real * grad_psi_real[j]
            - (xi * LL_det).imag() * psi_real[i] * psi_imag[j]
            // + grad_psi_real[i] * A_ij_imag * grad_psi_imag[j]
            ) * fe_values.JxW(q_point);
          
          // Imaginary equation on cell domain
          cell_matrix(i, j) += (
            (xi * LL_det).imag() * psi_imag[i] * psi_real[j]
            // - grad_psi_imag[i] * A_ij_imag * grad_psi_real[j]
            + (xi * LL_det).real() * psi_imag[i] * psi_imag[j]
            // - grad_psi_imag[i] * A_ij_real * grad_psi_imag[j]
            ) * fe_values.JxW(q_point);

            for (unsigned int d1 = 0; d1 < dim; ++d1)
            {
              for (unsigned int d2 = 0; d2 < dim; ++d2)
              {
                cell_matrix(i, j) += (
                  - A_ij_real[d1][d2] * grad_psi_real[i][d1] * grad_psi_real[j][d2]
                  + A_ij_imag[d1][d2] * grad_psi_real[i][d1] * grad_psi_imag[j][d2]
                ) * fe_values.JxW(q_point);
    
                cell_matrix(i, j) += (
                  - A_ij_imag[d1][d2] * grad_psi_imag[i][d1] * grad_psi_real[j][d2]
                  - A_ij_real[d1][d2] * grad_psi_imag[i][d1] * grad_psi_imag[j][d2]
                ) * fe_values.JxW(q_point);
              }
            }
        }
      }
    }

    for (const auto &face : cell->face_iterators())
      if ((face->at_boundary()) && (face->boundary_id() >= 20))
        {
          fe_face_values.reinit(cell, face);
          
          for (unsigned int q_point = 0; q_point < n_face_q_points; ++q_point)
            { 
              
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  psi_real[i] = fe_face_values[real].value(i, q_point);
                  psi_imag[i] = fe_face_values[imag].value(i, q_point);
                }
              
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  std::complex<double> robin_coeff(0.0, 0.0);
                  std::complex<double> neumann_coeff(0.0, 0.0);

                  if (face->boundary_id() == 20) // MMS Neumann boundary
                  {
                    neumann_coeff = (
                      mms_analytical_solution.gradient(fe_face_values.quadrature_point(q_point), 0) 
                      + i_complex * mms_analytical_solution.gradient(fe_face_values.quadrature_point(q_point), 1)
                      ) * fe_face_values.normal_vector(q_point);
                  }
                  // Sound hard boundary => neumann = 0 so nothing to do
                  else if (face->boundary_id() == 22)  // Transducer boundary top
                  {
                    std::complex<double> transducer_displacement(parameters.transducer_amplitude, 
                                                                  parameters.transducer_phase);

                    neumann_coeff = i_complex * parameters.angular_frequency * transducer_displacement;
                  }
                  else if (face->boundary_id() == 23)  // Transducer boundary top
                  {
                    std::complex<double> transducer_displacement(parameters.transducer_amplitude, 
                                                                  parameters.transducer_phase);

                    neumann_coeff = - i_complex * parameters.angular_frequency * transducer_displacement;
                  }
                  else if (face->boundary_id() == 30)  // MMS Robin boundary
                  {
                    std::complex<double> normal_gradient = (
                      mms_analytical_solution.gradient(fe_face_values.quadrature_point(q_point), 0) 
                      + i_complex * mms_analytical_solution.gradient(fe_face_values.quadrature_point(q_point), 1)
                      ) * fe_face_values.normal_vector(q_point);

                    std::complex<double> exact_value = 
                      mms_analytical_solution.value(fe_face_values.quadrature_point(q_point), 0)
                      + i_complex * mms_analytical_solution.value(fe_face_values.quadrature_point(q_point), 1);

                    robin_coeff = normal_gradient / exact_value;
                  }
                  else if (face->boundary_id() == 31)  // Air Impedance boundary
                  {
                    double air_density = 1.21; // kg/m^3
                    double air_speed_of_sound = 343.0; // m/s
                    double air_impedance = air_density * air_speed_of_sound;

                    robin_coeff = - density * speed_of_sound 
                      * speed_of_sound * xi / (i_complex 
                      * parameters.angular_frequency * air_impedance);
                  }
                  else if (face->boundary_id() == 32) // Impedance matching boundary
                  {
                    robin_coeff = - speed_of_sound
                      * xi / (i_complex * parameters.angular_frequency);
                  }
                  
                  // Robin contribution
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  {
                    // Real equation
                    cell_matrix(i, j) += (robin_coeff.real() * psi_real[i] * psi_real[j]
                                        - robin_coeff.imag() * psi_real[i] * psi_imag[j]
                                        ) * fe_face_values.JxW(q_point);
                    
                    // Imaginary equation
                    cell_matrix(i, j) += (robin_coeff.imag() * psi_imag[i] * psi_real[j]
                                        + robin_coeff.real() * psi_imag[i] * psi_imag[j]
                                        ) * fe_face_values.JxW(q_point);
                  }
                  
                  // Neumann contribution for real and imaginary equations
                  cell_rhs(i) += - psi_real[i] * neumann_coeff.real() 
                    * fe_face_values.JxW(q_point);
                  cell_rhs(i) += - psi_imag[i] * neumann_coeff.imag() 
                    * fe_face_values.JxW(q_point);
                }
            }
        }

    cell->get_dof_indices(local_dof_indices);
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
      {
        system_matrix.add(local_dof_indices[i], 
                          local_dof_indices[j], 
                          cell_matrix(i, j));
      }
      system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }
  }
  
  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler,
                                            types::boundary_id(10),
                                            mms_analytical_solution,
                                            boundary_values);

  MatrixTools::apply_boundary_values(boundary_values,
                                      system_matrix,
                                      solution,
                                      system_rhs);
}

template <int dim>
void 
HelmholtzProblem<dim>::solve()
{ 
  SparseDirectUMFPACK A_direct;
  A_direct.initialize(system_matrix);
  A_direct.vmult(solution, system_rhs);
}

template <int dim>
void
HelmholtzProblem<dim>::post_process_acoustic_potential()
{
  // Setup global vector for the acoustic potential
  DoFHandler<dim> dof_handler_potential(triangulation);
  FE_Q<dim> fe_potential(fe.degree);  // Scalar finite element
  dof_handler_potential.distribute_dofs(fe_potential);

  SparseMatrix<double> lhs_matrix;
  Vector<double> acoustic_potential;
  Vector<double> rhs_vector;

  DynamicSparsityPattern dsp_potential(dof_handler_potential.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler_potential, dsp_potential);
  SparsityPattern       sparsity_pattern_potential;
  sparsity_pattern_potential.copy_from(dsp_potential);

  lhs_matrix.reinit(sparsity_pattern_potential);
  acoustic_potential.reinit(dof_handler_potential.n_dofs());
  rhs_vector.reinit(dof_handler_potential.n_dofs());

  // Setup FEValues fo te scalar potiential linked to "solution"
  const QGauss<dim> quadrature_formula(fe.degree + 1);
  FEValues<dim>     fe_values(fe, 
                      quadrature_formula, 
                      update_values | update_gradients | 
                      update_quadrature_points | update_JxW_values);
  FEValues<dim>     fe_values_potential(fe_potential, 
                      quadrature_formula, 
                      update_values | update_gradients | 
                      update_quadrature_points | update_JxW_values);

  const unsigned int  n_q_points    = quadrature_formula.size();
  const unsigned int  dofs_per_cell = fe_potential.n_dofs_per_cell();

  FullMatrix<double>  cell_mass_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>      cell_rhs(dofs_per_cell);
  std::vector<double> acoustic_potential_q(n_q_points);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  auto cell_1 = dof_handler.begin_active();
  auto cell_2 = dof_handler_potential.begin_active();

  AssertThrow(&dof_handler_potential.get_triangulation() ==
              &dof_handler.get_triangulation(),
              ExcMessage("DoFHandlers must use the same triangulation"));

  for (; cell_1 != dof_handler.end() && cell_2 != dof_handler_potential.end(); ++cell_1, ++cell_2)
  {
    Assert(cell_1 != dof_handler.end(), ExcInternalError());
    Assert(cell_1->id() == cell_2->id(), ExcMessage("Cells don't match!"));

    fe_values.reinit(cell_1);
    fe_values_potential.reinit(cell_2);

    cell_mass_matrix = 0.;
    cell_rhs = 0.;

    double density = parameters.look_up_density[cell_2->material_id()];
    double speed_of_sound = parameters.look_up_speed_of_sound[cell_2->material_id()];
    std::complex<double> xi = parameters.look_up_xi[cell_2->material_id()];

    std::vector<Vector<double>> phi_values(fe_values.n_quadrature_points, 
                                                Vector<double>(2));
    std::vector<std::vector<Tensor<1,dim>>> phi_gradients(fe_values.n_quadrature_points, 
                                            std::vector<Tensor<1,dim>>(2));

    fe_values.get_function_values(solution, phi_values);
    fe_values.get_function_gradients(solution, phi_gradients);

    // Now we can use solution_values at each quadrature point
    for (unsigned int q = 0; q < n_q_points; ++q)
    {
      std::complex<double> i_complex(0., 1.);
      std::complex<double> phi(phi_values[q](0), phi_values[q](1));
      std::complex<double> pressure_q = density * speed_of_sound * speed_of_sound
                                      * xi * phi 
                                      / (parameters.angular_frequency * i_complex);
      
      Tensor<1, dim, std::complex<double>> velocity_q;
      for (unsigned int d = 0; d < dim; ++d)
      {
        velocity_q[d] = - phi_gradients[q][0][d] + i_complex * phi_gradients[q][1][d];
      }

      double time_averaged_pressure_square = 
        0.5 * std::real(pressure_q * std::conj(pressure_q));

      double time_averaged_velocity_square = 0.0;
      for (unsigned int d = 0; d < dim; ++d)
        time_averaged_velocity_square += 
          0.5 * std::real(velocity_q[d] * std::conj(velocity_q[d]));

      double isothermal_compressibility = 1. / 
        (density * speed_of_sound * speed_of_sound);

      double isothermal_compressibility_particle = 1. / 
        (parameters.particle_density * parameters.particle_speed_of_sound * parameters.particle_speed_of_sound);

      double compressibility_ratio = 
        isothermal_compressibility_particle / isothermal_compressibility;
      double density_ratio = parameters.particle_density / density;

      double f_1_coeff = 1 - compressibility_ratio;
      double f_2_coeff = 2 * (density_ratio - 1) / (2 * density_ratio + 1);

      acoustic_potential_q[q] = 
        4 * M_PI * parameters.particle_radius * parameters.particle_radius * parameters.particle_radius / 3
        * (0.5 * f_1_coeff * isothermal_compressibility * time_averaged_pressure_square
        - 0.75 * f_2_coeff * density * time_averaged_velocity_square);
      
      // If the cell is a wall, set the acoustic potential to zero
      if (cell_2->material_id() == 1)
        acoustic_potential_q[q] = 0.0;

      // Build mass matrix and rhs vector
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
        {
          cell_mass_matrix(i, j) += (
            fe_values_potential.shape_value(i, q) *
            fe_values_potential.shape_value(j, q)
            ) * fe_values_potential.JxW(q);
        }
        cell_rhs(i) += fe_values_potential.shape_value(i, q) 
                      * acoustic_potential_q[q] * fe_values_potential.JxW(q);
      }
    }
  
    cell_2->get_dof_indices(local_dof_indices);
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
      {
        lhs_matrix.add(local_dof_indices[i], 
                        local_dof_indices[j], 
                        cell_mass_matrix(i, j));
      }
      rhs_vector(local_dof_indices[i]) += cell_rhs(i);
    }

  }

  SparseDirectUMFPACK A_direct;
  A_direct.initialize(lhs_matrix);
  A_direct.vmult(acoustic_potential, rhs_vector);

  // Output the acoustic potential and the triangulation data
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler_potential);
  

  data_out.add_data_vector(acoustic_potential, "acoustic_potential");

  data_out.attach_triangulation(triangulation);
  // Add material_id as cell data
  Vector<float> material_ids(triangulation.n_active_cells());
  for (const auto &cell : triangulation.active_cell_iterators())
    material_ids[cell->active_cell_index()] = cell->material_id();

  data_out.add_data_vector(material_ids, "material_id");

  ComputeAcousticRadiationForce<dim> acoustic_radiation_force;
  data_out.add_data_vector(acoustic_potential, acoustic_radiation_force);

  ComputeSteadyVelocityField<dim> steady_velocity_field(parameters);
  data_out.add_data_vector(acoustic_potential, steady_velocity_field);

  // Output the triangulation data
  data_out.attach_triangulation(triangulation);

  data_out.build_patches();

  std::ofstream output(parameters.output_path + parameters.output_name + "_force.vtk");
  data_out.write_vtk(output);
}

template <int dim>
void
HelmholtzProblem<dim>::post_process_convergence(unsigned int cycle)
{
  Vector<float> difference_per_cell(triangulation.n_active_cells());

  // Compute the L2 norm of the error
  VectorTools::integrate_difference(dof_handler,
                                    solution,
                                    MMSAnalyticalSolution<dim>(parameters),
                                    difference_per_cell,
                                    QGauss<dim>(fe.degree + 1),
                                    VectorTools::L2_norm);
  double error_l2 = VectorTools::compute_global_error(triangulation,
                                                      difference_per_cell,
                                                      VectorTools::L2_norm);

  VectorTools::integrate_difference(dof_handler,
                                    solution,
                                    MMSAnalyticalSolution<dim>(parameters),
                                    difference_per_cell,
                                    QGauss<dim>(fe.degree + 1),
                                    VectorTools::H1_seminorm);
  double error_h1 = VectorTools::compute_global_error(triangulation,
                                                      difference_per_cell,
                                                      VectorTools::H1_seminorm);

  const QTrapezoid<1>  q_trapez;
  const QIterated<dim> q_iterated(q_trapez, fe.degree * 2 + 1);
  VectorTools::integrate_difference(dof_handler,
                                    solution,
                                    MMSAnalyticalSolution<dim>(parameters),
                                    difference_per_cell,
                                    q_iterated,
                                    VectorTools::Linfty_norm);
  double error_linf = VectorTools::compute_global_error(triangulation,
                                                        difference_per_cell,
                                                        VectorTools::Linfty_norm);

  convergence_table.add_value("cells", triangulation.n_global_active_cells());
  convergence_table.add_value("dofs", dof_handler.n_dofs());
  convergence_table.add_value("L2", error_l2);
  convergence_table.add_value("H1", error_h1);
  convergence_table.add_value("Linfty", error_linf);

  if (cycle == parameters.n_cycles - 1) // on last cycle format and print the table
    {
      convergence_table.set_precision("L2", 3);
      convergence_table.set_precision("H1", 3);
      convergence_table.set_precision("Linfty", 3);

      convergence_table.set_scientific("L2", true);
      convergence_table.set_scientific("H1", true);
      convergence_table.set_scientific("Linfty", true);

      convergence_table.set_tex_caption("cells", "\\# cells");
      convergence_table.set_tex_caption("dofs", "\\# dofs");
      convergence_table.set_tex_caption("L2", "$L^2$-error");
      convergence_table.set_tex_caption("H1", "$H^1$-error");
      convergence_table.set_tex_caption("Linfty", "$L^\\infty$-error");

      convergence_table.evaluate_convergence_rates(
        "L2", ConvergenceTable::reduction_rate_log2);
      convergence_table.evaluate_convergence_rates(
        "H1", ConvergenceTable::reduction_rate_log2);
      convergence_table.evaluate_convergence_rates(
        "Linfty", ConvergenceTable::reduction_rate_log2);
      
      convergence_table.write_text(std::cout);
      std::cout << std::endl;
    }
}

template <int dim>
void 
HelmholtzProblem<dim>::output_results() const
{
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.attach_triangulation(triangulation);

  std::vector<std::string> solution_names;
  solution_names.emplace_back("Re(phi)");
  solution_names.emplace_back("Im(phi)");

  data_out.add_data_vector(solution, solution_names);

  // Write exact solution
  Vector<float> l2_error(triangulation.n_active_cells());
  VectorTools::integrate_difference(dof_handler,
                                    solution,
                                    MMSAnalyticalSolution<dim>(parameters),
                                    l2_error,
                                    QGauss<dim>(fe.degree + 1),
                                    VectorTools::L2_norm);
  data_out.add_data_vector(l2_error, "L2_error");
  
  ComputePressureIntensity<dim> pressure_intensity(parameters);
  data_out.add_data_vector(solution, pressure_intensity);

  ComputePressurePhase<dim> pressure_phase(parameters);
  data_out.add_data_vector(solution, pressure_phase);

  ComputeVelocityIntensity<dim> velocity_intensity(parameters);
  data_out.add_data_vector(solution, velocity_intensity);

  ComputeVelocityPhase<dim> velocity_phase(parameters);
  data_out.add_data_vector(solution, velocity_phase);

  // Add material_id as cell data for easy visualization
  data_out.attach_triangulation(triangulation);
  Vector<float> material_ids(triangulation.n_active_cells());
  for (const auto &cell : triangulation.active_cell_iterators())
    material_ids[cell->active_cell_index()] = cell->material_id();
  data_out.add_data_vector(material_ids, "material_id");

  data_out.build_patches();

  std::ofstream output(parameters.output_path + parameters.output_name + ".vtk");
  data_out.write_vtk(output);

  double period = 2 * M_PI / parameters.angular_frequency;
  double dt = period / parameters.n_timesteps;

  ComputeTimeHarmonicPressure<dim>  time_harmonic_pressure(parameters);
  ComputeTimeHarmonicVelocity<dim>  time_harmonic_velocity(parameters);
  ComputeFluidDisplacement<dim>     fluid_displacement(parameters);
  
  std::vector<std::pair<double, std::string>> times_and_names;
  for (unsigned int timestep = 0; timestep < parameters.n_timesteps + 1; ++timestep)
  {
    double time = timestep * dt;
    time_harmonic_pressure.set_time(time);
    time_harmonic_velocity.set_time(time);
    fluid_displacement.set_time(time);

    DataOut<dim> data_out_time_harmonic;
    data_out_time_harmonic.attach_dof_handler(dof_handler);

    data_out_time_harmonic.add_data_vector(solution, time_harmonic_pressure);
    data_out_time_harmonic.add_data_vector(solution, time_harmonic_velocity);
    data_out_time_harmonic.add_data_vector(solution, fluid_displacement);

    // Add material_id as cell data for easy visualization
    data_out_time_harmonic.attach_triangulation(triangulation);
    Vector<float> material_ids(triangulation.n_active_cells());
    for (const auto &cell : triangulation.active_cell_iterators())
      material_ids[cell->active_cell_index()] = cell->material_id();
    data_out_time_harmonic.add_data_vector(material_ids, "material_id");
  
    data_out_time_harmonic.build_patches();
    
    std::string filename = parameters.output_name + Utilities::int_to_string(timestep, 3) + ".vtu";
    std::ofstream output_time_harmonic(parameters.output_path + filename);
    data_out_time_harmonic.write_vtu(output_time_harmonic);

    times_and_names.emplace_back (time, filename);
  }

  std::ofstream pvd_output(parameters.output_path + parameters.output_name + ".pvd");
  DataOutBase::write_pvd_record(pvd_output, times_and_names);

  if (parameters.convergence_study)
  {
    std::ofstream table_file(parameters.output_path + "convergence_table_" + parameters.output_name + ".txt");
    convergence_table.write_text(table_file);
  }
}

template <int dim>
void 
HelmholtzProblem<dim>::run()
{
  std::cout << "Running Helmholtz problem in " << dim << "D" << std::endl;

  make_grid();

  for (unsigned int cycle = 0; cycle < parameters.n_cycles; ++cycle)
  { 
    if (cycle != 0)
      triangulation.refine_global(1);


    std::cout << "Starting cycle " << cycle << std::endl;
    std::cout << "    Setting up system" << std::endl;
    setup_system();

    std::cout << "    Assembling system" << std::endl;
    assemble_system();

    std::cout << "    Solving system" << std::endl;
    solve();

    if (parameters.convergence_study)
    {
      post_process_convergence(cycle);
    }

  }

  std::cout << "    Post processing" << std::endl;
  post_process_acoustic_potential();
  output_results();
}

int main(int argc, char *argv[])
{
  try
    { 
      using namespace dealii;

      Settings parameters;
      if (!parameters.try_parse((argc > 1) ? (argv[1]) : ""))
        return 0;

      switch (parameters.dimensions)
      {
        case 1:
          {
            HelmholtzProblem<1> helmholtz_problem(parameters);
            helmholtz_problem.run();
            break;
          }
        case 2:
          {
            HelmholtzProblem<2> helmholtz_problem(parameters);
            helmholtz_problem.run();
            break;
          }
        case 3:
          {
            HelmholtzProblem<3> helmholtz_problem(parameters);
            helmholtz_problem.run();
            break;
          }
        default:
          AssertThrow(false, ExcNotImplemented());
      }

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