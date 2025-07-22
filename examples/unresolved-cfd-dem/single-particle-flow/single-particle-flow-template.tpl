# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# Listing of Parameters
#----------------------

set dimension = 3

#---------------------------------------------------
# Simulation Control
#---------------------------------------------------

subsection simulation control
  set method           = steady
  set output frequency = 1
  set time end         = 2
  set time step        = 0.005
  set output path      = ./{{case_folder}}/output/
end

#---------------------------------------------------
# FEM
#---------------------------------------------------

subsection FEM
  set velocity order = {{velocity_order}}
  set pressure order = {{pressure_order}}
end

#---------------------------------------------------
# Restart
#---------------------------------------------------

subsection restart
  set checkpoint = false
  set frequency  = 10
  set restart    = false
  set filename   = case
end

#---------------------------------------------------
# Physical Properties
#---------------------------------------------------

subsection physical properties
  subsection fluid 0
    set kinematic viscosity = {{kinematic_viscosity}}
    set density             = 1
  end
end

#---------------------------------------------------
# Mesh
#---------------------------------------------------

subsection mesh
  set type           = dealii
  set grid type      = subdivided_hyper_cube
  set grid arguments = {{divisions}} : -1 : 1 : true
  set initial refinement = {{refinement}}
end

#---------------------------------------------------
# Void Fraction
#---------------------------------------------------

subsection void fraction
  set mode                         = qcm
  set qcm sphere equal cell volume = true
  set read dem                     = true
  set dem file name                = {{case_folder}}/dem
  set l2 smoothing length          = {{l2_smoothing_length}}
  set quadrature rule              = gauss-lobatto
  set n quadrature points          = 5
end

#---------------------------------------------------
# CFD-DEM
#---------------------------------------------------

subsection cfd-dem
  set grad div                      = true
  set void fraction time derivative = false
  set drag force                    = true
  set buoyancy force                = false
  set shear force                   = false
  set pressure force                = false
  set drag model                    = rong
  set coupling frequency            = 100
  set grad-div length scale         = {{grad_div}}
  set vans model                    = modelA
end

#---------------------------------------------------
# Boundary Conditions
#---------------------------------------------------

subsection boundary conditions
  set number = 6
  subsection bc 0
      set type = function
    subsection u
      set Function expression = 1
    end
    subsection v
      set Function expression = 0
    end
    subsection w
      set Function expression = 0
    end
  end
  subsection bc 1
    set type = outlet
  end
  subsection bc 2
    set type = slip
  end
  subsection bc 3
    set type = slip
  end
  subsection bc 4
    set type = slip
  end
  subsection bc 5
    set type = slip
  end
end

#---------------------------------------------------
# Model parameters
#---------------------------------------------------

subsection model parameters
  subsection contact detection
    set contact detection method = dynamic
    set neighborhood threshold   = 1.3
  end
end

#---------------------------------------------------
# Physical Properties
#---------------------------------------------------

subsection lagrangian physical properties
  set g                        = 0.0, -9.81, 0.0
  set number of particle types = 1
  subsection particle type 0
    set size distribution type = uniform
    set diameter               = {{diameter}}
    set number of particles    = 1
    set density particles      = 1000
  end
end

#---------------------------------------------------
# Non-Linear Solver Control
#---------------------------------------------------

subsection non-linear solver
  subsection fluid dynamics
    set verbosity = quiet
  end
end

#---------------------------------------------------
# Linear Solver Control
#---------------------------------------------------

subsection linear solver
  subsection fluid dynamics
    set method                  = gmres
    set max iters               = 500
    set relative residual       = 1e-3
    set minimum residual        = 1e-10
    set preconditioner          = ilu
    set ilu preconditioner fill = 1
    set verbosity               = quiet
    set max krylov vectors      = 500
  end
end