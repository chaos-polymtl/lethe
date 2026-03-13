# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# Listing of Parameters
#----------------------

set dimension = 3

#---------------------------------------------------
# Simulation Control
#---------------------------------------------------

subsection simulation control
  set method           = bdf2
  set time end         = 3.0
  set time step        = 0.0013
  set output name      = CASE_NAME
  set output frequency = 20
  set group files      = 192
  set output path      = ./output/
end

#---------------------------------------------------
# Multiphysics
#---------------------------------------------------

subsection multiphysics
  set cls = true
end

#---------------------------------------------------
# CLS
#---------------------------------------------------

subsection CLS
  subsection interface regularization method
    set type      = REGULARIZATION_TYPE
    set frequency = REGULARIZATION_FREQUENCY
    set verbosity = extra verbose
    subsection projection-based interface sharpening
      set threshold           = 0.5
      set interface sharpness = 1.5
    end
    subsection geometric interface reinitialization
      # TODO after reinitialization initialization fix
      set max reinitialization distance = REGULARIZATION_DISTANCE
      set tanh thickness                = TANH_THICKNESS
      set transformation type           = tanh
    end
    subsection PDE-based interface reinitialization
      set output reinitialization steps = false
      set steady-state criterion        = 1e-4
      set max steps number              = 10000
      set diffusivity multiplier        = EPSILON
      set diffusivity power             = 1.0
      set reinitialization CFL          = 0.25
    end
  end
  subsection phase filtration
    set type = tanh
    set beta = 20
  end
  subsection surface tension force
    set enable                                    = true
    set phase indicator gradient diffusion factor = 4
    set curvature diffusion factor                = 1
    set output auxiliary fields                   = true
  end
end

#---------------------------------------------------
# Initial condition
#---------------------------------------------------

subsection initial conditions
  set type = nodal
  subsection uvwp
    set Function expression = 0; 0; 0; 0
  end
  subsection CLS
    set Function constants  = center=0.5
    set Function expression = 0.5 - 0.5*tanh((sqrt((x-center)*(x-center)+(y-center)*(y-center)+(z-center)*(z-center))-0.25)/(TANH_THICKNESS))
  end
end

#---------------------------------------------------
# Source term
#---------------------------------------------------

subsection source term
  subsection fluid dynamics
    set Function expression = 0; -0.98; 0; 0
  end
end

#---------------------------------------------------
# Physical Properties
#---------------------------------------------------

subsection physical properties
  set number of fluids = 2
  subsection fluid 0
    set density             = 1000
    set kinematic viscosity = 0.01
  end
  subsection fluid 1
    set density             = 100
    set kinematic viscosity = 0.01
  end
  set number of material interactions = 1
  subsection material interaction 0
    set type = fluid-fluid
    subsection fluid-fluid interaction
      set first fluid id              = 0
      set second fluid id             = 1
      set surface tension model       = constant
      set surface tension coefficient = 24.5
    end
  end
end

#---------------------------------------------------
# Post-processing
#---------------------------------------------------

subsection post-processing
  set calculate mass conservation = true
  set calculate barycenter        = true
  set barycenter name             = vof_barycenter_information
end

#---------------------------------------------------
# Mesh
#---------------------------------------------------

subsection mesh
  set type               = dealii
  set grid type          = subdivided_hyper_rectangle
  set grid arguments     = 1, 2, 1 : 0, 0, 0 : 1, 2, 1 : true
  set initial refinement = 4
end

#---------------------------------------------------
# Mesh Adaptation
#---------------------------------------------------

subsection mesh adaptation
  set type                     = adaptive
  set error estimator          = kelly
  set variable                 = phase
  set fraction type            = fraction
  set max refinement level     = 7
  set min refinement level     = 4
  set frequency                = 5
  set fraction refinement      = 0.99
  set fraction coarsening      = 0.001
  set initial refinement steps = 4
end

#---------------------------------------------------
# Timer
#---------------------------------------------------

subsection timer
  set type = end
end

#---------------------------------------------------
# Boundary Conditions
#---------------------------------------------------

subsection boundary conditions
  set number = 6
  subsection bc 0
    set id   = 0
    set type = noslip
  end
  subsection bc 1
    set id   = 1
    set type = noslip
  end
  subsection bc 2
    set id   = 2
    set type = noslip
  end
  subsection bc 3
    set id   = 3
    set type = noslip
  end
  subsection bc 4
    set id   = 4
    set type = noslip
  end
  subsection bc 5
    set id   = 5
    set type = noslip
  end
end

subsection boundary conditions CLS
  set number = 6
end

#---------------------------------------------------
# FEM
#---------------------------------------------------

subsection FEM
  set velocity order = 1
  set pressure order = 1
end

#---------------------------------------------------
# Non-Linear Solver Control
#---------------------------------------------------

subsection non-linear solver
  subsection fluid dynamics
    set tolerance      = 1e-5
    set max iterations = 20
    set verbosity      = verbose
  end
  subsection CLS
    set tolerance      = 1e-11
    set max iterations = 20
    set verbosity      = verbose
  end
  subsection CLS PDE-based interface reinitialization
    set tolerance      = 1e-11
    set max iterations = 20
    set verbosity      = verbose
  end
end

#---------------------------------------------------
# Linear Solver Control
#---------------------------------------------------

subsection linear solver
  subsection fluid dynamics
    set verbosity                             = verbose
    set method                                = gmres
    set max iters                             = 8000
    set relative residual                     = 1e-4
    set minimum residual                      = 1e-7
    set preconditioner                        = ilu
    set ilu preconditioner fill               = 0
    set ilu preconditioner absolute tolerance = 1e-12
    set ilu preconditioner relative tolerance = 1.00
    set max krylov vectors                    = 200
  end
  subsection CLS
    set verbosity                             = verbose
    set method                                = gmres
    set max iters                             = 8000
    set relative residual                     = 1e-9
    set minimum residual                      = 1e-12
    set preconditioner                        = ilu
    set ilu preconditioner fill               = 0
    set ilu preconditioner absolute tolerance = 1e-12
    set ilu preconditioner relative tolerance = 1.00
    set max krylov vectors                    = 200
  end
  subsection CLS PDE-based interface reinitialization
    set verbosity                             = verbose
    set method                                = gmres
    set max iters                             = 8000
    set relative residual                     = 1e-9
    set minimum residual                      = 1e-12
    set preconditioner                        = ilu
    set ilu preconditioner fill               = 1
    set ilu preconditioner absolute tolerance = 1e-12
    set ilu preconditioner relative tolerance = 1.00
    set max krylov vectors                    = 200
  end
end

# --------------------------------------------------
# Stabilization
#---------------------------------------------------

subsection stabilization
  set cls dcdd stabilization = false
end

# --------------------------------------------------
# Restart
#---------------------------------------------------

subsection restart
  set checkpoint = true
  set frequency  = 50
  set filename   = restart
  set restart    = false
end

#---------------------------------------------------
# Timer
#---------------------------------------------------

subsection timer
  set type = iteration
end
