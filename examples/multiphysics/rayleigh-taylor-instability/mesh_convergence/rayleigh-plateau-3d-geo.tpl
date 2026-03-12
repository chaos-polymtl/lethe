# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# Listing of Parameters
# ---------------------

set dimension        = 3
set print parameters = all

#---------------------------------------------------
# Simulation and IO Control
#---------------------------------------------------

subsection simulation control
  set method           = bdf2
  set time end         = 0.08
  set time step        = MAX_TIME_STEP
  set output name      = CASE_NAME
  set output frequency = 200
  set output path      = ./output/
  set adapt            = true
  set max cfl          = 0.25
  set max time step    = MAX_TIME_STEP
  set group files      = 96
  set subdivision      = 1
  #set output time interval = 0.05, 0.08
end

#---------------------------------------------------
# Multiphysics
#---------------------------------------------------

subsection multiphysics
  set CLS = true
end

#---------------------------------------------------
# CLS
#---------------------------------------------------

subsection CLS
  subsection surface tension force
    set enable                  = true
    set output auxiliary fields = true
  end
  subsection interface regularization method
    set type      = REGULARIZATION_TYPE
    set frequency = REGULARIZATION_FREQUENCY
    set verbosity = extra verbose
    subsection projection-based interface sharpening
      set interface sharpness = 1.5
    end
    subsection pde-based interface reinitialization
      set steady-state criterion = 1e-4
      set max steps number       = 10000
      set diffusivity multiplier = DIFFUSIVITY_MULT
      set diffusivity power      = 1.0
      set reinitialization CFL   = 0.25
    end
    subsection geometric interface reinitialization
      set max reinitialization distance = REGULARIZATION_DISTANCE
      set tanh thickness                = TANH_THICKNESS
      set transformation type           = tanh
    end
  end
  subsection phase filtration
    set type = tanh
    set beta = 20
  end
end

#---------------------------------------------------
# Initial condition
#---------------------------------------------------

subsection initial conditions
  set type = nodal
  subsection uvwp
    set Function constants  = U=1.569
    set Function expression = if((y^2 + z^2) <= 1.3110e-6, U, 0); 0; 0; 0
  end
  subsection CLS
    set Function constants  = r=1.145e-3
    set Function expression = 0.5 - 0.5 * tanh(sign(sqrt(y*y + z*z)-r)*min(abs(sqrt(y*y + z*z)-r), REGULARIZATION_DISTANCE) / TANH_THICKNESS)
    set smoothing type      = none
  end
end

#---------------------------------------------------
# Physical Properties
#---------------------------------------------------

subsection physical properties
  set number of fluids = 2
  subsection fluid 0
    set density             = 1.196
    set kinematic viscosity = 2.54e-4
  end
  subsection fluid 1
    set density             = 1196
    set kinematic viscosity = 2.54e-5
  end
  set number of material interactions = 1
  subsection material interaction 0
    set type = fluid-fluid
    subsection fluid-fluid interaction
      set surface tension coefficient = 0.0674
    end
  end
end

#---------------------------------------------------
# Mesh
#---------------------------------------------------

subsection mesh
  set type               = dealii
  set grid type          = subdivided_hyper_rectangle
  set grid arguments     = 8, 2, 2: -0.0458, -0.01145, -0.01145 : 0.0458, 0.01145, 0.01145 : true
  set initial refinement = INITIAL_REFINEMENT
end

#---------------------------------------------------
# Mesh Adaptation
#---------------------------------------------------

subsection mesh adaptation
  set type                     = adaptive
  set error estimator          = kelly
  set variable                 = phase
  set fraction type            = fraction
  set max refinement level     = MAX_REFINEMENT
  set min refinement level     = MIN_REFINEMENT
  set fraction refinement      = 0.99
  set fraction coarsening      = 0.001
  set initial refinement steps = 4
  set frequency                = 10
end

# --------------------------------------------------
# Boundary Conditions
#---------------------------------------------------

subsection boundary conditions
  set number         = 6
  set time dependent = true
  subsection bc 0
    set id   = 0
    set type = function
    subsection u
      set Function constants  = U=1.569, delta=0.3, kappa=0.7, r=1.145e-3
      set Function expression = if ((y^2 + z^2) <= 1.3110e-6, U*(1 + delta*sin(kappa*U*t/r)), 0)
    end
  end
  subsection bc 1
    set id   = 1
    set type = outlet
    set beta = 0
  end
  subsection bc 2
    set id   = 2
    set type = outlet
    set beta = 0
  end
  subsection bc 3
    set id   = 3
    set type = outlet
    set beta = 0
  end
  subsection bc 4
    set id   = 4
    set type = outlet
    set beta = 0
  end
  subsection bc 5
    set id   = 5
    set type = outlet
    set beta = 0
  end
end

# --------------------------------------------------
# Boundary Conditions CLS
#---------------------------------------------------

subsection boundary conditions CLS
  set number         = 6
  set time dependent = true
  subsection bc 0
    set id   = 0
    set type = dirichlet
    subsection dirichlet
      set Function constants  = r=1.145e-3
      set Function expression = 0.5 - 0.5 * tanh(sign(sqrt(y*y + z*z)-r)*min(abs(sqrt(y*y + z*z)-r), REGULARIZATION_DISTANCE) / TANH_THICKNESS)
    end
  end
  subsection bc 1
    set id = 1
  end
  subsection bc 2
    set id = 2
  end
  subsection bc 3
    set id = 3
  end
  subsection bc 4
    set id = 4
  end
  subsection bc 5
    set id = 5
  end
end

# --------------------------------------------------
# Non-Linear Solver Control
#---------------------------------------------------

subsection non-linear solver
  subsection fluid dynamics
    set tolerance      = 1e-7
    set max iterations = 20
  end
  subsection CLS
    set tolerance      = 1e-10
    set max iterations = 20
  end
end

# --------------------------------------------------
# Linear Solver Control
#---------------------------------------------------

subsection linear solver
  subsection fluid dynamics
    set relative residual       = 5e-4
    set minimum residual        = 1e-8
    set ilu preconditioner fill = 0
  end
  subsection CLS
    set relative residual = 1e-8
    set minimum residual  = 1e-11
  end
end

# --------------------------------------------------
# Timer
#---------------------------------------------------

subsection timer
  set type = iteration
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

subsection stabilization
  set cls dcdd stabilization = false
end
