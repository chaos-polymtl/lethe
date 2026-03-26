# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#----------------------
# Listing of Parameters
#----------------------

set dimension = 3

# ---------------------
# Simulation and IO Control
#---------------------------------------------------

subsection simulation control
  set method           = bdf2
  set time end         = 3.0
  set time step        = 0.0013
  set max time step    = 0.0013
  set output name      = CASE_NAME
  set output frequency = 100
  set group files      = 192
  set output path      = ./output/
end

#---------------------------------------------------
# Multiphysics
#---------------------------------------------------

subsection multiphysics
  set cls           = true
  set heat transfer = true
end

#---------------------------------------------------
# CLS
#---------------------------------------------------

subsection CLS
  subsection phase filtration
    set type = tanh
    set beta = 20
  end
  subsection interface reinitialization method
    set type      = REINITIALIZATION_TYPE
    set frequency = REINITIALIZATION_FREQUENCY
    set verbosity = verbose
    subsection projection-based interface sharpening
      set interface sharpness = 1.5
      set type                = constant
      set max iterations      = 50
      set tolerance           = 1e-7
    end
    subsection PDE-based interface reinitialization
      set output reinitialization steps = false
      set steady-state criterion        = 1e-4
      set max steps number              = 10000
      set diffusivity multiplier        = EPSILON
      set diffusivity power             = 1.0
      set reinitialization CFL          = 0.25
    end
    subsection geometric interface reinitialization
      set max reinitialization distance = REINITIALIZATION_DISTANCE
      set transformation type           = tanh
      set tanh thickness                = TANH_THICKNESS
    end
  end
  subsection surface tension force
    set enable                                    = true
    set phase indicator gradient diffusion factor = 4
    set curvature diffusion factor                = 1
    set output auxiliary fields                   = true
    set enable marangoni effect                   = true
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
    set Function constants  = r=0.25
    set Function expression = 0.5 - 0.5 * tanh((sqrt(x*x + y*y + z*z)-r) / TANH_THICKNESS)
  end
  subsection temperature
    set Function expression = x+1.5
  end
end

#---------------------------------------------------
# Physical Properties
#---------------------------------------------------

subsection physical properties
  set number of fluids = 2
  subsection fluid 1
    set density              = 1
    set kinematic viscosity  = 1
    set thermal conductivity = 1000000
  end
  subsection fluid 0
    set density              = 1
    set kinematic viscosity  = 1
    set thermal conductivity = 1000000
  end
  set number of material interactions = 1
  subsection material interaction 0
    set type = fluid-fluid
    subsection fluid-fluid interaction
      set first fluid id                              = 0
      set second fluid id                             = 1
      set surface tension model                       = linear
      set surface tension coefficient                 = 3.0
      set reference state temperature                 = 0.0
      set temperature-driven surface tension gradient = -1.0
    end
  end
end

#---------------------------------------------------
# Mesh
#---------------------------------------------------

subsection mesh
  set type               = dealii
  set grid type          = subdivided_hyper_rectangle
  set grid arguments     = 3, 3, 3: -1.5, -1.5, -1.5: 1.5, 1.5, 1.5: true
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
  set frequency                = 20
  set fraction refinement      = 0.99
  set fraction coarsening      = 0.001
  set initial refinement steps = 4
end

# --------------------------------------------------
# Boundary Conditions
#---------------------------------------------------

subsection boundary conditions
  set number = 6
  subsection bc 0
    set id   = 0
    set type = function
    subsection u
      set Function expression = 0
    end
    subsection v
      set Function expression = 0
    end
    subsection w
      set Function expression = 0
    end
  end
  subsection bc 1
    set id   = 1
    set type = function
    subsection u
      set Function expression = 0
    end
    subsection v
      set Function expression = 0
    end
    subsection w
      set Function expression = 0
    end
  end
  subsection bc 2
    set id   = 2
    set type = function
    subsection u
      set Function expression = 0
    end
    subsection v
      set Function expression = 0
    end
    subsection w
      set Function expression = 0
    end
  end
  subsection bc 3
    set id   = 3
    set type = function
    subsection u
      set Function expression = 0
    end
    subsection v
      set Function expression = 0
    end
    subsection w
      set Function expression = 0
    end
  end
  subsection bc 4
    set id   = 4
    set type = function
    subsection u
      set Function expression = 0
    end
    subsection v
      set Function expression = 0
    end
    subsection w
      set Function expression = 0
    end
  end
  subsection bc 5
    set id   = 5
    set type = function
    subsection u
      set Function expression = 0
    end
    subsection v
      set Function expression = 0
    end
    subsection w
      set Function expression = 0
    end
  end
end

subsection boundary conditions heat transfer
  set number = 6
  subsection bc 0
    set id   = 0
    set type = temperature
    subsection value
      set Function expression = 0
    end
  end
  subsection bc 1
    set id   = 1
    set type = temperature
    subsection value
      set Function expression = 3
    end
  end
  subsection bc 2
    set id   = 2
    set type = noflux
  end
  subsection bc 3
    set id   = 3
    set type = noflux
  end
  subsection bc 4
    set id   = 4
    set type = noflux
  end
  subsection bc 5
    set id   = 5
    set type = noflux
  end
end

#---------------------------------------------------
# Post-processing
#---------------------------------------------------

subsection post-processing
  set calculate barycenter     = true
  set calculate kinetic energy = true
end

#---------------------------------------------------
# FEM
#---------------------------------------------------

subsection FEM
  set velocity order    = 1
  set pressure order    = 1
  set cls order         = 1
  set temperature order = 1
end

# --------------------------------------------------
# Non-Linear Solver Control
#---------------------------------------------------

subsection non-linear solver
  subsection CLS
    set tolerance      = 1e-8
    set max iterations = 20
    set verbosity      = verbose
  end
  subsection heat transfer
    set tolerance      = 1e-6
    set max iterations = 20
    set verbosity      = verbose
  end
  subsection fluid dynamics
    set tolerance      = 1e-6
    set max iterations = 20
    set verbosity      = verbose
  end
  subsection CLS PDE-based interface reinitialization
    set tolerance      = 1e-8
    set max iterations = 20
    set verbosity      = verbose
  end
end

# --------------------------------------------------
# Boundary Conditions CLS
#---------------------------------------------------

subsection boundary conditions CLS
  set number = 6
end

# --------------------------------------------------
# Linear Solver Control
#---------------------------------------------------

subsection linear solver
  subsection fluid dynamics
    set verbosity                             = verbose
    set method                                = gmres
    set max iters                             = 8000
    set relative residual                     = 1e-8
    set minimum residual                      = 1e-8
    set ilu preconditioner fill               = 0
    set ilu preconditioner absolute tolerance = 1e-12
    set ilu preconditioner relative tolerance = 1.00
    set max krylov vectors                    = 200
  end
  subsection CLS
    set verbosity                             = verbose
    set method                                = gmres
    set max iters                             = 8000
    set relative residual                     = 1e-8
    set minimum residual                      = 1e-11
    set ilu preconditioner fill               = 0
    set ilu preconditioner absolute tolerance = 1e-12
    set ilu preconditioner relative tolerance = 1.00
    set max krylov vectors                    = 200
  end
  subsection heat transfer
    set verbosity                             = verbose
    set method                                = gmres
    set max iters                             = 8000
    set relative residual                     = 1e-8
    set minimum residual                      = 1e-8
    set ilu preconditioner fill               = 0
    set ilu preconditioner absolute tolerance = 1e-12
    set ilu preconditioner relative tolerance = 1.00
    set max krylov vectors                    = 200
  end
  subsection CLS PDE-based interface reinitialization
    set verbosity                             = verbose
    set method                                = gmres
    set max iters                             = 8000
    set relative residual                     = 1e-8
    set minimum residual                      = 1e-11
    set ilu preconditioner fill               = 0
    set ilu preconditioner absolute tolerance = 1e-12
    set ilu preconditioner relative tolerance = 1.00
    set max krylov vectors                    = 200
  end
end

# --------------------------------------------------
# Restart
#---------------------------------------------------

subsection restart
  set checkpoint = true
  set frequency  = 100
  set filename   = restart
  set restart    = false
end

# --------------------------------------------------
# Stabilization
#---------------------------------------------------

subsection stabilization
  set cls dcdd stabilization = false
end

# --------------------------------------------------
# Timer
#---------------------------------------------------

subsection timer
  set type = iteration
end
