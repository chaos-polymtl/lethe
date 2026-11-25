# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# Template of the parameter file of the Reynolds number sweep of the backward-facing
# step. The cases are generated with generate_cases_locally.py, which replaces the
# kinematic viscosity by its value for every Reynolds number (Re = 2/nu).

# Listing of Parameters
#----------------------

set dimension = 2

#---------------------------------------------------
# Simulation Control
#---------------------------------------------------

subsection simulation control
  set method            = steady
  set number mesh adapt = 10
  set output name       = backward_facing_step_output
  set output frequency  = 1
end

#---------------------------------------------------
# Physical Properties
#---------------------------------------------------

subsection physical properties
  set number of fluids = 1
  subsection fluid 0
    set kinematic viscosity = {{ kinematic_viscosity }} # Re = {{ reynolds }}
  end
end

#---------------------------------------------------
# Initial Conditions
#---------------------------------------------------

# The steady solution is obtained by a continuation on the kinematic viscosity,
# starting from Re = 100, so that every Reynolds number is solved in the same way

subsection initial conditions
  set type = ramp
  subsection ramp
    subsection kinematic viscosity
      set initial kinematic viscosity = 0.02
      set iterations                  = 8
      set alpha                       = 0.5
    end
  end
end

#---------------------------------------------------
# Mesh
#---------------------------------------------------

subsection mesh
  set type               = gmsh
  set file name          = ../../backward-facing-step.msh
  set initial refinement = 1
end

#---------------------------------------------------
# Mesh Adaptation Control
#---------------------------------------------------

subsection mesh adaptation
  set variable             = velocity
  set type                 = adaptive
  set error estimator      = kelly
  set min refinement level = 1
  set fraction refinement  = 0.1
  set fraction coarsening  = 0
end

#---------------------------------------------------
# Post-processing
#---------------------------------------------------

subsection post-processing
  set output velocity gradient = true
end

#---------------------------------------------------
# FEM
#---------------------------------------------------

subsection FEM
  set pressure degree = 1
  set velocity degree = 1
end

#---------------------------------------------------
# Boundary Conditions
#---------------------------------------------------

subsection boundary conditions
  set number = 3
  subsection bc 0
    set id   = 0
    set type = noslip
  end
  subsection bc 1
    set id   = 1
    set type = function
    subsection u
      set Function expression = 6*(y-1)*(2-y)
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
    set type = outlet
  end
end

#---------------------------------------------------
# Non-Linear Solver Control
#---------------------------------------------------

subsection non-linear solver
  subsection fluid dynamics
    set verbosity      = verbose
    set tolerance      = 1e-6
    set max iterations = 20
  end
end

#---------------------------------------------------
# Linear Solver Control
#---------------------------------------------------

subsection linear solver
  subsection fluid dynamics
    set verbosity                   = verbose
    set method                      = gmres
    set max iters                   = 500
    set max krylov vectors          = 500
    set relative residual           = 1e-4
    set minimum residual            = 1e-9
    set preconditioner              = amg
    set amg preconditioner ilu fill = 1
  end
end
