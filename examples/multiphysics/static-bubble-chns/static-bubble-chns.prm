# SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# Listing of Parameters
#----------------------

set dimension = 2

#---------------------------------------------------
# Simulation Control
#---------------------------------------------------

subsection simulation control
  set method           = bdf1
  set time end         = 1.5
  set time step        = 0.025
  set output name      = static-bubble-chns
  set output frequency = 3
  set output path      = ./
end

#---------------------------------------------------
# Multiphysics
#---------------------------------------------------

subsection multiphysics
  set cahn hilliard = true
end

#---------------------------------------------------
# Cahn Hilliard
#---------------------------------------------------

subsection cahn hilliard
  set potential smoothing coefficient = 0
  subsection epsilon
    set method    = automatic
    set verbosity = verbose
  end
end

#---------------------------------------------------
# Mesh
#---------------------------------------------------

subsection mesh
  set type               = dealii
  set grid type          = hyper_rectangle
  set grid arguments     = -2.5, -2.5 : 2.5, 2.5 : true
  set initial refinement = 9
end

#---------------------------------------------------
# Initial conditions
#---------------------------------------------------

subsection initial conditions
  subsection cahn hilliard
    set Function expression = -tanh(({{R}} - sqrt((x) * (x) + (y) * (y)))/(sqrt(2)*0.0138107)); 0
  end
end

#---------------------------------------------------
# Physical Properties
#---------------------------------------------------

subsection physical properties
  set number of fluids = 2
  subsection fluid 1
    set density             = 10
    set kinematic viscosity = 0.1
  end
  subsection fluid 0
    set density             = 10
    set kinematic viscosity = 0.1
  end
  set number of material interactions = 1
  subsection material interaction 0
    subsection fluid-fluid interaction
      set surface tension coefficient     = 1
      set cahn hilliard mobility constant = {{mobility}}
    end
  end
end

#---------------------------------------------------
# Boundary Conditions
#---------------------------------------------------

subsection boundary conditions
  set number = 4
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
end

subsection boundary conditions cahn hilliard
  set number = 4
  subsection bc 0
    set id   = 0
    set type = noflux
  end
  subsection bc 1
    set id   = 1
    set type = noflux
  end
  subsection bc 2
    set id   = 2
    set type = noflux
  end
  subsection bc 3
    set id   = 3
    set type = noflux
  end
end

#---------------------------------------------------
# Analytical solution
#---------------------------------------------------

subsection analytical solution
  set enable = true
  subsection uvwp
    set Function expression = 0; 0; if (x^2 + y^2 < {{R}}^2 , 2, 0)
  end
end

#---------------------------------------------------
# Post-processing
#---------------------------------------------------

subsection post-processing
  set verbosity        = verbose
  set output frequency = 1

  set calculate barycenter       = true
  set calculate phase statistics = true
end

#---------------------------------------------------
# Non-Linear Solver Control
#---------------------------------------------------

subsection non-linear solver
  subsection cahn hilliard
    set tolerance = 1e-6
  end
  subsection fluid dynamics
    set tolerance = 1e-4
  end
end

#---------------------------------------------------
# Linear Solver Control
#---------------------------------------------------

subsection linear solver
  subsection fluid dynamics
    set relative residual = 1e-3
    set minimum residual  = 1e-5
  end
  subsection cahn hilliard
    set relative residual = 1e-3
    set minimum residual  = 1e-7
  end
end
