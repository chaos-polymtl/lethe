# Listing of Parameters
# ---------------------

set dimension = 2

#---------------------------------------------------
# Simulation and IO Control
#---------------------------------------------------

subsection simulation control
  set method           = bdf2
  set time end         = ENDTIME
  set time step        = TIMESTEP
  set output name      = capillary-wave-TSM-TIMESTEPMULTIPLIER
  set output frequency = OUTPUTFREQUENCY
  set output path      = ./output-TSM-TIMESTEPMULTIPLIER/
end

#---------------------------------------------------
# Multiphysics
#---------------------------------------------------

subsection multiphysics
  set VOF = true
end

#---------------------------------------------------
# VOF
#---------------------------------------------------

subsection VOF
  subsection surface tension force
    set enable                                   = true
    set phase fraction gradient diffusion factor = 4
    set curvature diffusion factor               = 1
    set output auxiliary fields                  = true
  end
  subsection phase filtration
    set type      = tanh
    set beta      = 10
  end
end

#---------------------------------------------------
# Initial condition
#---------------------------------------------------

subsection initial conditions
  set type = nodal
  subsection uvwp
    set Function expression = 0; 0; 0
  end
  subsection VOF
    set Function expression = if (y<=1e-6*cos(2*3.14159/1e-4*x), min(0.5-(y-1e-6*cos(2*3.14159/1e-4*x))/1e-6,1), max(0.5-(y-1e-6*cos(2*3.14159/1e-4*x))/1e-6,0))
    subsection projection step
      set enable           = true
      set diffusion factor = 1
    end
  end
end

#---------------------------------------------------
# Source term
#---------------------------------------------------

subsection source term
  subsection xyz
    set Function expression = 0; 0; 0
  end
end

#---------------------------------------------------
# Physical Properties
#---------------------------------------------------

subsection physical properties
  set number of fluids = 2
  subsection fluid 1
    set density             = 1
    set kinematic viscosity = 5e-6
  end
  subsection fluid 0
    set density             = 1
    set kinematic viscosity = 5e-6
  end
  set number of material interactions = 1
  subsection material interaction 0
    set type = fluid-fluid
    subsection fluid-fluid interaction
      set first fluid id              = 0
      set second fluid id             = 1
      set surface tension model       = constant
      set surface tension coefficient = 0.01
    end
  end
end

#---------------------------------------------------
# Mesh
#---------------------------------------------------

subsection mesh
  set type               = dealii
  set grid type          = subdivided_hyper_rectangle
  set grid arguments     = 4, 12 : -5e-5, -1.5e-4 : 5e-5, 1.5e-4 : true
  set initial refinement = 4
end

#---------------------------------------------------
# Mesh Adaptation
#---------------------------------------------------

subsection mesh adaptation
  set type                     = kelly
  set variable                 = phase
  set fraction type            = fraction
  set max refinement level     = 5
  set min refinement level     = 3
  set frequency                = 1
  set fraction refinement      = 0.95
  set fraction coarsening      = 0.05
  set initial refinement steps = 4
end


# --------------------------------------------------
# Boundary Conditions
#---------------------------------------------------

subsection boundary conditions
  set number = 3
  subsection bc 0
    set id                 = 0
    set type               = periodic
    set periodic_id        = 1
    set periodic_direction = 0
  end
  subsection bc 1
    set id   = 2
    set type = noslip
  end
  subsection bc 2
    set id   = 3
    set type = noslip
  end
end

#---------------------------------------------------
# FEM
#---------------------------------------------------

subsection FEM
  set velocity order = 1
  set pressure order = 1
  set VOF order      = 1
end

# --------------------------------------------------
# Non-Linear Solver Control
#---------------------------------------------------

subsection non-linear solver
  subsection fluid dynamics
    set tolerance      = 1e-9
    set max iterations = 20
    set verbosity      = quiet
  end
  subsection VOF
    set tolerance      = 1e-10
    set max iterations = 2
    set verbosity      = quiet
  end
end

# --------------------------------------------------
# Linear Solver Control
#---------------------------------------------------

subsection linear solver
  subsection fluid dynamics
    set method            = gmres
    set relative residual = 1e-3
    set minimum residual  = 1e-10
    set preconditioner    = ilu
  end
  subsection VOF
    set method            = gmres
    set relative residual = 1e-4
    set minimum residual  = 1e-12
    set preconditioner    = ilu
  end
end
