# Listing of Parameters
#----------------------

set dimension = 3

#---------------------------------------------------
# Simulation Control
#---------------------------------------------------

subsection simulation control
  set time step         = 1e-7
  set time end          = 8
  set log frequency     = 100000
  set output frequency  = 100000
  set output path       = ./out_{{ER}}/
  set output name       = out
  set output boundaries = true
end

#---------------------------------------------------
# Model parameters
#---------------------------------------------------

subsection model parameters
  subsection contact detection
    set contact detection method                = dynamic
    set dynamic contact search size coefficient = 0.9
    set neighborhood threshold                  = 20
  end
  subsection load balancing
    set load balance method = once
    set step                = 150000
  end
  set particle particle contact force method = linear
  set particle wall contact force method     = linear
  set integration method                     = velocity_verlet
end

#---------------------------------------------------
# Physical Properties
#---------------------------------------------------

subsection lagrangian physical properties
  set g                        = 0.0, 0.0, -9.81
  set number of particle types = 1
  subsection particle type 0
    set size distribution type            = uniform
    set diameter                          = 0.2
    set number of particles               = 1
    set density particles                 = 2600
    set young modulus particles           = {{YP}}
    set poisson ratio particles           = 0.3
    set restitution coefficient particles = {{er}}
    set friction coefficient particles    = 0.0
  end
  set young modulus wall           = 1000000000000
  set poisson ratio wall           = 0.30
  set restitution coefficient wall = {{er}}
  set friction coefficient wall    = 0.
end

#---------------------------------------------------
# Insertion Info
#---------------------------------------------------

subsection insertion info
  set insertion method    = list
  set insertion frequency = 10000
  set list x              = 1.
  set list y              = 1.
  set list z              = 0.5
  set list velocity x     = 0.0
  set list velocity y     = 0.0
  set list velocity z     = 0.0
  set list omega x        = 0.0
  set list omega y        = 0.0
  set list omega z        = 0.0
  set list diameters      = 0.2
end

#---------------------------------------------------
# Mesh
#---------------------------------------------------

subsection mesh
  set type                                = dealii
  set grid type                           = hyper_cube
  set grid arguments                      = 0 : 2 : false
  set initial refinement                  = 4
  set expand particle-wall contact search = false
end

#---------------------------------------------------
# Boundary Condition
#---------------------------------------------------

subsection DEM boundary conditions
  set number of boundary conditions = 1
  subsection boundary condition 0
    set boundary id         = 0
    set type                = fixed_wall
  end
end
