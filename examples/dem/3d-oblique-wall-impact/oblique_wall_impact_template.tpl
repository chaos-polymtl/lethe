# Listing of Parameters
#----------------------

set dimension = 3

#---------------------------------------------------
# Simulation Control
#---------------------------------------------------

subsection simulation control
  set time step         = 2e-10
  set time end          = 0.01
  set log frequency     = 50000000
  set output frequency  = 50000000
  set output path       = ./{{FN}}/
  set output name       = out
  set output boundaries = true
end

#---------------------------------------------------
# Model parameters
#---------------------------------------------------

subsection model parameters
  set particle particle contact force method = hertz_mindlin_limit_overlap
  set particle wall contact force method     = nonlinear
  set integration method                     = velocity_verlet
  set rolling resistance torque method       = no_resistance

end

#---------------------------------------------------
# Physical Properties
#---------------------------------------------------

subsection lagrangian physical properties
  set g                        = 0.0, 0.0, 0.0
  set number of particle types = 1
  subsection particle type 0
    set size distribution type                 = uniform
    set diameter                               = 0.005
    set number of particles                    = 1
    set density particles                      = 4000
    set young modulus particles                = 380e9
    set poisson ratio particles                = 0.23
    set restitution coefficient particles      = 1.0
    set friction coefficient particles         = 0.092

  end
  set young modulus wall           = 70e9
  set poisson ratio wall           = 0.25
  set restitution coefficient wall = 1.0
  set friction coefficient wall    = 0.092
end

#---------------------------------------------------
# Insertion Info
#---------------------------------------------------

subsection insertion info
  set insertion method    = list
  set insertion frequency = 10000
  set list x              = 1.
  set list y              = 1.
  set list z              = 0.005
  set list velocity x     = 0.0
  set list velocity y     = {{vy}}
  set list velocity z     = {{vz}}
  set list omega x        = 0.0
  set list omega y        = 0.0
  set list omega z        = 0.0
  set list diameters      = 0.005
end

#---------------------------------------------------
# Mesh
#---------------------------------------------------

subsection mesh
  set type                                = dealii
  set grid type                           = hyper_cube
  set grid arguments                      = 0 : 2 : false
  set initial refinement                  = 0
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
