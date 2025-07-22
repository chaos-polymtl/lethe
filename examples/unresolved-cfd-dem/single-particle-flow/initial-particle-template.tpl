# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# Listing of Parameters
#----------------------

set dimension = 3

#---------------------------------------------------
# Simulation Control
#---------------------------------------------------

subsection simulation control
  set time step        = 0.000001
  set time end         = 0.000001
  set output frequency = 1
  set output path      = ./{{case_folder}}/output_dem/
end

#---------------------------------------------------
# Timer
#---------------------------------------------------

subsection timer
  set type = none
end

#---------------------------------------------------
# Restart
#---------------------------------------------------

subsection restart
  set checkpoint = true
  set frequency  = 1
  set filename   = {{case_folder}}/dem
end

#---------------------------------------------------
# Physical Properties
#---------------------------------------------------

subsection lagrangian physical properties
  set g                        = 0.0, 0.0, 0.0
  set number of particle types = 1
  subsection particle type 0
    set size distribution type = uniform
    set diameter               = {{diameter}}
    set number of particles    = 1
    set density particles      = 1000
  end
end

#---------------------------------------------------
# Insertion Info
#---------------------------------------------------

subsection insertion info
  set insertion frequency = 1
  set insertion method    = list
  set list x              = {{x}}
  set list y              = {{y}}
  set list z              = {{z}}
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