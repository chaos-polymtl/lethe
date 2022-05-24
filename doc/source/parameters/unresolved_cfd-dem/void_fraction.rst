***********************************************
Void Fraction
***********************************************
In this subsection, all parameters required for the calculation of the void fraction are introduced.

.. code-block:: text

   subsection void fraction
      set mode = pcm
      set read dem = true
      set dem file name = dem
      set l2 smoothing factor = 0
      set l2 lower bound = 0
      set l2 upper bound = 1
      set bound void fraction = false
      set particle refinement factor = 0
   end

* The ``mode`` parameter allows the user to choose the method for void fraction calculation. Currently, there are two methods implemented. The first one is to calculate the void fraction from `function`. In this case, an additional subsection is required to insert the function such as:

.. code-block:: text

   subsection function
   set Function expression = 0.5 + 0.25*sin(pi*x)*sin(pi*y)
   end
     
If the ``mode`` chosen is ``pcm``, then the void fraction is calculated using the Particle Centered Method. If it is set to ``qcm``, then the void fraction is calculated using the Quadrature Centered Method. If it is set to ``spm``, then the void fraction is calculated using the satellite point method (divided approach). In these methods, the remaining parameters are required:

* The ``read dem`` allows us to read an already existing dem simulation result which can be obtained from checkpointing the Lethe-DEM simulation. This is important as the gls_vans solver requires reading an initial dem triangulation and particle information to simulate flows in the presence of particles. 
* The ``dem_file_name`` parameter specifies the prefix of the dem files that must be read.
* The ``l2 smoothing factor`` is a smoothing length used for smoothing the L2 projection of the void fraction to avoid sharp discontinuities which can lead to instabilities in the simulation.
* The ``l2 lower bound`` and ``l2 upper bound`` are the minimum and maximum values around which the void fraction is bounded. This is important especially for upper bounds as the void fraction can sometimes slightly exceed a value of 1 when projected.
* The ``bound void fraction`` parameter determines whether or not to bound the void fraction between the lower and upper bounds specified in the previous two parameters. As the void fraction is calculated and then projected using L2 projection, it can sometimes exceeds the maximum value of 1. In order to prevent this, we use an active set method to bound the void fraction.
* The ``particle refinement factor`` is only required for the ``spm``. It allows to determine the number of pseudo-particles that we want to divide our particle into. By default, it is set to 0 refinements, and results in no refinement of the original meshed particle (division into 7 particles in 3D). Every additional refinement results in a :math:`2^{dim}` times more particles. The figure below shows how the number of pseudo-particles change with every refinement. Every cell in the particle mesh represents a pseudo-particle in the satellite point method

.. image:: images/refinement.png
