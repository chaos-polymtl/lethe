=====================
Particle Ray Tracing
=====================

In this subsection, the parameters used by the ``lethe-particles-ray-tracing`` application are explained. This application simulates the propagation of photons (or rays) through the simulation domain to reconstruct the surface formed by particles, using the same principles as profilometry. The parameters defined here control the initial conditions for photon insertion, including their starting positions, directions, and any random offsets applied to their trajectories.

.. code-block:: text

  subsection particle ray tracing
    set starting photon insertion position                  = 0., 0., 0.
    set insertion unit tensors                              = 1.,0.,0. : 0., 1., 0. : 0., 0., 1.
    set number of inserted photons per direction            = 1 : 1 : 1
    set distance between photons on insertion per direction = 1 : 1 : 1
    set reference displacement vector                       = 0., 0., 1.
    set photon insertion maximum offset                     = 0.
    set photon insertion prn seed                           = 0
    set photon maximum angular offset                       = 0.
    set photon angular offset prn seed                      = 1
  end

-  ``starting photon insertion position`` is the location of the first photon being inserted, given as a 3D coordinate (x,y,z).

-  ``insertion unit tensors`` are three vectors defining the insertion grid directions relative to the ``starting photon position``. By default, these correspond to the Cartesian axes (x-y-z). The vectors do not need to be orthogonal or normalized, but they must not be non-zero vectors. The insertion grid is formed by creating insertion points relative to the ``starting photon insertion position`` using these three vectors.

-  ``number of inserted photons per direction`` are the number of photons to insert along each of the three ``insertion unit tensors`` . Example: ``5 : 10 : 1`` inserts 5 × 10 × 1 = 50 photons.

-  ``distance between photons on insertion per direction`` is the spacing between photons along each insertion direction when the ``photon insertion maximum offset`` is zero. All values must be strictly greater than zero.

-  ``reference displacement vector`` is the prescribed displacement vector that defines the nominal direction of photon propagation. If the ``photon maximum angular offset`` is zero, photons travel exactly in this direction.

-  ``photon insertion maximum offset`` is the maximum positional offset applied randomly to the initial location of each photon in each direction defined by the ``insertion unit tensors``. If set to zero, photons are perfectly aligned.

-  ``photon insertion prn seed`` is the pseudo random seed used to generate the insertion offsets.

-  ``photon maximum angular offset`` is the maximum angular deviation allowed between each photon’s propagation direction and the reference displacement vector. A value of zero means all photons move exactly along the reference vector otherwise, photons are scattered randomly within the specified angle defined in radians.

-  ``photon angular offset prn seed`` is the pseudo random seed used to generate the angular offsets.
