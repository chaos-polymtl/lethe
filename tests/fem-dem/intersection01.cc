#include <fem-dem/particle_projector.h>

int
main()
{
  Point<3> c_sphere(0,0,0);
  double r_sphere = 0.5;
  Tensor<1,3> normal_vector;
  normal_vector[0] = 0.0; normal_vector[1] = 0.0; normal_vector[2] = 1.0;
  Point<3> p_plane(0,0,0.4);

  double volume = plane_sphere_intersection (c_sphere, r_sphere, normal_vector, p_plane);

  std::cout << "Intersection volume of 3D sphere of center (" << c_sphere << ") and radius " << r_sphere << " with the plane going through point (" << p_plane << ") and normal to the vector (" << normal_vector << ") : " << volume << std::endl;

Point<2> c_sphere_2(0,0);
  Tensor<1,2> normal_vector_2;
  normal_vector_2[0] = 0.0; normal_vector_2[1] = 1.0;
  Point<2> p_plane_2(0,0.4);

volume = plane_sphere_intersection (c_sphere_2, r_sphere, normal_vector_2, p_plane_2);

 std::cout << "Intersection area of 2D circle of center (" << c_sphere_2 << ") and radius " << r_sphere << " with the line going through point (" << p_plane_2 << ") and normal to the vector (" << normal_vector_2 << ") : " << volume << std::endl;
}
