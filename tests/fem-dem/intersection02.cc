#include <fem-dem/particle_projector.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>


using namespace dealii;

int
main()
{
  const unsigned int fe_degree            = 1;
  const unsigned int mapping_degree       = 1;
  const unsigned int dim                  = 3;
  const unsigned int n_global_refinements = 2;
  Point<3> c_sphere(1,0,0);
  double r_sphere = 0.5;

  double V_sphere_out = 0;
  
  Triangulation<dim> tria;
  DoFHandler<dim> dof_handler(tria);
  const FE_Q<3>    fe(1);
  std::shared_ptr<dealii::Mapping<dim>> mapping;

  mapping = std::make_shared<dealii::MappingQ1<dim>>();

  GridGenerator::hyper_cube(tria);
  tria.refine_global(n_global_refinements);

  dof_handler.distribute_dofs(fe);

  std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;

  for (const auto &cell : dof_handler.cell_iterators())
      {
        if (cell-> at_boundary()){
          std::cout << "Volume of intersection with the boundary for cell with center at (" << cell->center() << ") is : " << sphere_boundary_intersection (mapping, cell, c_sphere, r_sphere);
          std::cout << std::endl;
          std::cout << "The face center as calculated is at : " << cell->face(0)->center() << std::endl;
        }   
      }
}
