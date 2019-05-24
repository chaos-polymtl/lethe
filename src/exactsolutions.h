#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>

using namespace dealii;
#ifndef LETHE_EXACTSOLUTIONS_H
#define LETHE_EXACTSOLUTIONS_H

template<int dim>
class ExactSolutionMMS : public Function<dim>
{
public:
    ExactSolutionMMS() : Function<dim>(3) {}
    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &values) const;
private:
    double time;
};
template<int dim>
void ExactSolutionMMS<dim>::vector_value(const Point<dim> &p,
                                         Vector<double> &values) const
{
    assert(dim==2);
    const double a = M_PI;
    double x = p[0];
    double y = p[1];
    values(0) = sin(a*x)*sin(a*x)*cos(a*y)*sin(a*y);
    values(1) = -cos(a*x)*sin(a*x)*sin(a*y)*sin(a*y);
}


#endif




