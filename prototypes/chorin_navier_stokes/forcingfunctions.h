#include <deal.II/base/function.h>

// Finally, this is as in previous programs:
using namespace dealii;

template<int dim>
class MMSSineForcingFunction : public Function<dim>
{
public:
    MMSSineForcingFunction() : Function<dim>(3) {};
    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &values) const;
};
template<int dim>
void MMSSineForcingFunction<dim>::vector_value(const Point<dim> &p,
                                                     Vector<double> &values) const
{
    assert(dim==2);
    const double a = M_PI;

    double x = p[0];
    double y = p[1];
    values(0) = (2*a*a*(-sin(a*x)*sin(a*x) +
                       cos(a*x)*(cos(a*x)))*sin(a*y)*cos(a*y)
            - 4*a*a*sin(a*x)*sin(a*x)*sin(a*y)*cos(a*y)
                 - 2.0*x)*(-1.)
            + a*std::pow(sin(a*x),3.) * std::pow(sin(a*y),2.) * std::cos(a*x);
    values(1) = (2*a*a*(sin(a*y)*(sin(a*y)) - cos(a*y)*cos(a*y))
            *sin(a*x)*cos(a*x) + 4*a*a*sin(a*x)*sin(a*y)*sin(a*y)
            *cos(a*x) - 2.0*y)*(-1)
            + a*std::pow(sin(a*x),2.) * std::pow(sin(a*y),3.) * std::cos(a*y);


}

template<int dim>
class NoForce : public Function<dim>
{
public:
    NoForce() : Function<dim>(3) {};
    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &values) const;
};
template<int dim>
void NoForce<dim>::vector_value(const Point<dim> &/*p*/,
                                Vector<double> &values) const
{
    assert(dim==2);
    values(0) = 0.;
    values(1) = 0.;

}
