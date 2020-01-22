#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>

using namespace dealii;

template<int dim>
class ExactSolutionMMS : public Function<dim>
{
public:
    ExactSolutionMMS() : Function<dim>(3) {}
    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &values) const;
};
template<int dim>
void ExactSolutionMMS<dim>::vector_value(const Point<dim> &p,
                                                    Vector<double> &values) const
{
    const double a = M_PI;
    double x = p[0];
    double y = p[1];
    values(0) = sin(a*x)*sin(a*x)*cos(a*y)*sin(a*y);
    values(1) = -cos(a*x)*sin(a*x)*sin(a*y)*sin(a*y);
    values(2) = -2 + x*x + y*y;
}

template<int dim>
class ExactSolutionTaylorCouette : public Function<dim>
{
public:
    ExactSolutionTaylorCouette() : Function<dim>(3)
    {
        eta_=0.25;
        ri_=0.25;
    }
    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &values) const;

private:
    double eta_;
    double ri_=0.25;
};
template<int dim>
void ExactSolutionTaylorCouette<dim>::vector_value(const Point<dim> &p,
                                                    Vector<double> &values) const
{
    const double a = M_PI;
    double x = p[0];
    double y = p[1];
    double r= std::sqrt(x*x+y*y);
    double theta= std::atan2(y,x);
    double A= -(eta_*eta_)/(1.-eta_*eta_);
    double B= ri_ * ri_ / (1.-eta_*eta_);
    double utheta= A*r + B/r;
    values(0) = -std::sin(theta)*utheta;
    values(1) = std::cos(theta)*utheta;
    values(2) = 0.;
}
