#include <deal.II/base/function.h>
#include <deal.II/base/parsed_function.h>

using namespace dealii;

template <int dim>
class FunctionDefined : public Function<dim>
{
private:
  Functions::ParsedFunction<dim> *u;
  Functions::ParsedFunction<dim> *v;
  Functions::ParsedFunction<dim> *w;

public:
  FunctionDefined (Functions::ParsedFunction<dim> *p_u, Functions::ParsedFunction<dim> *p_v,Functions::ParsedFunction<dim> *p_w) :
    Function<dim>(dim+1),
    u(p_u),
    v(p_v),
    w(p_w)
  {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component ) const;
};

template <int dim>
double FunctionDefined<dim>::value (const Point<dim> &p,
                                      const unsigned int component) const
{
  Assert (component < this->n_components,
          ExcIndexRange (component, 0, this->n_components));
  if (component==0)
    {
      return u->value(p);
    }
  else if(component==1)
    {
      return v->value(p);
    }
  else if(component==2)
    {
      return w->value(p);
    }
  return 0.;
}


template <int dim>
class PoiseuilleInlet : public Function<dim>
{
public:
  PoiseuilleInlet () : Function<dim>(dim+1)
  {
      y2_=1.;
      y1_=0.;
      dy_= 0.5*(y2_ + y1_);
      vmax_=1./dy_/dy_;
  };

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component ) const;
private:
  double y2_;
  double y1_;
  double dy_;
  double vmax_;
};


template <int dim>
double PoiseuilleInlet<dim>::value (const Point<dim> &p,
                                    const unsigned int component) const
{
    Assert (component < this->n_components,
            ExcIndexRange (component, 0, this->n_components));

    double y=p[1];

    if (component==0)
    {
        return vmax_ *(y-y1_)*(y2_-y);
    }
    else if(component==1)
        return 0.;
    return 0.;
}

template <int dim>
class ConstantXInlet : public Function<dim>
{
public:
  ConstantXInlet () : Function<dim>(dim+1)
  {
      value_=1.;
  };

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component ) const;
private:
    double value_;
};


template <int dim>
double ConstantXInlet<dim>::value (const Point<dim> &/*p*/,
                                    const unsigned int component) const
{
    Assert (component < this->n_components,
            ExcIndexRange (component, 0, this->n_components));

    if (component==0)
    {
        return value_;
    }
    else if(component==1)
        return 0.;
    return 0.;
}

template <int dim>
class ConstantXSlip : public Function<dim>
{
public:
  ConstantXSlip () : Function<dim>(dim+1)
  {
      value_=1.;
  };

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component ) const;
private:
    double value_;
};


template <int dim>
double ConstantXSlip<dim>::value (const Point<dim> &p,
                                    const unsigned int component) const
{
    Assert (component < this->n_components,
            ExcIndexRange (component, 0, this->n_components));

    if (component==0)
    {
        return value_;
    }
    else if(component==1)
        return 0.;
    return 0.;
}
