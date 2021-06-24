#include <deal.II/base/function.h>
using namespace dealii;

template <int dim>
class RotatingWall : public Function<dim>
{
public:
  RotatingWall () : Function<dim>(dim+1) {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component ) const;
};


template <int dim>
double RotatingWall<dim>::value (const Point<dim> &p,
                                   const unsigned int component) const
{
    Assert (component < this->n_components,
            ExcIndexRange (component, 0, this->n_components))

    if (component==0)
        return -p[1];
    else if(component==1)
        return p[0];
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
