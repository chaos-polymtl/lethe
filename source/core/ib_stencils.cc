//
// Created by lucka on 2021-05-21.
//

#include <core/ib_stencils.h>





template<int dim>
unsigned int
IBStencilsS1<dim>::nb_points()
{
    return 1;
}

template<int dim>
std::vector<double>
IBStencilsS1<dim>::coefficients()
{

    std::vector<double> coef(2);
    coef[0] = 9;
    coef[1] = -8;

    return coef;

}

template<int dim>
std::tuple<Point<dim>,std::vector<Point<dim>>>
IBStencilsS1<dim>::points(IBParticle<dim> p,Point<dim> dof_point)
{

    Tensor<1, dim, double> vect_ib = (dof_point - p.position -
                                          p.radius * (dof_point - p.position) / (dof_point - p.position).norm());


    Point<dim, double> ib_point(dof_point - vect_ib);

    Point<dim, double> interpolation_point_1(dof_point + vect_ib * 1 / 8);

    std::vector<Point < dim>>interpolation_points(1);
    interpolation_points[0]=interpolation_point_1;

    return {ib_point,interpolation_points};

}


template <int dim>
unsigned int
IBStencilsS2<dim>::nb_points()
{
    return 2;
}

template <int dim>
std::vector<double>
IBStencilsS2<dim>::coefficients()
{
    std::vector<double> coef(3);
    coef[0]=153;;
    coef[1]=-288;
    coef[2]=136;


    return coef;
}

template <int dim>
std::tuple<Point<dim>,std::vector<Point<dim>>>
IBStencilsS2<dim>::points(IBParticle<dim> p,Point<dim> dof_point)
{
    Tensor<1, dim, double> vect_ib =(dof_point - p.position -p.radius *(dof_point-p.position)/(dof_point-p.position).norm());


    Point<dim, double> ib_point(dof_point - vect_ib);

    Point<dim, double> interpolation_point_1(dof_point+vect_ib* 1./16.);

    Point<dim, double> interpolation_point_2(dof_point+vect_ib* 1./8.);

    std::vector<Point<dim>> interpolation_points(2);
    interpolation_points[0]=interpolation_point_1;
    interpolation_points[1]=interpolation_point_2;

    return {ib_point,interpolation_points};
}

template <int dim>
unsigned int
IBStencilsS3<dim>::nb_points()
{
    return 3;
}

template <int dim>
std::vector<double>
IBStencilsS3<dim>::coefficients()
{
    std::vector<double> coef(4);
    coef[0]=2925;
    coef[1]=-8424;
    coef[2]=8100;
    coef[3]=-2600;


    return coef;
}

template <int dim>
std::tuple<Point<dim>,std::vector<Point<dim>>>
IBStencilsS3<dim>::points(IBParticle<dim> p,Point<dim> dof_point)
{
    Tensor<1, dim, double> vect_ib =(dof_point - p.position -p.radius *(dof_point-p.position)/(dof_point-p.position).norm());


    Point<dim, double> ib_point(dof_point - vect_ib);

    Point<dim, double> interpolation_point_1(dof_point+vect_ib* 1./24.);

    Point<dim, double> interpolation_point_2(dof_point+vect_ib* 1/12.);

    Point<dim, double> interpolation_point_3(dof_point+vect_ib* 1./8.);

    std::vector<Point<dim>> interpolation_points(3);
    interpolation_points[0]=interpolation_point_1;
    interpolation_points[1]=interpolation_point_2;
    interpolation_points[2]=interpolation_point_3;

    return {ib_point,interpolation_points};
}



template <int dim>
unsigned int
IBStencilsS4<dim>::nb_points()
{
    return 4;
}

template <int dim>
std::vector<double>
IBStencilsS4<dim>::coefficients()
{
    std::vector<double> coef(5);
    coef[0]=58905;
    coef[1]=-228480;
    coef[2]=332640;
    coef[3]=-215424;
    coef[4]=52360;


    return coef;
}

template <int dim>
std::tuple<Point<dim>,std::vector<Point<dim>>>
IBStencilsS4<dim>::points(IBParticle<dim> p,Point<dim> dof_point)
{
    Tensor<1, dim, double> vect_ib =(dof_point - p.position -p.radius *(dof_point-p.position)/(dof_point-p.position).norm());


    Point<dim, double> ib_point(dof_point - vect_ib);

    Point<dim, double> interpolation_point_1(dof_point+vect_ib* 1./32.);

    Point<dim, double> interpolation_point_2(dof_point+vect_ib* 1/16.);

    Point<dim, double> interpolation_point_3(dof_point+vect_ib* 3./32.);

    Point<dim, double> interpolation_point_4(dof_point+vect_ib* 1./8.);

    std::vector<Point<dim>> interpolation_points(4);
    interpolation_points[0]=interpolation_point_1;
    interpolation_points[1]=interpolation_point_2;
    interpolation_points[2]=interpolation_point_3;
    interpolation_points[3]=interpolation_point_4;

    return {ib_point,interpolation_points};
}


template <int dim>
unsigned int
IBStencilsCell<dim>::nb_points()
{
    return 1;
}

template <int dim>
std::vector<double>
IBStencilsCell<dim>::coefficients()
{
    std::vector<double> coef(1);
    coef[0]=1;



    return coef;
}

template <int dim>
std::tuple<Point<dim>,std::vector<Point<dim>>>
IBStencilsCell<dim>::points(IBParticle<dim> p,Point<dim> dof_point)
{
    Tensor<1, dim, double> vect_ib =(dof_point - p.position -p.radius *(dof_point-p.position)/(dof_point-p.position).norm());


    Point<dim, double> ib_point(dof_point - vect_ib);

    Point<dim, double> interpolation_point_1(dof_point+vect_ib* 1./8.);

    std::vector<Point<dim>> interpolation_points(1);
    interpolation_points[0]=interpolation_point_1;

    return {ib_point,interpolation_points};
}







template class IBStencilsS1<2>;
template class IBStencilsS1<3>;
template class IBStencilsS2<2>;
template class IBStencilsS2<3>;
template class IBStencilsS3<2>;
template class IBStencilsS3<3>;
template class IBStencilsS4<2>;
template class IBStencilsS4<3>;
template class IBStencilsCell<2>;
template class IBStencilsCell<3>;