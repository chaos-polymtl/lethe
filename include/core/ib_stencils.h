//
// Created by lucka on 2021-05-21.
//


#include <deal.II/base/table_handler.h>
#include <deal.II/base/tensor.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping_manifold.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>

#include <core/ib_particle.h>


#ifndef LETHE_IB_STENCILS_H
#define LETHE_IB_STENCILS_H



using namespace dealii;



template<int dim>
class IBStencilsS1 {

public:

    virtual unsigned int nb_points();

    virtual std::tuple<Point<dim>,std::vector<Point<dim>>> points(IBParticle<dim> p,Point<dim> dof_point);

    virtual std::vector<double> coefficients();




};


template<int dim>
class IBStencilsS2 {

public:

    virtual unsigned int nb_points();

    virtual std::tuple<Point<dim>,std::vector<Point<dim>>> points(IBParticle<dim> p,Point<dim> dof_point);

    virtual std::vector<double> coefficients();




};

template<int dim>
class IBStencilsS3{

public:

    virtual unsigned int nb_points();

    virtual std::tuple<Point<dim>,std::vector<Point<dim>>> points(IBParticle<dim> p,Point<dim> dof_point);

    virtual std::vector<double> coefficients();




};

template<int dim>
class IBStencilsS4{

public:

    virtual unsigned int nb_points();

    virtual std::tuple<Point<dim>,std::vector<Point<dim>>> points(IBParticle<dim> p,Point<dim> dof_point);

    virtual std::vector<double> coefficients();




};

template<int dim>
class IBStencilsCell{

public:

    virtual unsigned int nb_points();

    virtual std::tuple<Point<dim>,std::vector<Point<dim>>> points(IBParticle<dim> p,Point<dim> dof_point);

    virtual std::vector<double> coefficients();





};

#endif //LETHE_IB_STENCILS_H