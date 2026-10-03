// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_time_harmonic_maxwell_scratch_data_h
#define lethe_time_harmonic_maxwell_scratch_data_h

#include <core/physical_property_model.h>

#include <solvers/physical_properties_manager.h>
#include <solvers/physics_scratch_data.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/table.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/types.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping.h>

#include <array>
#include <complex>
#include <map>
#include <memory>
#include <vector>

DeclException1(
  TimeHarmonicMaxwellDimensionNotSupported,
  int,
  << "The time-harmonic Maxwell solver does not support dimension: " << arg1
  << ". Currently, only 3D problems are supported as the 2D version of curls and cross products have completely different definitions than their 3D counterparts.");

/**
 * This helper function projects a 3D tensor onto the tangential plane
 * defined by the given normal vector. Mathematically, it computes:
 * \f[
 * \mathbf{t} = \mathbf{n} \times (\mathbf{v} \times \mathbf{n})
 * \f]
 * which removes the normal component of the vector \f$\mathbf{v}\f$,
 * keeping only the tangential part.
 *
 * @tparam dim Spatial dimension.
 * @tparam Number Type of the components of the vector (real or complex).
 * @param tensor Input vector (field value) to be projected.
 * @param normal Unit normal vector defining the face orientation.
 * @return The tangential component of the input vector.
 *
 * @note This operation is used to obtain traces in H^{-1/2}(curl) spaces,
 * where only tangential components are continuous across interfaces. The
 * inline keyword is enforced to make sure that the tensor operations are
 * efficient as they are called frequently during assembly.
 */
template <int dim, typename Number>
DEAL_II_ALWAYS_INLINE inline Tensor<1, dim, Number>
map_H12(const Tensor<1, dim, Number> &tensor, const Tensor<1, dim> &normal)
{
  if (dim != 3)
    {
      AssertThrow(
        false,
        ExcMessage(
          "The map_H12 function is only implemented for 3D problems."));
    }
  return cross_product_3d(normal, cross_product_3d(tensor, normal));
}

/**
 * @brief Compute the effective electromagnetic properties of a material at a
 * set of points. The effective electric permittivity is
 * \f$\varepsilon' + i(\varepsilon'' + \sigma)\f$ and the effective magnetic
 * permeability is \f$\mu' + i\mu''\f$, where each property is evaluated with
 * the physical property models of the material.
 *
 * @note The time-harmonic Maxwell equations do not support multiple fluids, so
 * the properties only depend on the material id of the cell since the interface
 * between two materials is conformed to the mesh.
 *
 *  @param[in] physical_properties_manager The object that manages the
 * physical properties of the problem and provides them at any given position
 * of the domain.
 * @param[in] field_values A map containing the values of the fields at the
 * current position. This is used to compute the material properties that
 * depend on the fields.
 *  @param[in] material_id The material id of the current position, used to
 * determine the appropriate material properties from the input parameters.
 *  @param[out] effective_electric_permittivity Effective electric
 * permittivity at the current position. This value is updated in place by the
 * function.
 *  @param[out] effective_magnetic_permeability Effective magnetic
 * permeability at the current position. This value is updated in place by the
 * function.
 */
void
compute_effective_electromagnetic_properties(
  const PhysicalPropertiesManager            &physical_properties_manager,
  const std::map<field, std::vector<double>> &field_values_vectors,
  const unsigned int                          material_id,
  std::vector<std::complex<double>>          &effective_electric_permittivities,
  std::vector<std::complex<double>> &effective_magnetic_permeabilities);


/**
 * @brief Class that stores the information required by the assembly of the
 * Discontinuous Petrov-Galerkin (DPG) ultraweak formulation of the
 * time-harmonic Maxwell equations.
 *
 * The three finite element spaces (test, interior trial and skeleton trial)
 * are [E_real, E_imag, H_real, H_imag] systems whose four fields use the same
 * base element. Every shape function is therefore a real base function
 * \f$\boldsymbol{\varphi}_a\f$ of the base element multiplied by 1 (real part)
 * or \f$i\f$ (imaginary part), and the same base functions describe the
 * electric and magnetic fields. Consequently, the scratch only evaluates the
 * real base functions of each space (and the curl of the test base functions)
 * at the quadrature points. The assemblers combine them with the complex
 * coefficients of the equations and fill the real 2x2 blocks associated with
 * the real and imaginary parts of each pair of base functions.
 *
 * To avoid querying the finite elements during the assembly, each shape
 * function of the three spaces is classified once, at construction, with the
 * ShapeFunctionType flags, and the local dof index of each base function is
 * stored for each of the four fields. For each face of the reference cell, the
 * scratch also stores the base functions associated with the face, which are
 * the only ones with a non-zero tangential trace on the face.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions. Only
 * dim = 3 is supported.
 *
 * @ingroup solvers
 */
template <int dim>
class TimeHarmonicMaxwellScratchData : public PhysicsScratchDataBase
{
public:
  /**
   * @brief Bit flags identifying to which field a shape function of the
   * [E_real, E_imag, H_real, H_imag] finite element systems belongs. The
   * composite flags is_electric and is_magnetic allow to check whether a shape
   * function belongs to the electric or magnetic field, regardless of whether
   * it is the real or imaginary part.
   */
  enum ShapeFunctionType : unsigned char
  {
    electric_real = 1u << 0,
    electric_imag = 1u << 1,
    magnetic_real = 1u << 2,
    magnetic_imag = 1u << 3,

    is_electric = electric_real | electric_imag,
    is_magnetic = magnetic_real | magnetic_imag
  };

  /**
   * Index of each field in the arrays of local dof indices of the base
   * functions (e.g. test_base_dofs[E_real_field][a] is the local dof index of
   * the base function a for the real part of the electric field).
   */
  static constexpr unsigned int E_real_field = 0;
  static constexpr unsigned int E_imag_field = 1;
  static constexpr unsigned int H_real_field = 2;
  static constexpr unsigned int H_imag_field = 3;

  /**
   * @brief Constructor. Creates the FEValues and FEFaceValues objects of the
   * interior trial, skeleton trial and test spaces and allocates the memory of
   * the scratch. It does not do any evaluation, since this is done at the cell
   * level.
   *
   * @param[in] properties_manager Manager of the physical property models.
   *
   * @param[in] fe_trial_interior Finite element of the interior trial space.
   *
   * @param[in] fe_trial_skeleton Finite element of the skeleton trial space.
   *
   * @param[in] fe_test Finite element of the test space.
   *
   * @param[in] cell_quadrature Quadrature used for the cell integrals.
   *
   * @param[in] face_quadrature Quadrature used for the face integrals.
   *
   * @param[in] mapping Mapping of the domain.
   */
  TimeHarmonicMaxwellScratchData(
    const PhysicalPropertiesManager &properties_manager,
    const FiniteElement<dim>        &fe_trial_interior,
    const FiniteElement<dim>        &fe_trial_skeleton,
    const FiniteElement<dim>        &fe_test,
    const Quadrature<dim>           &cell_quadrature,
    const Quadrature<dim - 1>       &face_quadrature,
    const Mapping<dim>              &mapping)
    : properties_manager(properties_manager)
    , fe_values_trial_interior(mapping,
                               fe_trial_interior,
                               cell_quadrature,
                               update_values | update_quadrature_points |
                                 update_JxW_values)
    , fe_values_test(mapping,
                     fe_test,
                     cell_quadrature,
                     update_values | update_gradients)
    , fe_face_values_trial_skeleton(mapping,
                                    fe_trial_skeleton,
                                    face_quadrature,
                                    update_values | update_quadrature_points |
                                      update_normal_vectors | update_JxW_values)
    , fe_face_values_test(mapping, fe_test, face_quadrature, update_values)
    , gather_temperature(false)
  {
    allocate();
  }

  /**
   * @brief Copy constructor. Same as the main constructor. It only uses the
   * other scratch to build the FEValues objects, it does not copy the content
   * of the other scratch since, by definition of the WorkStream mechanism, the
   * content of the scratch is reset on a cell basis.
   *
   * @param[in] sd The scratch data to be copied.
   */
  TimeHarmonicMaxwellScratchData(const TimeHarmonicMaxwellScratchData<dim> &sd)
    : properties_manager(sd.properties_manager)
    , fe_values_trial_interior(sd.fe_values_trial_interior.get_mapping(),
                               sd.fe_values_trial_interior.get_fe(),
                               sd.fe_values_trial_interior.get_quadrature(),
                               update_values | update_quadrature_points |
                                 update_JxW_values)
    , fe_values_test(sd.fe_values_test.get_mapping(),
                     sd.fe_values_test.get_fe(),
                     sd.fe_values_test.get_quadrature(),
                     update_values | update_gradients)
    , fe_face_values_trial_skeleton(
        sd.fe_face_values_trial_skeleton.get_mapping(),
        sd.fe_face_values_trial_skeleton.get_fe(),
        sd.fe_face_values_trial_skeleton.get_quadrature(),
        update_values | update_quadrature_points | update_normal_vectors |
          update_JxW_values)
    , fe_face_values_test(sd.fe_face_values_test.get_mapping(),
                          sd.fe_face_values_test.get_fe(),
                          sd.fe_face_values_test.get_quadrature(),
                          update_values)
    , gather_temperature(false)
  {
    allocate();
    if (sd.gather_temperature)
      enable_temperature(sd.fe_values_temperature->get_fe(),
                         sd.fe_values_temperature->get_quadrature(),
                         sd.fe_face_values_temperature->get_quadrature(),
                         sd.fe_values_temperature->get_mapping());
  }

  /**
   * @brief Allocate the memory of all the members of the scratch. It also
   * classifies the shape functions of the three finite element spaces with the
   * ShapeFunctionType flags, stores the local dof indices of their base
   * functions for each field and builds, for each face of the reference cell,
   * the lists of base functions associated with the face.
   */
  void
  allocate() override;

  /**
   * @brief Classify the shape functions of a finite element space of the
   * time-harmonic Maxwell equations according to the field they belong to and
   * store the local dof index of each base function for each field. Each shape
   * function of the [E_real, E_imag, H_real, H_imag] finite element system
   * receives the ShapeFunctionType flag of its field, which is determined with
   * FiniteElement::shape_function_belongs_to(). Its base function index is
   * given by its position in the base element (and the copy of the base
   * element for vector-valued fields made of scalar elements, such as FE_DGQ),
   * so the same index refers to the same base function in the four fields.
   *
   * @param[in] fe Finite element system of the [E_real, E_imag, H_real,
   * H_imag] fields (interior trial, skeleton trial or test space).
   *
   * @param[out] shape_function_type ShapeFunctionType flag of each shape
   * function of the finite element.
   *
   * @param[out] base_dofs Local dof index of each base function for each
   * field: base_dofs[field][a].
   *
   * @param[out] base_function_of_dof Base function index of each shape
   * function of the finite element.
   */
  void
  classify_shape_functions(
    const FiniteElement<dim>                 &fe,
    std::vector<unsigned char>               &shape_function_type,
    std::array<std::vector<unsigned int>, 4> &base_dofs,
    std::vector<unsigned int>                &base_function_of_dof) const;

  /**
   * @brief Gather, for each face of the reference cell, the base functions of
   * a finite element space that are associated with the face (i.e., whose
   * degrees of freedom are located on the edges or the interior of the face).
   * For Nedelec elements, these are the only base functions whose tangential
   * trace does not vanish on the face.
   *
   * @param[in] fe Finite element system of the [E_real, E_imag, H_real,
   * H_imag] fields.
   *
   * @param[in] base_function_of_dof Base function index of each shape function
   * of the finite element (see classify_shape_functions()).
   *
   * @param[out] face_base_functions Base functions associated with each face:
   * face_base_functions[face_no].
   */
  void
  gather_face_base_functions(
    const FiniteElement<dim>               &fe,
    const std::vector<unsigned int>        &base_function_of_dof,
    std::vector<std::vector<unsigned int>> &face_base_functions) const;

  /**
   * @brief Enable the evaluation of the temperature field, on which the
   * electromagnetic properties may depend.
   *
   * @param[in] fe_temperature Finite element of the heat transfer physics.
   *
   * @param[in] cell_quadrature Quadrature used for the cell integrals.
   *
   * @param[in] face_quadrature Quadrature used for the face integrals.
   *
   * @param[in] mapping Mapping of the domain.
   */
  void
  enable_temperature(const FiniteElement<dim>  &fe_temperature,
                     const Quadrature<dim>     &cell_quadrature,
                     const Quadrature<dim - 1> &face_quadrature,
                     const Mapping<dim>        &mapping);

  /**
   * @brief Reinitialize the content of the scratch for a cell. It evaluates
   * the real base functions of the test space and their curl, and the real
   * base functions of the interior trial space at the quadrature points. The
   * temperature values are set to zero and are only overwritten by
   * reinit_temperature when the material of the cell depends on the
   * temperature.
   *
   * @param[in] cell The cell of the interior trial space DoFHandler over which
   * the assembly is carried out.
   *
   * @param[in] cell_test The same cell, but for the test space DoFHandler.
   */
  void
  reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
         const typename DoFHandler<dim>::active_cell_iterator &cell_test);

  /**
   * @brief Reinitialize the content of the scratch for a face of the current
   * cell. It evaluates the real base functions of the test space and the
   * tangential traces of the real base functions of the skeleton trial space
   * at the face quadrature points. The skeleton base functions are only
   * evaluated for the base functions associated with the face. The test base
   * functions are evaluated for all the base functions on boundary faces,
   * where the Robin boundary conditions need their full trace, and only for
   * the base functions associated with the face on interior faces. The face
   * temperature values are set to zero and are only overwritten by
   * reinit_face_temperature when the material of the cell depends on the
   * temperature.
   *
   * @param[in] cell_skeleton The current cell, for the skeleton trial space
   * DoFHandler.
   *
   * @param[in] cell_test The current cell, for the test space DoFHandler.
   *
   * @param[in] face_no Index of the face within the cell.
   */
  void
  reinit_face(
    const typename DoFHandler<dim>::active_cell_iterator &cell_skeleton,
    const typename DoFHandler<dim>::active_cell_iterator &cell_test,
    const unsigned int                                    face_no);

  /**
   * @brief Evaluate the temperature field at the cell quadrature points.
   *
   * @tparam VectorType The vector type of the temperature solution.
   *
   * @param[in] cell_temperature The current cell, for the heat transfer
   * DoFHandler.
   *
   * @param[in] temperature_solution Temperature solution of the heat transfer
   * physics.
   */
  template <typename VectorType>
  void
  reinit_temperature(
    const typename DoFHandler<dim>::active_cell_iterator &cell_temperature,
    const VectorType                                     &temperature_solution)
  {
    fe_values_temperature->reinit(cell_temperature);
    fe_values_temperature->get_function_values(temperature_solution,
                                               this->temperature_values);
  }

  /**
   * @brief Evaluate the temperature field at the quadrature points of a face
   * of the current cell.
   *
   * @tparam VectorType The vector type of the temperature solution.
   *
   * @param[in] cell_temperature The current cell, for the heat transfer
   * DoFHandler.
   *
   * @param[in] face_no Index of the face within the cell.
   *
   * @param[in] temperature_solution Temperature solution of the heat transfer
   * physics.
   */
  template <typename VectorType>
  void
  reinit_face_temperature(
    const typename DoFHandler<dim>::active_cell_iterator &cell_temperature,
    const unsigned int                                    face_no,
    const VectorType                                     &temperature_solution)
  {
    fe_face_values_temperature->reinit(cell_temperature, face_no);
    fe_face_values_temperature->get_function_values(
      temperature_solution, this->face_temperature_values);
  }

  /**
   * @brief Compute the effective electromagnetic properties at the cell
   * quadrature points from the fields gathered by the scratch.
   */
  void
  calculate_physical_properties();

  /**
   * @brief Compute the effective electromagnetic properties at the face
   * quadrature points from the fields gathered by the scratch.
   */
  void
  calculate_face_physical_properties();

  // Physical properties
  const PhysicalPropertiesManager properties_manager;
  unsigned int                    material_id;
  bool                            cell_material_needs_temperature;

  // FEValues of the three DPG spaces. Because all spaces share the same
  // triangulation and mapping, the quadrature points, normal vectors and JxW
  // values are only computed with the trial spaces.
  FEValues<dim>     fe_values_trial_interior;
  FEValues<dim>     fe_values_test;
  FEFaceValues<dim> fe_face_values_trial_skeleton;
  FEFaceValues<dim> fe_face_values_test;

  // Extractors of the [E_real, E_imag, H_real, H_imag] finite element systems
  const FEValuesExtractors::Vector extractor_E_real{0};
  const FEValuesExtractors::Vector extractor_E_imag{dim};
  const FEValuesExtractors::Vector extractor_H_real{2 * dim};
  const FEValuesExtractors::Vector extractor_H_imag{3 * dim};

  // Sizes
  unsigned int n_q_points;
  unsigned int n_face_q_points;
  unsigned int n_dofs_test;
  unsigned int n_dofs_trial_interior;
  unsigned int n_dofs_trial_skeleton;

  // Number of base functions of each space. Each field of a space has this
  // number of shape functions.
  unsigned int n_base_functions_test;
  unsigned int n_base_functions_trial_interior;
  unsigned int n_base_functions_trial_skeleton;

  // Classification of the shape functions of each space with the
  // ShapeFunctionType flags
  std::vector<unsigned char> shape_function_type_test;
  std::vector<unsigned char> shape_function_type_trial_interior;
  std::vector<unsigned char> shape_function_type_trial_skeleton;

  // Local dof index of each base function for each field:
  // test_base_dofs[field][a]. The electric test functions (F) are the E_real
  // and E_imag fields of the test space and the magnetic ones (I) are the
  // H_real and H_imag fields.
  std::array<std::vector<unsigned int>, 4> test_base_dofs;
  std::array<std::vector<unsigned int>, 4> trial_interior_base_dofs;
  std::array<std::vector<unsigned int>, 4> trial_skeleton_base_dofs;

  // Base functions associated with each face of the reference cell:
  // face_test_base_functions[face_no]
  std::vector<std::vector<unsigned int>> face_test_base_functions;
  std::vector<std::vector<unsigned int>> face_trial_skeleton_base_functions;

  // List of all the test base functions (0, 1, ..., n_base_functions_test - 1)
  std::vector<unsigned int> all_test_base_functions;

  // Cell quadrature
  std::vector<double>     JxW;
  std::vector<Point<dim>> quadrature_points;

  // Temperature and effective properties at the cell quadrature points
  std::vector<double>                  temperature_values;
  std::map<field, std::vector<double>> fields;
  std::vector<std::complex<double>>    effective_electric_permittivities;
  std::vector<std::complex<double>>    effective_magnetic_permeabilities;

  // Real base functions of the test space and their curl, and real base
  // functions of the interior trial space at the cell quadrature points,
  // indexed [a][q] so that the loops on the quadrature points are contiguous
  Table<2, Tensor<1, dim>> phi_test;
  Table<2, Tensor<1, dim>> curl_phi_test;
  Table<2, Tensor<1, dim>> phi_trial_interior;

  // Current face
  unsigned int       face_no;
  bool               face_at_boundary;
  types::boundary_id face_boundary_id;

  // Face quadrature
  std::vector<double>         face_JxW;
  std::vector<Point<dim>>     face_quadrature_points;
  std::vector<Tensor<1, dim>> face_normals;

  // Temperature and effective properties at the face quadrature points
  std::vector<double>                  face_temperature_values;
  std::map<field, std::vector<double>> face_fields;
  std::vector<std::complex<double>>    face_effective_electric_permittivities;
  std::vector<std::complex<double>>    face_effective_magnetic_permeabilities;

  // Real base functions of the test space and their cross product with the
  // normal at the face quadrature points, indexed [a][q]
  Table<2, Tensor<1, dim>> phi_test_face;
  Table<2, Tensor<1, dim>> n_cross_phi_test_face;

  // Tangential traces (map_H12) of the real base functions of the skeleton
  // trial space and their cross product with the normal at the face
  // quadrature points, indexed [d][q]
  Table<2, Tensor<1, dim>> tangential_phi_trial_skeleton;
  Table<2, Tensor<1, dim>> n_cross_tangential_phi_trial_skeleton;

  // Temperature coupling
  bool                               gather_temperature;
  std::shared_ptr<FEValues<dim>>     fe_values_temperature;
  std::shared_ptr<FEFaceValues<dim>> fe_face_values_temperature;
};

#endif
