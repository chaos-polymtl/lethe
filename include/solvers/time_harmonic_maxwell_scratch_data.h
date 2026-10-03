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
 * @param tensor Input vector (field value) to be projected.
 * @param normal Unit normal vector defining the face orientation.
 * @return The tangential component of the input vector.
 *
 * @note This operation is used to obtain traces in H^{-1/2}(curl) spaces,
 * where only tangential components are continuous across interfaces. The
 * inline keyword is enforced to make sure that the tensor operations are
 * efficient as they are called frequently during assembly.
 */
template <int dim>
DEAL_II_ALWAYS_INLINE inline Tensor<1, dim, std::complex<double>>
map_H12(const Tensor<1, dim, std::complex<double>> &tensor,
        const Tensor<1, dim>                       &normal)
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
 * time-harmonic Maxwell equations. For each cell, it evaluates the complex
 * shape functions of the test space (\f$\mathbf{F}\f$ for the electric test
 * functions and \f$\mathbf{I}\f$ for the magnetic ones, with their curls and
 * complex conjugates) and of the interior trial space (\f$\mathbf{E}\f$ and
 * \f$\mathbf{H}\f$) and the effective material properties at the quadrature
 * points. For each face of the cell, it evaluates
 * the test functions and the tangential traces of the skeleton trial space
 * (\f$\hat{\mathbf{E}}\f$ and \f$\hat{\mathbf{H}}\f$) at the face quadrature
 * points. The complex values are stored with their complex conjugates because
 * the conjugate of a complex tensor is not implemented in deal.II.
 *
 * To avoid querying the finite elements during the assembly, each shape
 * function of the three spaces is classified once, at construction, with the
 * ShapeFunctionType flags, and the shape functions of each field are gathered
 * in lists of local dof indices. The shape functions of a field are only
 * evaluated for the dofs of that field, since they vanish for the others.
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
   * ShapeFunctionType flags and builds the lists of local dof indices of each
   * field.
   */
  void
  allocate() override;

  /**
   * @brief Classify the shape functions of a finite element space of the
   * time-harmonic Maxwell equations according to the field they belong to.
   * Each shape function of the [E_real, E_imag, H_real, H_imag] finite element
   * system receives the ShapeFunctionType flag of its field, which is
   * determined with FiniteElement::shape_function_belongs_to(). Since every
   * shape function belongs to exactly one field, the local indices of the
   * electric and magnetic shape functions are then gathered in two lists.
   *
   * @param[in] fe Finite element system of the [E_real, E_imag, H_real,
   * H_imag] fields (interior trial, skeleton trial or test space).
   *
   * @param[out] shape_function_type ShapeFunctionType flag of each shape
   * function of the finite element.
   *
   * @param[out] electric_dofs Local indices of the shape functions that belong
   * to the electric field (E_real or E_imag).
   *
   * @param[out] magnetic_dofs Local indices of the shape functions that belong
   * to the magnetic field (H_real or H_imag).
   */
  void
  classify_shape_functions(const FiniteElement<dim>   &fe,
                           std::vector<unsigned char> &shape_function_type,
                           std::vector<unsigned int>  &electric_dofs,
                           std::vector<unsigned int>  &magnetic_dofs) const;

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
   * the shape functions of the test and interior trial spaces at the
   * quadrature points. The temperature values are set to zero and are only
   * overwritten by reinit_temperature when the material of the cell depends on
   * the temperature.
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
   * cell. It evaluates the test functions and the tangential traces of the
   * skeleton trial space at the face quadrature points. The face temperature
   * values are set to zero and are only overwritten by reinit_face_temperature
   * when the material of the cell depends on the temperature.
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

  // Classification of the shape functions of each space with the
  // ShapeFunctionType flags
  std::vector<unsigned char> shape_function_type_test;
  std::vector<unsigned char> shape_function_type_trial_interior;
  std::vector<unsigned char> shape_function_type_trial_skeleton;

  // Local dof indices of each field: electric (F) and magnetic (I) test
  // functions, interior electric (E) and magnetic (H) trial functions and
  // skeleton electric (E_hat) and magnetic (H_hat) trial functions
  std::vector<unsigned int> test_dofs_electric;
  std::vector<unsigned int> test_dofs_magnetic;
  std::vector<unsigned int> trial_interior_dofs_electric;
  std::vector<unsigned int> trial_interior_dofs_magnetic;
  std::vector<unsigned int> trial_skeleton_dofs_electric;
  std::vector<unsigned int> trial_skeleton_dofs_magnetic;

  // Cell quadrature
  std::vector<double>     JxW;
  std::vector<Point<dim>> quadrature_points;

  // Temperature and effective properties at the cell quadrature points
  std::vector<double>                  temperature_values;
  std::map<field, std::vector<double>> fields;
  std::vector<std::complex<double>>    effective_electric_permittivities;
  std::vector<std::complex<double>>    effective_magnetic_permeabilities;

  // Test functions at the cell quadrature points, indexed [q][k]
  Table<2, Tensor<1, dim, std::complex<double>>> phi_F;
  Table<2, Tensor<1, dim, std::complex<double>>> phi_F_conj;
  Table<2, Tensor<1, dim, std::complex<double>>> curl_phi_F;
  Table<2, Tensor<1, dim, std::complex<double>>> curl_phi_F_conj;
  Table<2, Tensor<1, dim, std::complex<double>>> phi_I;
  Table<2, Tensor<1, dim, std::complex<double>>> phi_I_conj;
  Table<2, Tensor<1, dim, std::complex<double>>> curl_phi_I;
  Table<2, Tensor<1, dim, std::complex<double>>> curl_phi_I_conj;

  // Interior trial functions at the cell quadrature points, indexed [q][k]
  Table<2, Tensor<1, dim, std::complex<double>>> phi_E;
  Table<2, Tensor<1, dim, std::complex<double>>> phi_H;

  // Current face
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

  // Test functions at the face quadrature points, indexed [q][k]
  Table<2, Tensor<1, dim, std::complex<double>>> phi_F_face;
  Table<2, Tensor<1, dim, std::complex<double>>> phi_F_face_conj;
  Table<2, Tensor<1, dim, std::complex<double>>> phi_I_face_conj;
  Table<2, Tensor<1, dim, std::complex<double>>> n_cross_phi_I_face;
  Table<2, Tensor<1, dim, std::complex<double>>> n_cross_phi_I_face_conj;

  // Tangential traces of the skeleton trial functions at the face quadrature
  // points, indexed [q][k]
  Table<2, Tensor<1, dim, std::complex<double>>> phi_E_hat;
  Table<2, Tensor<1, dim, std::complex<double>>> n_cross_phi_E_hat;
  Table<2, Tensor<1, dim, std::complex<double>>> n_cross_phi_H_hat;

  // Temperature coupling
  bool                               gather_temperature;
  std::shared_ptr<FEValues<dim>>     fe_values_temperature;
  std::shared_ptr<FEFaceValues<dim>> fe_face_values_temperature;
};

#endif
