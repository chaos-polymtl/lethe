// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_time_harmonic_maxwell_scratch_data_h
#define lethe_time_harmonic_maxwell_scratch_data_h

#include <core/physical_property_model.h>

#include <solvers/physical_properties_manager.h>
#include <solvers/physics_scratch_data.h>

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

/**
 * @brief Compute the effective electromagnetic properties of a material at a
 * set of points. The effective electric permittivity is
 * \f$\varepsilon' + i(\varepsilon'' + \sigma)\f$ and the effective magnetic
 * permeability is \f$\mu' + i\mu''\f$, where each property is evaluated with
 * the physical property models of the material. If the material is an
 * effective medium, the properties are given by the effective medium models of
 * its fluid-dispersed phase interaction. The conductivity of both phases is
 * then already included in the imaginary part of the effective permittivity.
 *
 * @note The time-harmonic Maxwell equations do not support multiple
 * fluids so the fluid_id is not necessary to determine the material
 * properties. The material properties only depend on the material_id, which
 * is used to identify the different materials in the input file and assign
 * them their corresponding properties.
 *
 *  @param[in] physical_properties_manager The object that manages the
 * physical properties of the problem and provides them at any given position
 * of the domain.
 * @param[in] field_values_vector A map containing the vector values of the
 * fields at the current position. This is used to compute the material
 * properties that depend on the fields.
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
 * @brief Check if at least one electromagnetic property of a material depends
 * on a field. If the material is an effective medium, the dependencies are
 * those of the effective medium models of its fluid-dispersed phase
 * interaction. Otherwise, they are those of the electric permittivity, electric
 * conductivity and magnetic permeability models of the material.
 *
 * @param[in] physical_properties_manager The object that manages the
 * physical properties of the problem.
 * @param[in] material_id The material id of the cell.
 * @param[in] field_id The field for which the dependency is checked.
 * @return true if at least one electromagnetic property of the material
 * depends on the field, false otherwise.
 */
bool
electromagnetic_properties_depend_on(
  const PhysicalPropertiesManager &physical_properties_manager,
  const unsigned int               material_id,
  const field                      field_id);


/**
 * @brief Class that stores the information required by the assembly of the
 * Discontinuous Petrov-Galerkin (DPG) ultraweak formulation of the
 * time-harmonic Maxwell equations. For each cell, it evaluates the complex
 * shape functions of the test space (\f$\mathbf{F}\f$ for the electric test
 * functions and \f$\mathbf{I}\f$ for the magnetic ones, with their curls and
 * complex conjugates) and of the interior trial space (\f$\mathbf{E}\f$ and
 * \f$\mathbf{H}\f$) and the effective material properties at the quadrature
 * points. For each face of the cell, it evaluates the test functions and the
 * tangential traces of the skeleton trial space (\f$\hat{\mathbf{E}}\f$ and
 * \f$\hat{\mathbf{H}}\f$) at the face quadrature points. The complex values are
 * stored with their complex conjugates because the conjugate of a complex
 * tensor is not implemented in deal.II.
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
  {
    gather_temperature   = false;
    gather_void_fraction = false;
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
  {
    gather_temperature   = sd.gather_temperature;
    gather_void_fraction = sd.gather_void_fraction;

    allocate();

    if (sd.gather_temperature)
      enable_temperature(sd.fe_values_temperature->get_fe(),
                         sd.fe_values_temperature->get_quadrature(),
                         sd.fe_face_values_temperature->get_quadrature(),
                         sd.fe_values_temperature->get_mapping());

    if (sd.gather_void_fraction)
      enable_void_fraction(sd.fe_values_void_fraction->get_fe(),
                           sd.fe_values_void_fraction->get_quadrature(),
                           sd.fe_face_values_void_fraction->get_quadrature(),
                           sd.fe_values_void_fraction->get_mapping());
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
   * @brief Enable the evaluation of the void fraction field, on which the
   * electromagnetic properties of an effective medium may depend.
   *
   * @param[in] fe_void_fraction Finite element of the void fraction.
   *
   * @param[in] cell_quadrature Quadrature used for the cell integrals.
   *
   * @param[in] face_quadrature Quadrature used for the face integrals.
   *
   * @param[in] mapping Mapping of the domain.
   */
  void
  enable_void_fraction(const FiniteElement<dim>  &fe_void_fraction,
                       const Quadrature<dim>     &cell_quadrature,
                       const Quadrature<dim - 1> &face_quadrature,
                       const Mapping<dim>        &mapping);

  /**
   * @brief Reinitialize the content of the scratch for a cell. It evaluates
   * the shape functions of the test and interior trial spaces at the
   * quadrature points. The temperature values are set to zero and are only
   * overwritten by reinit_temperature when the material of the cell depends on
   * the temperature. Likewise, the void fraction values are set to one (pure
   * fluid) and are only overwritten by reinit_void_fraction when the material
   * of the cell depends on the void fraction.
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
   * when the material of the cell depends on the temperature. Likewise, the
   * face void fraction values are set to one (pure fluid) and are only
   * overwritten by reinit_face_void_fraction when the material of the cell
   * depends on the void fraction.
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
   * @brief Evaluate the void fraction field at the cell quadrature points.
   *
   * @tparam VectorType The vector type of the void fraction solution.
   *
   * @param[in] cell_void_fraction The current cell, for the void fraction
   * DoFHandler.
   *
   * @param[in] void_fraction_solution Void fraction solution.
   */
  template <typename VectorType>
  void
  reinit_void_fraction(
    const typename DoFHandler<dim>::active_cell_iterator &cell_void_fraction,
    const VectorType &void_fraction_solution)
  {
    fe_values_void_fraction->reinit(cell_void_fraction);
    fe_values_void_fraction->get_function_values(void_fraction_solution,
                                                 this->void_fraction_values);
  }

  /**
   * @brief Evaluate the void fraction field at the quadrature points of a face
   * of the current cell.
   *
   * @tparam VectorType The vector type of the void fraction solution.
   *
   * @param[in] cell_void_fraction The current cell, for the void fraction
   * DoFHandler.
   *
   * @param[in] face_no Index of the face within the cell.
   *
   * @param[in] void_fraction_solution Void fraction solution.
   */
  template <typename VectorType>
  void
  reinit_face_void_fraction(
    const typename DoFHandler<dim>::active_cell_iterator &cell_void_fraction,
    const unsigned int                                    face_no,
    const VectorType &void_fraction_solution)
  {
    fe_face_values_void_fraction->reinit(cell_void_fraction, face_no);
    fe_face_values_void_fraction->get_function_values(
      void_fraction_solution, this->face_void_fraction_values);
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
  // ID of the material of the current cell, used to determine the appropriate
  // material properties from the input parameters.
  unsigned int material_id;
  // Flag to indicate whether the material of the cell depends on the
  // temperature.
  bool cell_material_needs_temperature;
  // Flag to indicate whether the material of the cell depends on the void
  // fraction.
  bool cell_material_needs_void_fraction;

  // FEValues of the interior trial space
  FEValues<dim> fe_values_trial_interior;
  // FEValues of the test space
  FEValues<dim> fe_values_test;
  // FEFaceValues of the skeleton trial space
  FEFaceValues<dim> fe_face_values_trial_skeleton;
  // FEFaceValues of the test space
  FEFaceValues<dim> fe_face_values_test;

  // Extractor for the FESystem of the real component of the electric field
  const FEValuesExtractors::Vector extractor_E_real{0};
  // Extractor for the FESystem of the imaginary component of the electric field
  const FEValuesExtractors::Vector extractor_E_imag{dim};
  // Extractor for the FESystem of the real component of the magnetic field
  const FEValuesExtractors::Vector extractor_H_real{2 * dim};
  // Extractor for the FESystem of the imaginary component of the magnetic field
  const FEValuesExtractors::Vector extractor_H_imag{3 * dim};

  // Number of quadrature points
  unsigned int n_q_points;
  // Number of face quadrature points
  unsigned int n_face_q_points;
  // Number of dofs of the test space
  unsigned int n_dofs_test;
  // Number of dofs of the interior trial space
  unsigned int n_dofs_trial_interior;
  // Number of dofs of the skeleton trial space
  unsigned int n_dofs_trial_skeleton;

  // Classification of the shape functions of the test space
  std::vector<unsigned char> shape_function_type_test;
  // Classification of the shape functions of the interior trial space
  std::vector<unsigned char> shape_function_type_trial_interior;
  // Classification of the shape functions of the skeleton trial space
  std::vector<unsigned char> shape_function_type_trial_skeleton;

  // Local dof indices of the test functions associated with the electric field
  // (F)
  std::vector<unsigned int> test_dofs_electric;
  // Local dof indices of the test functions associated with the magnetic field
  // (I)
  std::vector<unsigned int> test_dofs_magnetic;
  // Local dof indices of the interior trial functions associated with the
  // electric field (E)
  std::vector<unsigned int> trial_interior_dofs_electric;
  // Local dof indices of the interior trial functions associated with the
  // magnetic field (H)
  std::vector<unsigned int> trial_interior_dofs_magnetic;
  // Local dof indices of the skeleton trial functions associated with the
  // electric field (E_hat)
  std::vector<unsigned int> trial_skeleton_dofs_electric;
  // Local dof indices of the skeleton trial functions associated with the
  // magnetic field (H_hat)
  std::vector<unsigned int> trial_skeleton_dofs_magnetic;

  // Cell container for the volume integrals weights
  std::vector<double> JxW;
  // Cell container for the quadrature points
  std::vector<Point<dim>> quadrature_points;

  // Temperature values at the cell quadrature points
  std::vector<double> temperature_values;
  // Void fraction values at the cell quadrature points
  std::vector<double> void_fraction_values;
  // Map of the field values at the cell quadrature points, indexed by field
  std::map<field, std::vector<double>> fields;
  // Effective electric permittivities at the cell quadrature points
  std::vector<std::complex<double>> effective_electric_permittivities;
  // Effective magnetic permeabilities at the cell quadrature points
  std::vector<std::complex<double>> effective_magnetic_permeabilities;

  // Cell electric test functions
  Table<2, Tensor<1, dim, std::complex<double>>> phi_F;
  // Cell complex conjugate of the electric test functions
  Table<2, Tensor<1, dim, std::complex<double>>> phi_F_conj;
  // Cell curl of the electric test functions
  Table<2, Tensor<1, dim, std::complex<double>>> curl_phi_F;
  // Cell complex conjugate of the curl of the electric test functions
  Table<2, Tensor<1, dim, std::complex<double>>> curl_phi_F_conj;
  // Cell magnetic test functions
  Table<2, Tensor<1, dim, std::complex<double>>> phi_I;
  // Cell complex conjugate of the magnetic test functions
  Table<2, Tensor<1, dim, std::complex<double>>> phi_I_conj;
  // Cell curl of the magnetic test functions
  Table<2, Tensor<1, dim, std::complex<double>>> curl_phi_I;
  // Cell complex conjugate of the curl of the magnetic test functions
  Table<2, Tensor<1, dim, std::complex<double>>> curl_phi_I_conj;

  // Cell interior trial electric functions
  Table<2, Tensor<1, dim, std::complex<double>>> phi_E;
  // Cell interior trial magnetic functions
  Table<2, Tensor<1, dim, std::complex<double>>> phi_H;

  // Flag to indicate whether the current face is at the boundary of the domain
  bool face_at_boundary;
  // Boundary ID of the current face, used to determine the appropriate boundary
  // conditions from the input parameters
  types::boundary_id face_boundary_id;

  // Face container for the face integrals weights
  std::vector<double> face_JxW;
  // Face container for the face quadrature points
  std::vector<Point<dim>> face_quadrature_points;
  // Face container for the face normal vectors
  std::vector<Tensor<1, dim>> face_normals;

  // Face temperature values at the face quadrature points
  std::vector<double> face_temperature_values;
  // Face void fraction values at the face quadrature points
  std::vector<double> face_void_fraction_values;
  // Map of the field values at the face quadrature points, indexed by field
  std::map<field, std::vector<double>> face_fields;
  // Effective electric permittivities at the face quadrature points
  std::vector<std::complex<double>> face_effective_electric_permittivities;
  // Effective magnetic permeabilities at the face quadrature points
  std::vector<std::complex<double>> face_effective_magnetic_permeabilities;

  // Face electric test functions at the face quadrature points
  Table<2, Tensor<1, dim, std::complex<double>>> phi_F_face;
  // Face complex conjugate of the electric test functions at the face
  // quadrature points
  Table<2, Tensor<1, dim, std::complex<double>>> phi_F_face_conj;
  // Face magnetic test functions at the face quadrature points
  Table<2, Tensor<1, dim, std::complex<double>>> phi_I_face_conj;
  // Face normal cross product of the magnetic test functions at the face
  // quadrature points
  Table<2, Tensor<1, dim, std::complex<double>>> n_cross_phi_I_face;
  // Face complex conjugate of the normal cross product of the magnetic test
  // functions at the face quadrature points
  Table<2, Tensor<1, dim, std::complex<double>>> n_cross_phi_I_face_conj;

  // Face skeleton trial electric functions at the face quadrature points
  Table<2, Tensor<1, dim, std::complex<double>>> phi_E_hat;
  // Face normal cross product of the skeleton trial electric functions at the
  // face quadrature points
  Table<2, Tensor<1, dim, std::complex<double>>> n_cross_phi_E_hat;
  // Face normal cross product of the complex conjugate of the skeleton trial
  // magnetic functions at the face quadrature points
  Table<2, Tensor<1, dim, std::complex<double>>> n_cross_phi_H_hat;

  // Flag to indicate whether the temperature field needs to be gathered and
  // evaluated at the quadrature points
  bool gather_temperature;
  // FEValues of the temperature field
  std::shared_ptr<FEValues<dim>> fe_values_temperature;
  // FEFaceValues of the temperature field
  std::shared_ptr<FEFaceValues<dim>> fe_face_values_temperature;

  // Flag to indicate whether the void fraction field needs to be gathered and
  // evaluated at the quadrature points
  bool gather_void_fraction;
  // FEValues of the void fraction field
  std::shared_ptr<FEValues<dim>> fe_values_void_fraction;
  // FEFaceValues of the void fraction field
  std::shared_ptr<FEFaceValues<dim>> fe_face_values_void_fraction;
};

#endif
