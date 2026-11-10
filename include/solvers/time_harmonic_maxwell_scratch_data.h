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

#include <array>
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
 * time-harmonic Maxwell equations. For each cell, it evaluates the base
 * functions of the test space (\f$\mathbf{F}\f$ and \f$\mathbf{I}\f$) with
 * their curls, the base functions of the interior trial space (\f$\mathbf{E}\f$
 * and \f$\mathbf{H}\f$) and the effective material properties at the quadrature
 * points. Additionnaly, for each face of the cell, it evaluates the base
 * functions of the tangential traces of the test space and the base functions
 * of the skeleton trial space (\f$\hat{\mathbf{E}}\f$ and
 * \f$\hat{\mathbf{H}}\f$), with their cross products with the normal. Since the
 * Nedelec base functions whose degrees of freedom are not located on a face
 * have a zero tangential trace on this face, only the base functions associated
 * with the face are evaluated.
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
   * @brief Fields of the finite element systems, numbered in the order in which
   * they appear in the systems following the convention : [E_real, E_imag,
   * H_real, H_imag].
   */
  enum ElectromagneticFieldComponents : unsigned int
  {
    E_real = 0,
    E_imag = 1,
    H_real = 2,
    H_imag = 3
  };

  // Number of fields of the finite element systems
  static constexpr unsigned int n_electromagnetic_fields = 4;

  /**
   * @brief Constructor. Creates the FEValues and FEFaceValues objects of the
   * interior trial, skeleton trial and test spaces and allocates the memory of
   * the scratch. It does not do any evaluation, since this is done at the cell
   * level.
   *
   * @param[in] properties_manager Manager of the physical property models.
   *
   * @param[in] omega Angular frequency of the electromagnetic fields.
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
    const double                     omega,
    const FiniteElement<dim>        &fe_trial_interior,
    const FiniteElement<dim>        &fe_trial_skeleton,
    const FiniteElement<dim>        &fe_test,
    const Quadrature<dim>           &cell_quadrature,
    const Quadrature<dim - 1>       &face_quadrature,
    const Mapping<dim>              &mapping)
    : properties_manager(properties_manager)
    , omega(omega)
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
    , omega(sd.omega)
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
   * @brief Allocate the memory of all the members of the scratch. It builds the lists
   * of the local dof index of each base function for each of the four
   * fields of the three finite element spaces for the reference cell. It also
   * builds the lists of the base functions of the test and skeleton trial
   * spaces whose tangential trace does not vanish on each face of the
   * reference cell.
   */
  void
  allocate() override;

  /**
   * @brief Build the local dof index list of each base function for each field of
   * a finite element space of the time-harmonic Maxwell equations [E_real,
   * E_imag, H_real, H_imag].
   *
   * @param[in] fe Finite element system of the [E_real, E_imag, H_real,
   * H_imag] fields (interior trial, skeleton trial or test space).
   *
   * @param[out] local_dof_index_map Local dof index of each base function for
   * each field: local_dof_index_map[field][i].
   */
  void
  build_local_dof_index_map(
    const FiniteElement<dim> &fe,
    std::array<std::vector<unsigned int>, n_electromagnetic_fields>
      &local_dof_index_map) const;

  /**
   * @brief Identify, for each face of the reference cell, the base functions of
   * a finite element space whose tangential trace does not vanish on the face
   * (the relevant ones for the face integrals, this would not be needed if the
   * FE_FaceNedelec existed). For Nedelec elements, these are the base functions
   * whose degrees of freedom are located on the edges or the interior of the
   * face.
   *
   * @param[in] fe Finite element system of the [E_real, E_imag, H_real,
   * H_imag] fields.
   *
   * @param[in] local_dof_index_map Local dof index of each base function for
   * each field (see compute_local_dof_index_map()).
   *
   * @param[out] tangential_base_functions Base functions whose tangential
   * trace does not vanish on each face: tangential_base_functions[face_no].
   */
  void
  identify_tangential_base_functions(
    const FiniteElement<dim> &fe,
    const std::array<std::vector<unsigned int>, n_electromagnetic_fields>
                                           &local_dof_index_map,
    std::vector<std::vector<unsigned int>> &tangential_base_functions) const;

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
   * @brief Compute the coefficients of the time-harmonic Maxwell equations of
   * the cell at the quadrature points.
   */
  void
  calculate_cell_coefficients();

  /**
   * @brief Compute the effective electromagnetic properties at the face
   * quadrature points from the fields gathered by the scratch.
   */
  void
  calculate_face_physical_properties();

  // Physical properties
  const PhysicalPropertiesManager properties_manager;
  // Angular frequency of the electromagnetic fields
  const double omega;
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
  const FEValuesExtractors::Vector extractor_E_real{E_real * dim};
  // Extractor for the FESystem of the imaginary component of the electric field
  const FEValuesExtractors::Vector extractor_E_imag{E_imag * dim};
  // Extractor for the FESystem of the real component of the magnetic field
  const FEValuesExtractors::Vector extractor_H_real{H_real * dim};
  // Extractor for the FESystem of the imaginary component of the magnetic field
  const FEValuesExtractors::Vector extractor_H_imag{H_imag * dim};

  // Number of quadrature points
  unsigned int n_q_points;
  // Number of face quadrature points
  unsigned int n_face_q_points;
  // Number of base functions of the test space
  unsigned int n_base_functions_test;
  // Number of base functions of the interior trial space
  unsigned int n_base_functions_trial_interior;
  // Number of base functions of the skeleton trial space
  unsigned int n_base_functions_trial_skeleton;

  // Local dof index of each base function of the test space for each field:
  // test_local_dof_index_map[field][i]. The electric test functions (F) are
  // the E_real and E_imag fields of the test space and the magnetic ones (I)
  // are the H_real and H_imag fields.
  std::array<std::vector<unsigned int>, n_electromagnetic_fields>
    test_local_dof_index_map;
  // Local dof index of each base function of the interior trial space for each
  // field: trial_interior_local_dof_index_map[field][i]
  std::array<std::vector<unsigned int>, n_electromagnetic_fields>
    trial_interior_local_dof_index_map;
  // Local dof index of each base function of the skeleton trial space for each
  // field: trial_skeleton_local_dof_index_map[field][i]
  std::array<std::vector<unsigned int>, n_electromagnetic_fields>
    trial_skeleton_local_dof_index_map;

  // Base functions of the test space whose tangential trace does not vanish on
  // each face of the reference cell: test_tangential_base_functions[face_no]
  std::vector<std::vector<unsigned int>> test_tangential_base_functions;
  // Base functions of the skeleton trial space whose tangential trace does not
  // vanish on each face of the reference cell:
  // trial_skeleton_tangential_base_functions[face_no]
  std::vector<std::vector<unsigned int>>
    trial_skeleton_tangential_base_functions;

  // Cell container for the volume integrals weights
  std::vector<double> JxW;
  // (1 + |iweps|^2) JxW, weight of (phi_i, phi_j) in the FF block of G
  std::vector<double> JxW_FF_container;
  // (1 + |iwmu|^2) JxW, weight of (phi_i, phi_j) in the II block of G
  std::vector<double> JxW_II_container;
  // iweps JxW
  std::vector<std::complex<double>> iweps_JxW_container;
  // iwmu JxW
  std::vector<std::complex<double>> iwmu_JxW_container;
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

  // Cell base functions of the test space
  Table<2, Tensor<1, dim>> phi_test;
  // Cell curl of the base functions of the test space
  Table<2, Tensor<1, dim>> curl_phi_test;
  // Cell base functions of the interior trial space
  Table<2, Tensor<1, dim>> phi_trial_interior;

  // Index of the current face within the cell
  unsigned int face_no;
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

  // Face base functions of the test space: phi_test_face[i][q]
  Table<2, Tensor<1, dim>> phi_test_face;
  // Face cross product of the normal with the base functions of the test space
  Table<2, Tensor<1, dim>> n_cross_phi_test_face;
  // Face tangential traces of the base functions of the skeleton trial space
  Table<2, Tensor<1, dim>> tangential_phi_trial_skeleton;
  // Face cross product of the normal with the tangential traces of the base
  // functions of the skeleton trial space
  Table<2, Tensor<1, dim>> n_cross_tangential_phi_trial_skeleton;

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
