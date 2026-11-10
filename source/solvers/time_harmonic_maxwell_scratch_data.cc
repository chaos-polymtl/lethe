// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/tensors_and_points_dimension_manipulation.h>

#include <solvers/time_harmonic_maxwell.h>
#include <solvers/time_harmonic_maxwell_scratch_data.h>

#include <deal.II/base/exceptions.h>

#include <algorithm>
#include <string>

void
compute_effective_electromagnetic_properties(
  const PhysicalPropertiesManager            &physical_properties_manager,
  const std::map<field, std::vector<double>> &field_values_vectors,
  const unsigned int                          material_id,
  std::vector<std::complex<double>>          &effective_electric_permittivities,
  std::vector<std::complex<double>>          &effective_magnetic_permeabilities)
{
  const unsigned int n_q_points =
    field_values_vectors.at(field::temperature).size();

  // Create temporary vectors to store the real and imaginary parts of the
  // electromagnetic properties
  std::vector<double> permittivity_real(n_q_points);
  std::vector<double> permittivity_imag(n_q_points);
  std::vector<double> conductivity(n_q_points);
  std::vector<double> permeability_real(n_q_points);
  std::vector<double> permeability_imag(n_q_points);

  // We use the fluid_id = 0 because the time-harmonic Maxwell solver only
  // supports a single fluid. The material_id is used to identify the different
  // materials in the input file and assign them their corresponding properties.
  if (physical_properties_manager.is_effective_medium(0, material_id))
    {
      // The material is a mixture of the fluid and of a dispersed phase. The
      // time-harmonic Maxwell solver supports a single fluid (id 0) and a
      // single dispersed phase (id 0).
      const unsigned int material_interaction_id =
        physical_properties_manager.get_material_interaction_id(
          material_interactions_type::fluid_dispersed_phase, 0, 0);

      physical_properties_manager
        .get_effective_electric_permittivity_real(material_interaction_id)
        ->vector_value(field_values_vectors, permittivity_real);
      physical_properties_manager
        .get_effective_electric_permittivity_imag(material_interaction_id)
        ->vector_value(field_values_vectors, permittivity_imag);
      physical_properties_manager
        .get_effective_magnetic_permeability_real(material_interaction_id)
        ->vector_value(field_values_vectors, permeability_real);
      physical_properties_manager
        .get_effective_magnetic_permeability_imag(material_interaction_id)
        ->vector_value(field_values_vectors, permeability_imag);
    }
  else
    {
      // The mixing law is applied to the complete complex permittivity
      // eps' + i(eps'' + sigma) of each phase, so the conductivity of both
      // phases is already included in the imaginary part of the effective
      // permittivity. The conductivity vector is therefore left at zero.
      physical_properties_manager
        .get_electric_permittivity_real(0, material_id)
        ->vector_value(field_values_vectors, permittivity_real);
      physical_properties_manager
        .get_electric_permittivity_imag(0, material_id)
        ->vector_value(field_values_vectors, permittivity_imag);
      physical_properties_manager.get_electric_conductivity(0, material_id)
        ->vector_value(field_values_vectors, conductivity);
      physical_properties_manager
        .get_magnetic_permeability_real(0, material_id)
        ->vector_value(field_values_vectors, permeability_real);
      physical_properties_manager
        .get_magnetic_permeability_imag(0, material_id)
        ->vector_value(field_values_vectors, permeability_imag);
    }

  // We resize the effective electromagnetic property vectors to match the
  // number of points because this function can be called with a different
  // number of points (e.g. cell and face quadratures)
  effective_electric_permittivities.resize(n_q_points);
  effective_magnetic_permeabilities.resize(n_q_points);

  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      effective_electric_permittivities[q] = {permittivity_real[q],
                                              permittivity_imag[q] +
                                                conductivity[q]};
      effective_magnetic_permeabilities[q] = {permeability_real[q],
                                              permeability_imag[q]};
    }
}

bool
electromagnetic_properties_depend_on(
  const PhysicalPropertiesManager &physical_properties_manager,
  const unsigned int               material_id,
  const field                      field_id)
{
  if (physical_properties_manager.is_effective_medium(0, material_id))
    {
      // The time-harmonic Maxwell solver supports a single fluid (id 0) and a
      // single dispersed phase (id 0). The dependencies of the effective
      // permittivity models include those of the conductivity models of the
      // phases used by the mixing law.
      const unsigned int material_interaction_id =
        physical_properties_manager.get_material_interaction_id(
          material_interactions_type::fluid_dispersed_phase, 0, 0);

      return physical_properties_manager
               .get_effective_electric_permittivity_real(
                 material_interaction_id)
               ->depends_on(field_id) ||
             physical_properties_manager
               .get_effective_electric_permittivity_imag(
                 material_interaction_id)
               ->depends_on(field_id) ||
             physical_properties_manager
               .get_effective_magnetic_permeability_real(
                 material_interaction_id)
               ->depends_on(field_id) ||
             physical_properties_manager
               .get_effective_magnetic_permeability_imag(
                 material_interaction_id)
               ->depends_on(field_id);
    }

  return physical_properties_manager.get_electric_conductivity(0, material_id)
           ->depends_on(field_id) ||
         physical_properties_manager
           .get_electric_permittivity_real(0, material_id)
           ->depends_on(field_id) ||
         physical_properties_manager
           .get_electric_permittivity_imag(0, material_id)
           ->depends_on(field_id) ||
         physical_properties_manager
           .get_magnetic_permeability_real(0, material_id)
           ->depends_on(field_id) ||
         physical_properties_manager
           .get_magnetic_permeability_imag(0, material_id)
           ->depends_on(field_id);
}

template <int dim>
void
TimeHarmonicMaxwellScratchData<dim>::allocate()
{
  const FiniteElement<dim> &fe_test = this->fe_values_test.get_fe();
  const FiniteElement<dim> &fe_trial_interior =
    this->fe_values_trial_interior.get_fe();
  const FiniteElement<dim> &fe_trial_skeleton =
    this->fe_face_values_trial_skeleton.get_fe();

  // Initialize the sizes
  this->n_q_points      = this->fe_values_test.get_quadrature().size();
  this->n_face_q_points = this->fe_face_values_test.get_quadrature().size();

  // Compute the local dof index of each base function for each field of each
  // space. This is only done once since it only depends on the finite
  // elements.
  std::vector<unsigned int> test_base_function_of_dof;
  std::vector<unsigned int> trial_interior_base_function_of_dof;
  std::vector<unsigned int> trial_skeleton_base_function_of_dof;
  compute_base_dofs(fe_test, this->test_base_dofs, test_base_function_of_dof);
  compute_base_dofs(fe_trial_interior,
                    this->trial_interior_base_dofs,
                    trial_interior_base_function_of_dof);
  compute_base_dofs(fe_trial_skeleton,
                    this->trial_skeleton_base_dofs,
                    trial_skeleton_base_function_of_dof);

  this->n_base_functions_test = this->test_base_dofs[E_real].size();
  this->n_base_functions_trial_interior =
    this->trial_interior_base_dofs[E_real].size();
  this->n_base_functions_trial_skeleton =
    this->trial_skeleton_base_dofs[E_real].size();

  // Base functions of the test and skeleton trial spaces associated with each
  // face of the reference cell
  compute_face_base_functions(fe_test,
                              test_base_function_of_dof,
                              this->face_test_base_functions);
  compute_face_base_functions(fe_trial_skeleton,
                              trial_skeleton_base_function_of_dof,
                              this->face_trial_skeleton_base_functions);

  // Cell quadrature
  this->JxW               = std::vector<double>(n_q_points);
  this->quadrature_points = std::vector<Point<dim>>(n_q_points);

  // Temperature, void fraction and effective properties at the cell quadrature
  // points. The temperature and the void fraction are the only fields on which
  // the electromagnetic properties can depend.
  this->temperature_values   = std::vector<double>(n_q_points, 0.);
  this->void_fraction_values = std::vector<double>(n_q_points, 1.);
  this->fields.clear();
  this->fields.insert(
    std::pair<field, std::vector<double>>(field::temperature, n_q_points));
  this->fields.insert(
    std::pair<field, std::vector<double>>(field::void_fraction, n_q_points));
  this->effective_electric_permittivities =
    std::vector<std::complex<double>>(n_q_points);
  this->effective_magnetic_permeabilities =
    std::vector<std::complex<double>>(n_q_points);

  // Coefficients of the cell terms of the DPG system
  this->JxW_FF        = std::vector<double>(n_q_points);
  this->JxW_II        = std::vector<double>(n_q_points);
  this->iweps_JxW     = std::vector<std::complex<double>>(n_q_points);
  this->conj_iwmu_JxW = std::vector<std::complex<double>>(n_q_points);
  this->iwmu_JxW      = std::vector<std::complex<double>>(n_q_points);

  // Test and interior trial functions at the cell quadrature points
  this->phi_test.reinit(n_base_functions_test, n_q_points);
  this->curl_phi_test.reinit(n_base_functions_test, n_q_points);
  this->phi_trial_interior.reinit(n_base_functions_trial_interior, n_q_points);

  // Face quadrature
  this->face_no                = 0;
  this->face_at_boundary       = false;
  this->face_boundary_id       = numbers::invalid_boundary_id;
  this->face_JxW               = std::vector<double>(n_face_q_points);
  this->face_quadrature_points = std::vector<Point<dim>>(n_face_q_points);
  this->face_normals           = std::vector<Tensor<1, dim>>(n_face_q_points);

  // Temperature, void fraction and effective properties at the face quadrature
  // points
  this->face_temperature_values   = std::vector<double>(n_face_q_points, 0.);
  this->face_void_fraction_values = std::vector<double>(n_face_q_points, 1.);
  this->face_fields.clear();
  this->face_fields.insert(
    std::pair<field, std::vector<double>>(field::temperature, n_face_q_points));
  this->face_fields.insert(
    std::pair<field, std::vector<double>>(field::void_fraction,
                                          n_face_q_points));
  this->face_effective_electric_permittivities =
    std::vector<std::complex<double>>(n_face_q_points);
  this->face_effective_magnetic_permeabilities =
    std::vector<std::complex<double>>(n_face_q_points);

  // Test and skeleton trial base functions at the face quadrature points
  this->phi_test_face.reinit(n_base_functions_test, n_face_q_points);
  this->n_cross_phi_test_face.reinit(n_base_functions_test, n_face_q_points);
  this->tangential_phi_trial_skeleton.reinit(n_base_functions_trial_skeleton,
                                             n_face_q_points);
  this->n_cross_tangential_phi_trial_skeleton.reinit(
    n_base_functions_trial_skeleton, n_face_q_points);
}

template <int dim>
void
TimeHarmonicMaxwellScratchData<dim>::compute_base_dofs(
  const FiniteElement<dim>                                        &fe,
  std::array<std::vector<unsigned int>, n_electromagnetic_fields> &base_dofs,
  std::vector<unsigned int> &base_function_of_dof) const
{
  // Since the four fields use the same base element, each field has a quarter
  // of the shape functions of the finite element system.
  const unsigned int n_base_functions =
    fe.n_dofs_per_cell() / n_electromagnetic_fields;

  // Extractor of each field, indexed by the ElectromagneticFieldComponents
  // enumeration
  const std::array<FEValuesExtractors::Vector, n_electromagnetic_fields>
    field_extractors = {
      {extractor_E_real, extractor_E_imag, extractor_H_real, extractor_H_imag}};

  base_function_of_dof.assign(fe.n_dofs_per_cell(),
                              numbers::invalid_unsigned_int);
  for (auto &field_base_dofs : base_dofs)
    field_base_dofs.assign(n_base_functions, numbers::invalid_unsigned_int);

  for (unsigned int k = 0; k < fe.n_dofs_per_cell(); ++k)
    {
      // Field of the shape function. Each shape function must belong to
      // exactly one field.
      unsigned int field_index    = numbers::invalid_unsigned_int;
      unsigned int n_owner_fields = 0;
      for (unsigned int candidate_field = 0;
           candidate_field < n_electromagnetic_fields;
           ++candidate_field)
        if (fe.shape_function_belongs_to(k, field_extractors[candidate_field]))
          {
            field_index = candidate_field;
            ++n_owner_fields;
          }
      AssertThrow(n_owner_fields == 1,
                  ExcMessage("The shape function " + std::to_string(k) +
                             " of the finite element " + fe.get_name() +
                             " does not belong to exactly one of the E_real, "
                             "E_imag, H_real and H_imag fields."));

      // Base function of the shape function k. The function
      // system_to_base_index identifies k by three indices:
      //  - the base element it belongs to (one per field);
      //  - the component of the field it represents. FE_DGQ^dim builds a
      //    vector field from the scalar FE_DGQ, which represents one component
      //    (x, y or z) of the field. FE_Nedelec is a vector element that
      //    represents the three components at once, so this index is always 0;
      //  - its index within the base element.
      // The base function combines the field component and the index within
      // the base element to number the shape functions of one field from 0 to
      // n_base_functions - 1. Since the four base elements are identical, a
      // base function is the same function in the four fields.
      const auto        &base_index            = fe.system_to_base_index(k);
      const unsigned int base_element          = base_index.first.first;
      const unsigned int field_component       = base_index.first.second;
      const unsigned int index_in_base_element = base_index.second;
      const unsigned int n_dofs_base_element =
        fe.base_element(base_element).n_dofs_per_cell();
      const unsigned int base_function =
        field_component * n_dofs_base_element + index_in_base_element;

      // Each base function must appear exactly once in each field. Since the
      // finite element has 4 * n_base_functions shape functions, checking that
      // no entry is written twice also guarantees that all the entries of
      // base_dofs are filled.
      AssertThrow(base_function < n_base_functions &&
                    base_dofs[field_index][base_function] ==
                      numbers::invalid_unsigned_int,
                  ExcMessage(
                    "The four fields of the finite element " + fe.get_name() +
                    " do not use the same base element. The time-harmonic "
                    "Maxwell solver requires the E_real, E_imag, H_real and "
                    "H_imag fields to be discretized with the same finite "
                    "element."));

      base_dofs[field_index][base_function] = k;
      base_function_of_dof[k]               = base_function;
    }
}

template <int dim>
void
TimeHarmonicMaxwellScratchData<dim>::compute_face_base_functions(
  const FiniteElement<dim>               &fe,
  const std::vector<unsigned int>        &base_function_of_dof,
  std::vector<std::vector<unsigned int>> &face_base_functions) const
{
  const unsigned int n_faces = fe.reference_cell().n_faces();
  face_base_functions.assign(n_faces, std::vector<unsigned int>());

  for (unsigned int face = 0; face < n_faces; ++face)
    {
      // The degrees of freedom located on the face belong to the four fields,
      // so each base function associated with the face appears four times. The
      // duplicates are removed after sorting the list.
      face_base_functions[face].reserve(fe.n_dofs_per_face(face));
      for (unsigned int i = 0; i < fe.n_dofs_per_face(face); ++i)
        face_base_functions[face].emplace_back(
          base_function_of_dof[fe.face_to_cell_index(i, face)]);

      std::ranges::sort(face_base_functions[face]);
      const auto duplicates = std::ranges::unique(face_base_functions[face]);
      face_base_functions[face].erase(duplicates.begin(), duplicates.end());
    }
}

template <int dim>
void
TimeHarmonicMaxwellScratchData<dim>::enable_temperature(
  const FiniteElement<dim>  &fe_temperature,
  const Quadrature<dim>     &cell_quadrature,
  const Quadrature<dim - 1> &face_quadrature,
  const Mapping<dim>        &mapping)
{
  this->gather_temperature = true;
  this->fe_values_temperature =
    std::make_shared<FEValues<dim>>(mapping,
                                    fe_temperature,
                                    cell_quadrature,
                                    update_values | update_quadrature_points);
  this->fe_face_values_temperature =
    std::make_shared<FEFaceValues<dim>>(mapping,
                                        fe_temperature,
                                        face_quadrature,
                                        update_values |
                                          update_quadrature_points);
}

template <int dim>
void
TimeHarmonicMaxwellScratchData<dim>::enable_void_fraction(
  const FiniteElement<dim>  &fe_void_fraction,
  const Quadrature<dim>     &cell_quadrature,
  const Quadrature<dim - 1> &face_quadrature,
  const Mapping<dim>        &mapping)
{
  this->gather_void_fraction = true;
  this->fe_values_void_fraction =
    std::make_shared<FEValues<dim>>(mapping,
                                    fe_void_fraction,
                                    cell_quadrature,
                                    update_values | update_quadrature_points);
  this->fe_face_values_void_fraction =
    std::make_shared<FEFaceValues<dim>>(mapping,
                                        fe_void_fraction,
                                        face_quadrature,
                                        update_values |
                                          update_quadrature_points);
}

template <int dim>
void
TimeHarmonicMaxwellScratchData<dim>::reinit(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const typename DoFHandler<dim>::active_cell_iterator &cell_test)
{
  // The curl of the test functions and the cross products are only defined in
  // 3D. The dimension is checked at compile time so the 3D operations are not
  // compiled when dim = 2.
  if constexpr (dim != 3)
    {
      (void)cell;
      (void)cell_test;
      AssertThrow(false, TimeHarmonicMaxwellDimensionNotSupported(dim));
    }
  else
    {
      this->material_id = cell->material_id();

      // We check if the physical properties depend on the temperature and
      // void fraction fields. If so, the solver will evaluate them at the
      // quadrature points using reinit_temperature and reinit_void_fraction.
      this->cell_material_needs_temperature =
        electromagnetic_properties_depend_on(properties_manager,
                                             material_id,
                                             field::temperature);
      this->cell_material_needs_void_fraction =
        electromagnetic_properties_depend_on(properties_manager,
                                             material_id,
                                             field::void_fraction);

      this->fe_values_trial_interior.reinit(cell);
      this->fe_values_test.reinit(cell_test);

      this->quadrature_points =
        this->fe_values_trial_interior.get_quadrature_points();

      for (unsigned int q = 0; q < n_q_points; ++q)
        this->JxW[q] = this->fe_values_trial_interior.JxW(q);

      // Base functions of the test space and their curl. Since the four
      // fields share the same base functions, they are evaluated with the
      // E_real shape functions only.
      for (unsigned int i = 0; i < n_base_functions_test; ++i)
        {
          const unsigned int k = this->test_base_dofs[E_real][i];
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              phi_test[i][q] = fe_values_test[extractor_E_real].value(k, q);
              curl_phi_test[i][q] = fe_values_test[extractor_E_real].curl(k, q);
            }
        }

      // Base functions of the interior trial space, evaluated in the same
      // way
      for (unsigned int j = 0; j < n_base_functions_trial_interior; ++j)
        {
          const unsigned int k = this->trial_interior_base_dofs[E_real][j];
          for (unsigned int q = 0; q < n_q_points; ++q)
            phi_trial_interior[j][q] =
              fe_values_trial_interior[extractor_E_real].value(k, q);
        }
    }
}

template <int dim>
void
TimeHarmonicMaxwellScratchData<dim>::reinit_face(
  const typename DoFHandler<dim>::active_cell_iterator &cell_skeleton,
  const typename DoFHandler<dim>::active_cell_iterator &cell_test,
  const unsigned int                                    face_no)
{
  if constexpr (dim != 3)
    {
      (void)cell_skeleton;
      (void)cell_test;
      (void)face_no;
      AssertThrow(false, TimeHarmonicMaxwellDimensionNotSupported(dim));
    }
  else
    {
      const auto face        = cell_skeleton->face(face_no);
      this->face_no          = face_no;
      this->face_at_boundary = face->at_boundary();
      this->face_boundary_id = this->face_at_boundary ?
                                 face->boundary_id() :
                                 numbers::invalid_boundary_id;
      std::ranges::fill(this->face_temperature_values, 0.);
      std::ranges::fill(this->face_void_fraction_values, 1.);

      this->fe_face_values_test.reinit(cell_test, face_no);
      this->fe_face_values_trial_skeleton.reinit(cell_skeleton, face_no);

      for (unsigned int q = 0; q < n_face_q_points; ++q)
        {
          this->face_JxW[q] = this->fe_face_values_trial_skeleton.JxW(q);
          this->face_quadrature_points[q] =
            this->fe_face_values_trial_skeleton.quadrature_point(q);
          this->face_normals[q] =
            this->fe_face_values_trial_skeleton.normal_vector(q);
        }

      // Base functions of the test space associated with the face, evaluated
      // with the E_real shape functions since the four fields share the same
      // base functions. The face terms of the bilinear forms, of the test norm
      // and of the load vector only involve the tangential traces of the test
      // functions, which vanish on this face for the other base functions.
      for (const unsigned int i : this->face_test_base_functions[face_no])
        {
          const unsigned int k = this->test_base_dofs[E_real][i];
          for (unsigned int q = 0; q < n_face_q_points; ++q)
            {
              phi_test_face[i][q] =
                fe_face_values_test[extractor_E_real].value(k, q);
              n_cross_phi_test_face[i][q] =
                cross_product_3d(face_normals[q], phi_test_face[i][q]);
            }
        }

      // Base functions of the skeleton trial space associated with the face.
      // To be in H^-1/2(curl), the fields need to have the tangential trace
      // mapping (n x (E x n)) which extracts the tangential component of the
      // field at the face.
      for (const unsigned int j :
           this->face_trial_skeleton_base_functions[face_no])
        {
          const unsigned int k = this->trial_skeleton_base_dofs[E_real][j];
          for (unsigned int q = 0; q < n_face_q_points; ++q)
            {
              tangential_phi_trial_skeleton[j][q] = tangential_component(
                fe_face_values_trial_skeleton[extractor_E_real].value(k, q),
                face_normals[q]);
              n_cross_tangential_phi_trial_skeleton[j][q] =
                cross_product_3d(face_normals[q],
                                 tangential_phi_trial_skeleton[j][q]);
            }
        }
    }
}

template <int dim>
void
TimeHarmonicMaxwellScratchData<dim>::calculate_physical_properties()
{
  set_field_vector(field::temperature, this->temperature_values, this->fields);
  set_field_vector(field::void_fraction,
                   this->void_fraction_values,
                   this->fields);
  compute_effective_electromagnetic_properties(
    this->properties_manager,
    this->fields,
    this->material_id,
    this->effective_electric_permittivities,
    this->effective_magnetic_permeabilities);
}

template <int dim>
void
TimeHarmonicMaxwellScratchData<dim>::calculate_cell_coefficients()
{
  static constexpr std::complex<double> imag{0., 1.};

  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      const std::complex<double> iweps =
        imag * omega * this->effective_electric_permittivities[q];
      const std::complex<double> iwmu =
        imag * omega * this->effective_magnetic_permeabilities[q];

      this->JxW_FF[q]        = (1. + std::norm(iweps)) * this->JxW[q];
      this->JxW_II[q]        = (1. + std::norm(iwmu)) * this->JxW[q];
      this->iweps_JxW[q]     = iweps * this->JxW[q];
      this->conj_iwmu_JxW[q] = std::conj(iwmu) * this->JxW[q];
      this->iwmu_JxW[q]      = iwmu * this->JxW[q];
    }
}

template <int dim>
void
TimeHarmonicMaxwellScratchData<dim>::calculate_face_physical_properties()
{
  set_field_vector(field::temperature,
                   this->face_temperature_values,
                   this->face_fields);
  set_field_vector(field::void_fraction,
                   this->face_void_fraction_values,
                   this->face_fields);
  compute_effective_electromagnetic_properties(
    this->properties_manager,
    this->face_fields,
    this->material_id,
    this->face_effective_electric_permittivities,
    this->face_effective_magnetic_permeabilities);
}

template class TimeHarmonicMaxwellScratchData<2>;
template class TimeHarmonicMaxwellScratchData<3>;
