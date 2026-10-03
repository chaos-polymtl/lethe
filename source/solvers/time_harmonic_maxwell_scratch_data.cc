// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/time_harmonic_maxwell_scratch_data.h>

#include <algorithm>
#include <ranges>

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

  physical_properties_manager.get_electric_permittivity_real(0, material_id)
    ->vector_value(field_values_vectors, permittivity_real);
  physical_properties_manager.get_electric_permittivity_imag(0, material_id)
    ->vector_value(field_values_vectors, permittivity_imag);
  physical_properties_manager.get_electric_conductivity(0, material_id)
    ->vector_value(field_values_vectors, conductivity);
  physical_properties_manager.get_magnetic_permeability_real(0, material_id)
    ->vector_value(field_values_vectors, permeability_real);
  physical_properties_manager.get_magnetic_permeability_imag(0, material_id)
    ->vector_value(field_values_vectors, permeability_imag);

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
  this->n_dofs_test     = fe_test.n_dofs_per_cell();
  this->n_dofs_trial_interior = fe_trial_interior.n_dofs_per_cell();
  this->n_dofs_trial_skeleton = fe_trial_skeleton.n_dofs_per_cell();

  // Classify the shape functions of each space according to the field they
  // belong to and store the local dof index of each base function for each
  // field. This is only done once since it only depends on the finite
  // elements.
  std::vector<unsigned int> test_base_function_of_dof;
  std::vector<unsigned int> trial_interior_base_function_of_dof;
  std::vector<unsigned int> trial_skeleton_base_function_of_dof;
  classify_shape_functions(fe_test,
                           this->shape_function_type_test,
                           this->test_base_dofs,
                           test_base_function_of_dof);
  classify_shape_functions(fe_trial_interior,
                           this->shape_function_type_trial_interior,
                           this->trial_interior_base_dofs,
                           trial_interior_base_function_of_dof);
  classify_shape_functions(fe_trial_skeleton,
                           this->shape_function_type_trial_skeleton,
                           this->trial_skeleton_base_dofs,
                           trial_skeleton_base_function_of_dof);

  this->n_base_functions_test = this->test_base_dofs[E_real_field].size();
  this->n_base_functions_trial_interior =
    this->trial_interior_base_dofs[E_real_field].size();
  this->n_base_functions_trial_skeleton =
    this->trial_skeleton_base_dofs[E_real_field].size();

  this->all_test_base_functions.resize(n_base_functions_test);
  for (unsigned int a = 0; a < n_base_functions_test; ++a)
    this->all_test_base_functions[a] = a;

  // Base functions associated with each face of the reference cell, for the
  // test and skeleton trial spaces
  gather_face_base_functions(fe_test,
                             test_base_function_of_dof,
                             this->face_test_base_functions);
  gather_face_base_functions(fe_trial_skeleton,
                             trial_skeleton_base_function_of_dof,
                             this->face_trial_skeleton_base_functions);

  // Cell quadrature
  this->JxW               = std::vector<double>(n_q_points);
  this->quadrature_points = std::vector<Point<dim>>(n_q_points);

  // Temperature and effective properties at the cell quadrature points. The
  // temperature is the only field on which the electromagnetic properties can
  // depend.
  this->temperature_values = std::vector<double>(n_q_points, 0.);
  this->fields.clear();
  this->fields.insert(
    std::pair<field, std::vector<double>>(field::temperature, n_q_points));
  this->effective_electric_permittivities =
    std::vector<std::complex<double>>(n_q_points);
  this->effective_magnetic_permeabilities =
    std::vector<std::complex<double>>(n_q_points);

  // Base functions at the cell quadrature points
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

  // Temperature and effective properties at the face quadrature points
  this->face_temperature_values = std::vector<double>(n_face_q_points, 0.);
  this->face_fields.clear();
  this->face_fields.insert(
    std::pair<field, std::vector<double>>(field::temperature, n_face_q_points));
  this->face_effective_electric_permittivities =
    std::vector<std::complex<double>>(n_face_q_points);
  this->face_effective_magnetic_permeabilities =
    std::vector<std::complex<double>>(n_face_q_points);

  // Base functions at the face quadrature points
  this->phi_test_face.reinit(n_base_functions_test, n_face_q_points);
  this->n_cross_phi_test_face.reinit(n_base_functions_test, n_face_q_points);
  this->tangential_phi_trial_skeleton.reinit(n_base_functions_trial_skeleton,
                                             n_face_q_points);
  this->n_cross_tangential_phi_trial_skeleton.reinit(
    n_base_functions_trial_skeleton, n_face_q_points);
}

template <int dim>
void
TimeHarmonicMaxwellScratchData<dim>::classify_shape_functions(
  const FiniteElement<dim>                 &fe,
  std::vector<unsigned char>               &shape_function_type,
  std::array<std::vector<unsigned int>, 4> &base_dofs,
  std::vector<unsigned int>                &base_function_of_dof) const
{
  // Since the four fields use the same base element, each field has a quarter
  // of the shape functions of the finite element system
  const unsigned int n_base_functions = fe.n_dofs_per_cell() / 4;

  shape_function_type.assign(fe.n_dofs_per_cell(), 0);
  base_function_of_dof.assign(fe.n_dofs_per_cell(), 0);
  for (auto &dofs : base_dofs)
    dofs.assign(n_base_functions, numbers::invalid_unsigned_int);

  for (unsigned int k = 0; k < fe.n_dofs_per_cell(); ++k)
    {
      if (fe.shape_function_belongs_to(k, extractor_E_real))
        shape_function_type[k] |= electric_real;
      if (fe.shape_function_belongs_to(k, extractor_E_imag))
        shape_function_type[k] |= electric_imag;
      if (fe.shape_function_belongs_to(k, extractor_H_real))
        shape_function_type[k] |= magnetic_real;
      if (fe.shape_function_belongs_to(k, extractor_H_imag))
        shape_function_type[k] |= magnetic_imag;

      // Field of the shape function
      unsigned int field_index = E_real_field;
      if (shape_function_type[k] == electric_imag)
        field_index = E_imag_field;
      else if (shape_function_type[k] == magnetic_real)
        field_index = H_real_field;
      else if (shape_function_type[k] == magnetic_imag)
        field_index = H_imag_field;

      // Base function index: position of the shape function in its base
      // element, and copy of the base element for vector fields made of
      // scalar elements (e.g., FE_DGQ^dim for the interior trial space)
      const auto [base_element_and_copy, index_in_base_element] =
        fe.system_to_base_index(k);
      const unsigned int base_function =
        base_element_and_copy.second *
          fe.base_element(base_element_and_copy.first).n_dofs_per_cell() +
        index_in_base_element;

      AssertIndexRange(base_function, n_base_functions);
      base_dofs[field_index][base_function] = k;
      base_function_of_dof[k]               = base_function;
    }

  for (const auto &dofs : base_dofs)
    for (const unsigned int k : dofs)
      {
        (void)k;
        Assert(k != numbers::invalid_unsigned_int,
               ExcMessage("The four fields of the time-harmonic Maxwell finite "
                          "element systems must use the same base element."));
      }
}

template <int dim>
void
TimeHarmonicMaxwellScratchData<dim>::gather_face_base_functions(
  const FiniteElement<dim>               &fe,
  const std::vector<unsigned int>        &base_function_of_dof,
  std::vector<std::vector<unsigned int>> &face_base_functions) const
{
  const unsigned int n_faces = fe.reference_cell().n_faces();
  face_base_functions.assign(n_faces, std::vector<unsigned int>());

  for (unsigned int face = 0; face < n_faces; ++face)
    {
      // The degrees of freedom located on the face belong to the four fields,
      // so each base function appears four times
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

      // We check if the physical properties depend on the temperature field.
      // If so, the solver will evaluate the temperature at the quadrature
      // points using reinit_temperature.
      this->cell_material_needs_temperature =
        properties_manager.get_electric_conductivity(0, material_id)
          ->depends_on(field::temperature) ||
        properties_manager.get_electric_permittivity_real(0, material_id)
          ->depends_on(field::temperature) ||
        properties_manager.get_electric_permittivity_imag(0, material_id)
          ->depends_on(field::temperature) ||
        properties_manager.get_magnetic_permeability_real(0, material_id)
          ->depends_on(field::temperature) ||
        properties_manager.get_magnetic_permeability_imag(0, material_id)
          ->depends_on(field::temperature);
      std::ranges::fill(this->temperature_values, 0.);

      this->fe_values_trial_interior.reinit(cell);
      this->fe_values_test.reinit(cell_test);

      this->quadrature_points =
        this->fe_values_trial_interior.get_quadrature_points();

      for (unsigned int q = 0; q < n_q_points; ++q)
        this->JxW[q] = this->fe_values_trial_interior.JxW(q);

      // Real base functions of the test space and their curl. They are
      // evaluated with the E_real shape functions since the four fields share
      // the same base functions.
      for (unsigned int a = 0; a < n_base_functions_test; ++a)
        {
          const unsigned int k = this->test_base_dofs[E_real_field][a];
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              this->phi_test(a, q) =
                this->fe_values_test[extractor_E_real].value(k, q);
              this->curl_phi_test(a, q) =
                this->fe_values_test[extractor_E_real].curl(k, q);
            }
        }

      // Real base functions of the interior trial space
      for (unsigned int c = 0; c < n_base_functions_trial_interior; ++c)
        {
          const unsigned int k =
            this->trial_interior_base_dofs[E_real_field][c];
          for (unsigned int q = 0; q < n_q_points; ++q)
            this->phi_trial_interior(c, q) =
              this->fe_values_trial_interior[extractor_E_real].value(k, q);
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

      // Real base functions of the test space. On interior faces, only the
      // tangential traces of the base functions associated with the face are
      // needed, while the Robin boundary conditions need the full trace of all
      // the base functions on boundary faces.
      const std::vector<unsigned int> &test_base_functions =
        this->face_at_boundary ? this->all_test_base_functions :
                                 this->face_test_base_functions[face_no];
      for (const unsigned int a : test_base_functions)
        {
          const unsigned int k = this->test_base_dofs[E_real_field][a];
          for (unsigned int q = 0; q < n_face_q_points; ++q)
            {
              const Tensor<1, dim> value =
                this->fe_face_values_test[extractor_E_real].value(k, q);
              this->phi_test_face(a, q) = value;
              this->n_cross_phi_test_face(a, q) =
                cross_product_3d(this->face_normals[q], value);
            }
        }

      // Real base functions of the skeleton trial space associated with the
      // face. To be in H^-1/2(curl), the fields need to have the tangential
      // trace mapping (n x (E x n)) which extracts the tangential component
      // of the field at the face. Strictly speaking, n x E_parallel = n x E,
      // and we would not need to use the map_H12 function for the cross
      // products, but we keep it for consistency.
      for (const unsigned int d :
           this->face_trial_skeleton_base_functions[face_no])
        {
          const unsigned int k =
            this->trial_skeleton_base_dofs[E_real_field][d];
          for (unsigned int q = 0; q < n_face_q_points; ++q)
            {
              const Tensor<1, dim> tangential_value = map_H12(
                this->fe_face_values_trial_skeleton[extractor_E_real].value(k,
                                                                            q),
                this->face_normals[q]);
              this->tangential_phi_trial_skeleton(d, q) = tangential_value;
              this->n_cross_tangential_phi_trial_skeleton(d, q) =
                cross_product_3d(this->face_normals[q], tangential_value);
            }
        }
    }
}

template <int dim>
void
TimeHarmonicMaxwellScratchData<dim>::calculate_physical_properties()
{
  set_field_vector(field::temperature, this->temperature_values, this->fields);
  compute_effective_electromagnetic_properties(
    this->properties_manager,
    this->fields,
    this->material_id,
    this->effective_electric_permittivities,
    this->effective_magnetic_permeabilities);
}

template <int dim>
void
TimeHarmonicMaxwellScratchData<dim>::calculate_face_physical_properties()
{
  set_field_vector(field::temperature,
                   this->face_temperature_values,
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
