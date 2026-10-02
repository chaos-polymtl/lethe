// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/time_harmonic_maxwell_scratch_data.h>

#include <algorithm>

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
  // belong to and gather their local indices in the corresponding lists. This
  // is only done once since it only depends on the finite elements. Since the
  // three spaces are [E_real, E_imag, H_real, H_imag] finite element systems,
  // every shape function belongs to exactly one field.
  classify_shape_functions(fe_test,
                           this->shape_function_type_test,
                           this->test_dofs_electric,
                           this->test_dofs_magnetic);
  classify_shape_functions(fe_trial_interior,
                           this->shape_function_type_trial_interior,
                           this->trial_interior_dofs_electric,
                           this->trial_interior_dofs_magnetic);
  classify_shape_functions(fe_trial_skeleton,
                           this->shape_function_type_trial_skeleton,
                           this->trial_skeleton_dofs_electric,
                           this->trial_skeleton_dofs_magnetic);

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

  // Test and interior trial functions at the cell quadrature points
  this->phi_F.reinit(n_q_points, n_dofs_test);
  this->phi_F_conj.reinit(n_q_points, n_dofs_test);
  this->curl_phi_F.reinit(n_q_points, n_dofs_test);
  this->curl_phi_F_conj.reinit(n_q_points, n_dofs_test);
  this->phi_I.reinit(n_q_points, n_dofs_test);
  this->phi_I_conj.reinit(n_q_points, n_dofs_test);
  this->curl_phi_I.reinit(n_q_points, n_dofs_test);
  this->curl_phi_I_conj.reinit(n_q_points, n_dofs_test);
  this->phi_E.reinit(n_q_points, n_dofs_trial_interior);
  this->phi_H.reinit(n_q_points, n_dofs_trial_interior);

  // Face quadrature
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

  // Test and skeleton trial functions at the face quadrature points
  this->phi_F_face.reinit(n_face_q_points, n_dofs_test);
  this->phi_F_face_conj.reinit(n_face_q_points, n_dofs_test);
  this->phi_I_face_conj.reinit(n_face_q_points, n_dofs_test);
  this->n_cross_phi_I_face.reinit(n_face_q_points, n_dofs_test);
  this->n_cross_phi_I_face_conj.reinit(n_face_q_points, n_dofs_test);
  this->phi_E_hat.reinit(n_face_q_points, n_dofs_trial_skeleton);
  this->n_cross_phi_E_hat.reinit(n_face_q_points, n_dofs_trial_skeleton);
  this->n_cross_phi_H_hat.reinit(n_face_q_points, n_dofs_trial_skeleton);
}

template <int dim>
void
TimeHarmonicMaxwellScratchData<dim>::classify_shape_functions(
  const FiniteElement<dim>   &fe,
  std::vector<unsigned char> &shape_function_type,
  std::vector<unsigned int>  &electric_dofs,
  std::vector<unsigned int>  &magnetic_dofs) const
{
  shape_function_type.assign(fe.n_dofs_per_cell(), 0);
  electric_dofs.clear();
  magnetic_dofs.clear();

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

      if (shape_function_type[k] & is_electric)
        electric_dofs.emplace_back(k);
      else if (shape_function_type[k] & is_magnetic)
        magnetic_dofs.emplace_back(k);
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
      static constexpr std::complex<double> imag{0., 1.};

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
        {
          this->JxW[q] = this->fe_values_trial_interior.JxW(q);

          // Electric test functions F
          for (const unsigned int k : this->test_dofs_electric)
            {
              phi_F[q][k] = fe_values_test[extractor_E_real].value(k, q) +
                            imag * fe_values_test[extractor_E_imag].value(k, q);
              phi_F_conj[q][k] =
                fe_values_test[extractor_E_real].value(k, q) -
                imag * fe_values_test[extractor_E_imag].value(k, q);

              curl_phi_F[q][k] =
                fe_values_test[extractor_E_real].curl(k, q) +
                imag * fe_values_test[extractor_E_imag].curl(k, q);
              curl_phi_F_conj[q][k] =
                fe_values_test[extractor_E_real].curl(k, q) -
                imag * fe_values_test[extractor_E_imag].curl(k, q);
            }

          // Magnetic test functions I
          for (const unsigned int k : this->test_dofs_magnetic)
            {
              phi_I[q][k] = fe_values_test[extractor_H_real].value(k, q) +
                            imag * fe_values_test[extractor_H_imag].value(k, q);
              phi_I_conj[q][k] =
                fe_values_test[extractor_H_real].value(k, q) -
                imag * fe_values_test[extractor_H_imag].value(k, q);

              curl_phi_I[q][k] =
                fe_values_test[extractor_H_real].curl(k, q) +
                imag * fe_values_test[extractor_H_imag].curl(k, q);
              curl_phi_I_conj[q][k] =
                fe_values_test[extractor_H_real].curl(k, q) -
                imag * fe_values_test[extractor_H_imag].curl(k, q);
            }

          // Interior electric (E) and magnetic (H) trial functions
          for (const unsigned int k : this->trial_interior_dofs_electric)
            {
              phi_E[q][k] =
                fe_values_trial_interior[extractor_E_real].value(k, q) +
                imag * fe_values_trial_interior[extractor_E_imag].value(k, q);
            }
          for (const unsigned int k : this->trial_interior_dofs_magnetic)
            {
              phi_H[q][k] =
                fe_values_trial_interior[extractor_H_real].value(k, q) +
                imag * fe_values_trial_interior[extractor_H_imag].value(k, q);
            }
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
      static constexpr std::complex<double> imag{0., 1.};

      const auto face        = cell_skeleton->face(face_no);
      this->face_at_boundary = face->at_boundary();
      this->face_boundary_id = this->face_at_boundary ?
                                 face->boundary_id() :
                                 numbers::invalid_boundary_id;
      std::ranges::fill(this->face_temperature_values, 0.);

      this->fe_face_values_test.reinit(cell_test, face_no);
      this->fe_face_values_trial_skeleton.reinit(cell_skeleton, face_no);

      for (unsigned int q = 0; q < n_face_q_points; ++q)
        {
          const Tensor<1, dim> &normal =
            this->fe_face_values_trial_skeleton.normal_vector(q);

          this->face_JxW[q] = this->fe_face_values_trial_skeleton.JxW(q);
          this->face_quadrature_points[q] =
            this->fe_face_values_trial_skeleton.quadrature_point(q);
          this->face_normals[q] = normal;

          // Electric test functions F
          for (const unsigned int k : this->test_dofs_electric)
            {
              phi_F_face[q][k] =
                fe_face_values_test[extractor_E_real].value(k, q) +
                imag * fe_face_values_test[extractor_E_imag].value(k, q);
              phi_F_face_conj[q][k] =
                fe_face_values_test[extractor_E_real].value(k, q) -
                imag * fe_face_values_test[extractor_E_imag].value(k, q);
            }

          // Magnetic test functions I
          for (const unsigned int k : this->test_dofs_magnetic)
            {
              phi_I_face_conj[q][k] =
                fe_face_values_test[extractor_H_real].value(k, q) -
                imag * fe_face_values_test[extractor_H_imag].value(k, q);

              n_cross_phi_I_face[q][k] = cross_product_3d(
                normal,
                fe_face_values_test[extractor_H_real].value(k, q) +
                  imag * fe_face_values_test[extractor_H_imag].value(k, q));
              n_cross_phi_I_face_conj[q][k] = cross_product_3d(
                normal,
                fe_face_values_test[extractor_H_real].value(k, q) -
                  imag * fe_face_values_test[extractor_H_imag].value(k, q));
            }

          // Skeleton trial functions. To be in H^-1/2(curl), the fields need
          // to have the tangential trace mapping (n x (E x n)) which extracts
          // the tangential component of the field at the face. Strictly
          // speaking, n x E_parallel = n x E, and we would not need to use the
          // map_H12 function for the cross products, but we keep it for
          // consistency.
          for (const unsigned int k : this->trial_skeleton_dofs_electric)
            {
              phi_E_hat[q][k] = map_H12(
                fe_face_values_trial_skeleton[extractor_E_real].value(k, q) +
                  imag *
                    fe_face_values_trial_skeleton[extractor_E_imag].value(k, q),
                normal);
              n_cross_phi_E_hat[q][k] =
                cross_product_3d(normal, phi_E_hat[q][k]);
            }
          for (const unsigned int k : this->trial_skeleton_dofs_magnetic)
            {
              n_cross_phi_H_hat[q][k] = cross_product_3d(
                normal,
                map_H12(
                  fe_face_values_trial_skeleton[extractor_H_real].value(k, q) +
                    imag *
                      fe_face_values_trial_skeleton[extractor_H_imag].value(k,
                                                                            q),
                  normal));
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
