// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/vof_algebraic_interface_reinitialization.h>

template <int dim>
void
VOFAlgebraicInterfaceReinitialization<dim>::assemble_system_matrix_and_rhs()
{
  // Reinitialize system matrix and right-hand side (rhs)
  this->system_matrix = 0;
  this->system_rhs    = 0;

  //

}