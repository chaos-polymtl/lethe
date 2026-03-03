// SPDX-FileCopyrightText: Copyright (c) 2019, 2021 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/source_terms.h>

// Pre-compile the 2D and 3D Navier-Stokes solver to ensure that the library is
// valid before we actually compile the solver This greatly helps with debugging
template class SourceTerms::SourceTerm<2>;
template class SourceTerms::SourceTerm<3>;
