// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/solid_objects_parameters.h>

namespace Parameters
{
  extern template class Nitsche<2>;
  extern template class Nitsche<3>;
  extern template class DEMSolidObjects<2>;
  extern template class DEMSolidObjects<3>;
} // namespace Parameters
