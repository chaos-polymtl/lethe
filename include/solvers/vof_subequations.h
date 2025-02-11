// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_vof_subequations_h
#define lethe_vof_subequations_h

/**
 * @brief IDs associated to the different subequations solved in Lethe.
 */
enum class VOFSubequationsID : unsigned int
{
  /// VOF phase fraction gradient L2 projection
  phase_gradient_projection = 0,
  /// VOF curvature L2 projection
  curvature_projection = 1,
  /// VOF algebraic interface reinitialization
  algebraic_interface_reinitialization = 2
};

/**
 * @brief Get the VOFSubequationsID associated to a specific VOF subequation
 * from a string.
 *
 * @param[in] vof_subequation_name String corresponding to a particular VOF
 * subequation.
 *
 * @return VOFSubequationsID corresponding to @p vof_subequation_name.
 */
inline VOFSubequationsID
get_vof_subequation_id(const std::string &vof_subequation_name)
{
  if (vof_subequation_name == "VOF phase gradient L2 projection")
    return VOFSubequationsID::phase_gradient_projection;
  else if (vof_subequation_name == "VOF curvature L2 projection")
    return VOFSubequationsID::curvature_projection;
  else if (vof_subequation_name == "VOF algebraic interface reinitialization")
    return VOFSubequationsID::algebraic_interface_reinitialization;
  else
    throw(std::invalid_argument("Invalid VOF subequation name. Options are: \n"
                                " <VOF phase gradient L2 projection>\n"
                                " <VOF curvature L2 projection>\n"
                                " <VOF algebraic interface reinitialization>"));
}

#endif
