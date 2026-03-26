// SPDX-FileCopyrightText: Copyright (c) 2025-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_cls_subequations_h
#define lethe_cls_subequations_h

/**
 * @brief IDs associated to the different subequations solved in Lethe.
 */
enum class CLSSubequationsID : unsigned int
{
  /// CLS phase fraction gradient L2 projection
  phase_gradient_projection = 0,
  /// CLS curvature L2 projection
  curvature_projection = 1,
  /// CLS algebraic interface reinitialization
  algebraic_interface_reinitialization = 2
};

/**
 * @brief Get the CLSSubequationsID associated to a specific CLS subequation
 * from a string.
 *
 * @param[in] cls_subequation_name String corresponding to a particular CLS
 * subequation.
 *
 * @return CLSSubequationsID corresponding to @p cls_subequation_name.
 */
inline CLSSubequationsID
get_cls_subequation_id(const std::string &cls_subequation_name)
{
  if (cls_subequation_name == "CLS phase gradient L2 projection")
    return CLSSubequationsID::phase_gradient_projection;
  else if (cls_subequation_name == "CLS curvature L2 projection")
    return CLSSubequationsID::curvature_projection;
  else if (cls_subequation_name == "CLS PDE-based interface reinitialization")
    return CLSSubequationsID::algebraic_interface_reinitialization;
  else
    throw(std::invalid_argument("Invalid CLS subequation name. Options are: \n"
                                " <CLS phase gradient L2 projection>\n"
                                " <CLS curvature L2 projection>\n"
                                " <CLS PDE-based interface reinitialization>"));
}

#endif
