// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_vtu_reader_h
#define lethe_vtu_reader_h

#include <map>
#include <string>
#include <vector>

/**
 * @brief Return the paths of the .vtu files that a .vtu or a .pvtu file
 * refers to.
 *
 * A .vtu file refers to itself. A .pvtu file refers to each of the .vtu files
 * listed as one of its pieces, whose paths are resolved relative to the
 * directory of the .pvtu file. The files are returned in the order in which
 * they are listed in the .pvtu file.
 *
 * @param[in] file_name Path to a .vtu or a .pvtu file.
 *
 * @return Paths of the .vtu files.
 */
std::vector<std::string>
get_vtu_file_names(const std::string &file_name);

/**
 * @brief Read the point coordinates and point data arrays of every piece of a
 * .vtu file written by deal.II, and append each of their components to a
 * vector of a map.
 *
 * The values of each piece are appended in the order of the pieces in the
 * file, so the vectors hold the points in the order in which they are stored
 * in the file. Only the formats written by deal.II are supported: the ascii
 * format and the zlib compressed and base64 encoded binary format of Float32
 * arrays. Pieces without points are skipped.
 *
 * @param[in,out] map Map in which the components are appended, using the keys
 * given by @p array_keys.
 *
 * @param[in] file_name Path to the .vtu file.
 *
 * @param[in] array_keys Map from the name of each array to read to the keys
 * under which its components are stored in @p map, one key per component. The
 * point coordinates are requested with the name "Points", which is the name of
 * the section that holds them in a .vtu file. An array may be given fewer keys
 * than it has components, in which case only its first components are read.
 * For instance, {"Points", {"x", "y"}} only reads the first two coordinates.
 */
void
fill_vectors_from_vtu(
  std::map<std::string, std::vector<double>>            &map,
  const std::string                                     &file_name,
  const std::map<std::string, std::vector<std::string>> &array_keys);

#endif
