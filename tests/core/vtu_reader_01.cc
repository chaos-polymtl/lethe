// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests get_vtu_file_names and fill_vectors_from_vtu, which read the
 * point coordinates and point data arrays of the .vtu and .pvtu files written
 * by deal.II.
 *
 * The files are written in the layout of the particle output of Lethe: points
 * that each carry a scalar diameter and a velocity vector. The test covers a
 * compressed binary file, an ascii file, a file with several pieces, one of
 * which is empty, as written by DataOutInterface::write_vtu_in_parallel, and a
 * .pvtu file stored in a sub-directory. It also covers reading fewer
 * components than an array has and requesting an array that the file does not
 * hold.
 */

// Deal.II
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>

// Lethe
#include <core/vtu_reader.h>

// Tests (with common definitions)
#include <../tests/tests.h>

#include <array>
#include <filesystem>
#include <fstream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

/**
 * @brief Build one patch per point, holding the diameter and the velocity of
 * the point.
 *
 * @param[in] points Coordinates of the points.
 * @param[in] diameters Diameter of each point.
 * @param[in] velocities Velocity of each point.
 *
 * @return Patches of the points.
 */
std::vector<DataOutBase::Patch<0, 3>>
build_patches(const std::vector<Point<3>>              &points,
              const std::vector<double>                &diameters,
              const std::vector<std::array<double, 3>> &velocities)
{
  std::vector<DataOutBase::Patch<0, 3>> patches(points.size());
  for (unsigned int i = 0; i < points.size(); ++i)
    {
      patches[i].vertices[0] = points[i];
      patches[i].patch_index = i;
      patches[i].data.reinit(4, 1);
      patches[i].data(0, 0) = diameters[i];
      for (unsigned int d = 0; d < 3; ++d)
        patches[i].data(1 + d, 0) = velocities[i][d];
    }
  return patches;
}

/**
 * @brief Log the vectors of a map, one key per line.
 *
 * @param[in] case_name Name of the case, used to identify the lines of output.
 * @param[in] map Map whose vectors are logged.
 */
void
log_map(const std::string                                &case_name,
        const std::map<std::string, std::vector<double>> &map)
{
  deallog << case_name << std::endl;
  for (const auto &[key, values] : map)
    {
      deallog << "  " << key << ":";
      for (const double value : values)
        deallog << " " << value;
      deallog << std::endl;
    }
}

void
test()
{
  // Values that are exactly represented as Float32 and with the 6 digits of
  // the ascii format, so that every format gives the same output.
  const std::vector<Point<3>>              points = {Point<3>(0.5, 0.25, -0.75),
                                                     Point<3>(1.5, -2., 0.125),
                                                     Point<3>(3., 0.0625, 4.)};
  const std::vector<double>                diameters  = {0.25, 0.5, 0.75};
  const std::vector<std::array<double, 3>> velocities = {{{1., -1., 0.5}},
                                                         {{0.25, 2., -0.125}},
                                                         {{-3., 0.75, 1.25}}};
  const std::vector<DataOutBase::Patch<0, 3>> patches =
    build_patches(points, diameters, velocities);

  // Names of the data arrays, the three components of the velocity being
  // written as a single vector array, as in the particle output of Lethe.
  const std::vector<std::string> data_names = {"diameter",
                                               "velocity",
                                               "velocity",
                                               "velocity"};
  const std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
    nonscalar_data_ranges = {
      {1,
       3,
       "velocity",
       DataComponentInterpretation::component_is_part_of_vector}};

  const std::map<std::string, std::vector<std::string>> array_keys = {
    {"Points", {"p_x", "p_y", "p_z"}},
    {"diameter", {"diameters"}},
    {"velocity", {"v_x", "v_y", "v_z"}}};

  const std::string directory_name = "vtu_files";
  std::filesystem::create_directory(directory_name);

  // Compressed binary file, which is the default format of deal.II.
  const std::string compressed_file_name = directory_name + "/compressed.vtu";
  {
    std::ofstream         file(compressed_file_name);
    DataOutBase::VtkFlags flags;
    DataOutBase::write_vtu(
      patches, data_names, nonscalar_data_ranges, flags, file);
  }
  {
    std::map<std::string, std::vector<double>> map;
    fill_vectors_from_vtu(map, compressed_file_name, array_keys);
    log_map("Compressed file", map);
  }

  // Ascii file.
  const std::string ascii_file_name = directory_name + "/ascii.vtu";
  {
    std::ofstream         file(ascii_file_name);
    DataOutBase::VtkFlags flags;
    flags.compression_level = DataOutBase::CompressionLevel::plain_text;
    DataOutBase::write_vtu(
      patches, data_names, nonscalar_data_ranges, flags, file);
  }
  {
    std::map<std::string, std::vector<double>> map;
    fill_vectors_from_vtu(map, ascii_file_name, array_keys);
    log_map("Ascii file", map);
  }

  // File with three pieces, the second one being empty, as written by
  // DataOutInterface::write_vtu_in_parallel when a process has no points.
  const std::string pieces_file_name = "pieces.vtu";
  {
    const std::vector<DataOutBase::Patch<0, 3>> first_piece(patches.begin(),
                                                            patches.begin() +
                                                              1);
    const std::vector<DataOutBase::Patch<0, 3>> empty_piece;
    const std::vector<DataOutBase::Patch<0, 3>> last_piece(patches.begin() + 1,
                                                           patches.end());

    std::ofstream         file(pieces_file_name);
    DataOutBase::VtkFlags flags;
    DataOutBase::write_vtu_header(file, flags);
    DataOutBase::write_vtu_main(
      first_piece, data_names, nonscalar_data_ranges, flags, file);
    DataOutBase::write_vtu_main(
      empty_piece, data_names, nonscalar_data_ranges, flags, file);
    DataOutBase::write_vtu_main(
      last_piece, data_names, nonscalar_data_ranges, flags, file);
    DataOutBase::write_vtu_footer(file);
  }
  {
    std::map<std::string, std::vector<double>> map;
    fill_vectors_from_vtu(map, pieces_file_name, array_keys);
    log_map("File with several pieces", map);
  }

  // A .vtu file refers to itself.
  for (const auto &file_name : get_vtu_file_names(pieces_file_name))
    deallog << "File referred to by the .vtu file: " << file_name << std::endl;

  // .pvtu file listing the compressed and ascii files. The listed files are
  // relative to the directory of the .pvtu file.
  const std::string pvtu_file_name = directory_name + "/pieces.pvtu";
  {
    std::ofstream         file(pvtu_file_name);
    DataOutBase::VtkFlags flags;
    DataOutBase::write_pvtu_record(file,
                                   {"compressed.vtu", "ascii.vtu"},
                                   data_names,
                                   nonscalar_data_ranges,
                                   flags);
  }
  {
    std::map<std::string, std::vector<double>> map;
    for (const auto &file_name : get_vtu_file_names(pvtu_file_name))
      {
        deallog << "File referred to by the .pvtu file: " << file_name
                << std::endl;
        fill_vectors_from_vtu(map, file_name, array_keys);
      }
    log_map("Files of the .pvtu file", map);
  }

  // Only the first two coordinates are read, as in a 2D simulation.
  {
    std::map<std::string, std::vector<double>> map;
    fill_vectors_from_vtu(map,
                          compressed_file_name,
                          {{"Points", {"p_x", "p_y"}}});
    log_map("First two coordinates", map);
  }

  // An array that the file does not hold. The message of the exception is not
  // logged since it holds the line of the source file that raised it.
  try
    {
      std::map<std::string, std::vector<double>> map;
      fill_vectors_from_vtu(map,
                            compressed_file_name,
                            {{"temperature", {"T"}}});
      deallog << "Missing array: no exception raised" << std::endl;
    }
  catch (const ExceptionBase &)
    {
      deallog << "Missing array: exception raised" << std::endl;
    }

  std::filesystem::remove_all(directory_name);
  std::filesystem::remove(pieces_file_name);
}

int
main()
{
  try
    {
      initlog();
      test();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  return 0;
}
