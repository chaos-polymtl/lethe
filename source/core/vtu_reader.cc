// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <deal.II/base/config.h>

#include <core/utilities.h>
#include <core/vtu_reader.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/utilities.h>

#ifdef DEAL_II_WITH_ZLIB
#  include <boost/iostreams/device/array.hpp>
#  include <boost/iostreams/filter/zlib.hpp>
#  include <boost/iostreams/filtering_stream.hpp>
#endif

#include <cctype>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <optional>
#include <sstream>
#include <string_view>

using namespace dealii;

namespace
{
  /**
   * @brief Opening tag and content of an XML element.
   */
  struct XMLElement
  {
    /// Text of the opening tag, from its '<' to its '>'.
    std::string_view tag;

    /// Text between the opening and the closing tags, which is empty for an
    /// element that closes itself.
    std::string_view content;

    /// Position right after the end of the element in the searched text.
    std::size_t end;
  };

  /**
   * @brief Read the whole content of a file in a string.
   *
   * @param[in] file_name Path to the file.
   * @param[in] file_description Role of the file, used in the error message
   * raised if the file can't be read.
   *
   * @return Content of the file.
   */
  std::string
  read_file_content(const std::string &file_name,
                    const std::string &file_description)
  {
    check_file_exists(file_name, file_description);

    const std::uintmax_t file_size = std::filesystem::file_size(file_name);
    std::string          content(file_size, '\0');
    std::ifstream        file(file_name, std::ios::binary);
    file.read(content.data(), static_cast<std::streamsize>(file_size));
    AssertThrow(file.gcount() == static_cast<std::streamsize>(file_size),
                ExcMessage("The " + file_description + " <" + file_name +
                           "> could not be read entirely."));
    return content;
  }

  /**
   * @brief Return the value of an attribute of an XML tag.
   *
   * @param[in] tag Text of the tag, from its '<' to its '>'.
   * @param[in] attribute_name Name of the attribute.
   *
   * @return Value of the attribute, or an empty string if the tag does not
   * have this attribute.
   */
  std::string
  get_attribute(const std::string_view tag, const std::string &attribute_name)
  {
    const std::string pattern  = attribute_name + "=\"";
    std::size_t       position = tag.find(pattern);

    // The attribute name must follow a white space, otherwise the pattern is
    // the end of a longer attribute name.
    while (position != std::string_view::npos &&
           (position == 0 ||
            !std::isspace(static_cast<unsigned char>(tag[position - 1]))))
      position = tag.find(pattern, position + 1);

    if (position == std::string_view::npos)
      return "";

    const std::size_t value_begin = position + pattern.size();
    const std::size_t value_end   = tag.find('"', value_begin);
    AssertThrow(value_end != std::string_view::npos,
                ExcMessage("The value of the attribute \"" + attribute_name +
                           "\" is not closed in the tag " + std::string(tag) +
                           "."));
    return std::string(tag.substr(value_begin, value_end - value_begin));
  }

  /**
   * @brief Find the next XML element with a given name in a text.
   *
   * Elements with the same name are assumed to be disjoint, as they are in
   * .vtu and .pvtu files.
   *
   * @param[in] text Text in which the element is searched.
   * @param[in] element_name Name of the element.
   * @param[in] start Position in @p text from which the element is searched.
   * @param[in] file_name Path to the file that holds the text, used in the
   * error messages.
   *
   * @return The element, or std::nullopt if there is no such element after
   * @p start.
   */
  std::optional<XMLElement>
  find_element(const std::string_view text,
               const std::string     &element_name,
               const std::size_t      start,
               const std::string     &file_name)
  {
    const std::string opening   = "<" + element_name;
    std::size_t       tag_begin = text.find(opening, start);

    // The element name must be followed by a white space, a '/' or a '>',
    // otherwise the pattern is the start of a longer element name.
    while (tag_begin != std::string_view::npos)
      {
        const std::size_t after_name = tag_begin + opening.size();
        if (after_name < text.size() &&
            (std::isspace(static_cast<unsigned char>(text[after_name])) ||
             text[after_name] == '/' || text[after_name] == '>'))
          break;
        tag_begin = text.find(opening, after_name);
      }

    if (tag_begin == std::string_view::npos)
      return std::nullopt;

    const std::size_t tag_end = text.find('>', tag_begin);
    AssertThrow(tag_end != std::string_view::npos,
                ExcMessage("A <" + element_name +
                           "> tag is not closed in the file <" + file_name +
                           ">."));

    XMLElement element;
    element.tag = text.substr(tag_begin, tag_end + 1 - tag_begin);

    if (text[tag_end - 1] == '/')
      {
        element.end = tag_end + 1;
        return element;
      }

    const std::string closing     = "</" + element_name + ">";
    const std::size_t content_end = text.find(closing, tag_end + 1);
    AssertThrow(content_end != std::string_view::npos,
                ExcMessage("A <" + element_name +
                           "> element is not closed in the file <" + file_name +
                           ">."));
    element.content = text.substr(tag_end + 1, content_end - tag_end - 1);
    element.end     = content_end + closing.size();
    return element;
  }

  /**
   * @brief Return the data array of the point coordinates and the point data
   * arrays of a piece of a .vtu file, by name.
   *
   * The data array of the point coordinates, which has no name, is returned
   * with the name "Points". The point data arrays are returned with their
   * own name.
   *
   * @param[in] piece Content of the piece.
   * @param[in] file_name Path to the .vtu file, used in the error messages.
   *
   * @return Map from the name of each data array to its element.
   */
  std::map<std::string, XMLElement>
  get_data_arrays(const std::string_view piece, const std::string &file_name)
  {
    std::map<std::string, XMLElement> data_arrays;

    const std::optional<XMLElement> points =
      find_element(piece, "Points", 0, file_name);
    if (points.has_value())
      {
        const std::optional<XMLElement> coordinates =
          find_element(points->content, "DataArray", 0, file_name);
        if (coordinates.has_value())
          data_arrays.emplace("Points", *coordinates);
      }

    const std::optional<XMLElement> point_data =
      find_element(piece, "PointData", 0, file_name);
    if (point_data.has_value())
      {
        std::optional<XMLElement> data_array =
          find_element(point_data->content, "DataArray", 0, file_name);
        while (data_array.has_value())
          {
            data_arrays.emplace(get_attribute(data_array->tag, "Name"),
                                *data_array);
            data_array = find_element(point_data->content,
                                      "DataArray",
                                      data_array->end,
                                      file_name);
          }
      }

    return data_arrays;
  }

  /**
   * @brief Decode the values of a data array of a .vtu file.
   *
   * @param[in] data_array Element of the data array.
   * @param[in] vtk_file_tag Tag of the VTKFile element of the file, which
   * states how binary data are compressed and the order of their bytes.
   * @param[in] n_values Number of values held by the array, which is its
   * number of tuples times its number of components.
   * @param[in] file_name Path to the .vtu file, used in the error messages.
   *
   * @return Values of the data array.
   */
  std::vector<double>
  decode_data_array(const XMLElement      &data_array,
                    const std::string_view vtk_file_tag,
                    const std::size_t      n_values,
                    const std::string     &file_name)
  {
    const std::string format = get_attribute(data_array.tag, "format");

    if (format == "ascii")
      {
        std::istringstream  stream{std::string(data_array.content)};
        std::vector<double> values;
        values.reserve(n_values);
        double value;
        while (stream >> value)
          values.push_back(value);
        AssertThrow(values.size() == n_values,
                    ExcMessage("A data array of the file <" + file_name +
                               "> holds " + std::to_string(values.size()) +
                               " values instead of " +
                               std::to_string(n_values) + "."));
        return values;
      }

    // deal.II only writes binary data in the zlib compressed format, as Float32
    // values in the byte order of the machine that wrote them.
    const std::string unsupported_format_message =
      "The data arrays of the file <" + file_name +
      "> are not in a supported format. Only the .vtu files written by "
      "Lethe, which are in the ascii format or in the zlib compressed binary "
      "format of Float32 values, are supported.";
    AssertThrow(format == "binary" &&
                  get_attribute(data_array.tag, "type") == "Float32" &&
                  get_attribute(vtk_file_tag, "compressor") ==
                    "vtkZLibDataCompressor",
                ExcMessage(unsupported_format_message));

    const std::string header_type = get_attribute(vtk_file_tag, "header_type");
    AssertThrow(header_type.empty() || header_type == "UInt32",
                ExcMessage(unsupported_format_message));

#ifdef DEAL_II_WORDS_BIGENDIAN
    const std::string machine_byte_order = "BigEndian";
#else
    const std::string machine_byte_order = "LittleEndian";
#endif
    AssertThrow(get_attribute(vtk_file_tag, "byte_order") == machine_byte_order,
                ExcMessage("The binary data of the file <" + file_name +
                           "> are not stored in the byte order of this "
                           "machine (" +
                           machine_byte_order + ")."));

#ifdef DEAL_II_WITH_ZLIB
    // Remove the white spaces around and within the base64 encoded data.
    std::string encoded;
    encoded.reserve(data_array.content.size());
    for (const char character : data_array.content)
      if (!std::isspace(static_cast<unsigned char>(character)))
        encoded.push_back(character);

    const std::string truncated_message =
      "A data array of the file <" + file_name + "> is truncated.";

    // The compressed data starts with a header made of the number of blocks,
    // the size of a block before compression, the size of the last block
    // before compression and the size of each block after compression, all
    // stored as 32 bits unsigned integers. The header is base64 encoded
    // separately from the blocks, so its size must be known to split the
    // encoded data. The number of blocks is read first from the 8 first
    // characters, which encode the 6 first bytes of the header.
    constexpr std::size_t n_characters_of_first_bytes = 8;
    AssertThrow(encoded.size() >= n_characters_of_first_bytes,
                ExcMessage(truncated_message));
    const std::vector<unsigned char> first_bytes =
      Utilities::decode_base64(encoded.substr(0, n_characters_of_first_bytes));
    std::uint32_t n_blocks;
    AssertThrow(first_bytes.size() >= sizeof(n_blocks),
                ExcMessage(truncated_message));
    std::memcpy(&n_blocks, first_bytes.data(), sizeof(n_blocks));

    // Base64 encodes each group of 3 bytes in 4 characters, padding the last
    // group.
    const std::size_t n_header_values     = 3 + std::size_t(n_blocks);
    const std::size_t header_size         = n_header_values * sizeof(n_blocks);
    const std::size_t n_header_characters = 4 * ((header_size + 2) / 3);
    AssertThrow(encoded.size() >= n_header_characters,
                ExcMessage(truncated_message));
    const std::vector<unsigned char> header_bytes =
      Utilities::decode_base64(encoded.substr(0, n_header_characters));
    AssertThrow(header_bytes.size() >= header_size,
                ExcMessage(truncated_message));
    std::vector<std::uint32_t> header(n_header_values);
    std::memcpy(header.data(), header_bytes.data(), header_size);

    // A size of the last block equal to zero means that the last block is
    // full.
    const std::size_t block_size = header[1];
    const std::size_t last_block_size =
      (header[2] == 0) ? block_size : std::size_t(header[2]);
    const std::size_t uncompressed_size =
      (n_blocks == 0) ? 0 : (n_blocks - 1) * block_size + last_block_size;
    AssertThrow(uncompressed_size == n_values * sizeof(float),
                ExcMessage("A data array of the file <" + file_name +
                           "> holds " + std::to_string(uncompressed_size) +
                           " bytes instead of " +
                           std::to_string(n_values * sizeof(float)) + "."));

    encoded.erase(0, n_header_characters);
    const std::vector<unsigned char> compressed_data =
      Utilities::decode_base64(encoded);

    // Each block is a separate zlib stream, which is uncompressed in place in
    // the vector of values.
    std::vector<float> values(n_values);
    char              *destination = reinterpret_cast<char *>(values.data());
    const char *source = reinterpret_cast<const char *>(compressed_data.data());
    std::size_t compressed_offset = 0;
    for (std::size_t block = 0; block < n_blocks; ++block)
      {
        const std::size_t compressed_block_size = header[3 + block];
        AssertThrow(compressed_offset + compressed_block_size <=
                      compressed_data.size(),
                    ExcMessage(truncated_message));

        const std::size_t uncompressed_block_size =
          (block + 1 == n_blocks) ? last_block_size : block_size;

        // The block is valid if it uncompresses to the size stated by the
        // header. The zlib decompressor throws if the block is corrupted.
        bool is_block_valid = false;
        try
          {
            boost::iostreams::filtering_istream block_stream;
            block_stream.push(boost::iostreams::zlib_decompressor());
            block_stream.push(
              boost::iostreams::array_source(source + compressed_offset,
                                             compressed_block_size));
            const auto expected_n_bytes =
              static_cast<std::streamsize>(uncompressed_block_size);
            block_stream.read(destination, expected_n_bytes);
            is_block_valid = block_stream.gcount() == expected_n_bytes;
          }
        catch (const std::exception &)
          {
            is_block_valid = false;
          }
        AssertThrow(is_block_valid,
                    ExcMessage("A data array of the file <" + file_name +
                               "> could not be uncompressed."));

        destination += uncompressed_block_size;
        compressed_offset += compressed_block_size;
      }

    return std::vector<double>(values.begin(), values.end());
#else
    AssertThrow(false,
                ExcMessage("The file <" + file_name +
                           "> holds compressed data, which can only be read "
                           "if deal.II was configured with zlib."));
    return {};
#endif
  }
} // namespace

std::vector<std::string>
get_vtu_file_names(const std::string &file_name)
{
  const std::filesystem::path path(file_name);
  if (path.extension() == ".vtu")
    return {file_name};

  AssertThrow(path.extension() == ".pvtu",
              ExcMessage("The file <" + file_name +
                         "> is neither a .vtu nor a .pvtu file."));

  const std::string content = read_file_content(file_name, ".pvtu file");

  // The pieces of a .pvtu file are the .vtu files, whose paths are relative
  // to the directory of the .pvtu file.
  std::vector<std::string>  vtu_file_names;
  std::optional<XMLElement> piece =
    find_element(content, "Piece", 0, file_name);
  while (piece.has_value())
    {
      const std::string source = get_attribute(piece->tag, "Source");
      AssertThrow(!source.empty(),
                  ExcMessage("A piece of the file <" + file_name +
                             "> has no source file."));
      vtu_file_names.push_back((path.parent_path() / source).string());
      piece = find_element(content, "Piece", piece->end, file_name);
    }

  AssertThrow(!vtu_file_names.empty(),
              ExcMessage("The file <" + file_name +
                         "> does not list any .vtu file."));
  return vtu_file_names;
}

void
fill_vectors_from_vtu(
  std::map<std::string, std::vector<double>>            &map,
  const std::string                                     &file_name,
  const std::map<std::string, std::vector<std::string>> &array_keys)
{
  const std::string content = read_file_content(file_name, ".vtu file");

  const std::optional<XMLElement> vtk_file =
    find_element(content, "VTKFile", 0, file_name);
  AssertThrow(vtk_file.has_value() &&
                get_attribute(vtk_file->tag, "type") == "UnstructuredGrid",
              ExcMessage("The file <" + file_name +
                         "> is not a .vtu file of an unstructured grid."));

  std::optional<XMLElement> piece =
    find_element(vtk_file->content, "Piece", 0, file_name);
  while (piece.has_value())
    {
      const std::string n_points_string =
        get_attribute(piece->tag, "NumberOfPoints");
      AssertThrow(!n_points_string.empty(),
                  ExcMessage("A piece of the file <" + file_name +
                             "> does not state its number of points."));
      const int n_points = Utilities::string_to_int(n_points_string);
      AssertThrow(n_points >= 0,
                  ExcMessage("A piece of the file <" + file_name +
                             "> has a negative number of points."));

      if (n_points > 0)
        {
          const std::map<std::string, XMLElement> data_arrays =
            get_data_arrays(piece->content, file_name);

          for (const auto &[array_name, keys] : array_keys)
            {
              const auto data_array = data_arrays.find(array_name);
              AssertThrow(data_array != data_arrays.end(),
                          ExcMessage(
                            "The array \"" + array_name +
                            "\" is missing from a piece of the file <" +
                            file_name + ">."));

              const std::string n_components_string =
                get_attribute(data_array->second.tag, "NumberOfComponents");
              const unsigned int n_components =
                n_components_string.empty() ?
                  1 :
                  Utilities::string_to_int(n_components_string);
              AssertThrow(keys.size() <= n_components,
                          ExcMessage(
                            "The array \"" + array_name + "\" of the file <" +
                            file_name + "> has fewer components than the " +
                            std::to_string(keys.size()) + " requested."));

              const std::vector<double> values =
                decode_data_array(data_array->second,
                                  vtk_file->tag,
                                  std::size_t(n_points) * n_components,
                                  file_name);

              for (unsigned int component = 0; component < keys.size();
                   ++component)
                {
                  std::vector<double> &component_values = map[keys[component]];
                  component_values.reserve(component_values.size() + n_points);
                  for (int point = 0; point < n_points; ++point)
                    component_values.push_back(
                      values[std::size_t(point) * n_components + component]);
                }
            }
        }

      piece = find_element(vtk_file->content, "Piece", piece->end, file_name);
    }
}
