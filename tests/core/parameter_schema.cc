// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the JSON schema of the parameters printed by
 * lethe-parameter-check --schema.
 *
 * The schema of the mesh parameters is printed and read back. It must be valid
 * JSON, identify the solver family and the dimension, describe each entry by
 * its value, default value, documentation and pattern in the order of its
 * declaration, and its values must be readable by a ParameterHandler. The
 * escaping of the strings written in the JSON format is also tested.
 */

// Deal.II
#include <deal.II/base/parameter_handler.h>

// Lethe
#include <core/parameters.h>
#include <core/utilities.h>

// Tests (with common definitions)
#include <../tests/tests.h>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <sstream>
#include <string>

/**
 * @brief Path in a property tree. The '/' separator is used, since the names
 * of the parameters may contain dots.
 */
boost::property_tree::ptree::path_type
tree_path(const std::string &path)
{
  return boost::property_tree::ptree::path_type(path, '/');
}

void
test_schema()
{
  deallog << "--- Schema ---" << std::endl;

  // A value different from the default is set, so that the value and the
  // default value of the entry can be told apart
  ParameterHandler prm;
  Parameters::Mesh<2>::declare_parameters(prm);
  prm.enter_subsection("mesh");
  prm.set("type", "gmsh");
  prm.leave_subsection();

  std::ostringstream schema;
  print_parameter_schema(schema, prm, "lethe-fluid", 2);

  // The schema must be valid JSON
  boost::property_tree::ptree tree;
  std::istringstream          schema_stream(schema.str());
  boost::property_tree::read_json(schema_stream, tree);

  // The versions depend on the build, so only their presence is reported
  deallog << "Lethe version present: "
          << tree.get_optional<std::string>("lethe_version").has_value()
          << std::endl;
  deallog << "deal.II version present: "
          << tree.get_optional<std::string>("dealii_version").has_value()
          << std::endl;
  deallog << "Solver family: " << tree.get<std::string>("solver_family")
          << std::endl;
  deallog << "Dimension: " << tree.get<unsigned int>("dimension") << std::endl;

  const boost::property_tree::ptree &mesh_type =
    tree.get_child(tree_path("parameters/mesh/type"));
  deallog << "mesh/type value: " << mesh_type.get<std::string>("value")
          << std::endl;
  deallog << "mesh/type default value: "
          << mesh_type.get<std::string>("default_value") << std::endl;
  deallog << "mesh/type documentation: "
          << mesh_type.get<std::string>("documentation") << std::endl;
  deallog << "mesh/type pattern: "
          << mesh_type.get<std::string>("pattern_description") << std::endl;

  // The entries are listed in the order of their declaration
  deallog << "First entries of the mesh subsection:";
  unsigned int n_entries = 0;
  for (const auto &[name, entry] :
       tree.get_child(tree_path("parameters/mesh")))
    {
      if (n_entries++ == 3)
        break;
      deallog << " \"" << name << "\"";
    }
  deallog << std::endl;

  // The values of the schema can be read back by a ParameterHandler
  std::ostringstream parameters_json;
  boost::property_tree::write_json(parameters_json,
                                   tree.get_child("parameters"));

  ParameterHandler prm_read;
  Parameters::Mesh<2>::declare_parameters(prm_read);
  std::istringstream parameters_stream(parameters_json.str());
  prm_read.parse_input_from_json(parameters_stream);
  prm_read.enter_subsection("mesh");
  deallog << "mesh/type read back: " << prm_read.get("type") << std::endl;
  prm_read.leave_subsection();
}

void
test_json_string()
{
  deallog << "--- JSON strings ---" << std::endl;
  deallog << to_json_string("plain") << std::endl;
  deallog << to_json_string("a \"quoted\" word") << std::endl;
  deallog << to_json_string("a\\path") << std::endl;
  deallog << to_json_string("two\nlines\tand a tab") << std::endl;
  deallog << to_json_string(std::string("bell\a")) << std::endl;
}

int
main()
{
  try
    {
      initlog();

      test_schema();
      test_json_string();
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
