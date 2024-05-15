/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 -  by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------
 */

#ifndef lethe_mu_parser_internal_h
#define lethe_mu_parser_internal_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/thread_local_storage.h>

#include <memory>
#include <string>
#include <vector>

namespace internal
{
  namespace FunctionParserCustom
  {
    int
    mu_round(double val);

    double
    mu_if(double condition, double thenvalue, double elsevalue);

    double
    mu_or(double left, double right);

    double
    mu_and(double left, double right);

    double
    mu_int(double value);

    double
    mu_ceil(double value);

    double
    mu_floor(double value);

    double
    mu_cot(double value);

    double
    mu_csc(double value);

    double
    mu_sec(double value);

    double
    mu_log(double value);

    double
    mu_pow(double a, double b);

    double
    mu_erf(double value);

    double
    mu_erfc(double value);

    double
    mu_rand_seed(double seed);

    double
    mu_rand();


    std::vector<std::string>
    get_function_names();

    DeclException2(ExcParseError,
                   int,
                   std::string,
                   << "Parsing Error at Column " << arg1
                   << ". The parser said: " << arg2);



    class muParserBase
    {
    public:
      virtual ~muParserBase() = default;
    };

    struct ParserData
    {
      ParserData() = default;

      ParserData(const ParserData &) = delete;

      std::vector<double> vars;

      std::vector<std::unique_ptr<muParserBase>> parsers;
    };

    template <int n_components>
    class ParserImplementation
    {
    public:
      ParserImplementation();

      virtual ~ParserImplementation() = default;

      virtual void
      initialize(const std::string                   &vars,
                 const std::vector<std::string>      &expressions,
                 const std::map<std::string, double> &constants,
                 const bool                           time_dependent = false);

      void
      init_muparser() const;

      double
      do_value(const dealii::Tensor<1, n_components> &p,
               const double                           time,
               unsigned int                           component) const;

      void
      do_all_values(const dealii::Tensor<1, n_components> &p,
                    const double                           time,
                    dealii::ArrayView<double>             &values) const;

      std::vector<std::string> expressions;

    private:
      mutable dealii::Threads::ThreadLocalStorage<
        ::internal::FunctionParserCustom::ParserData>
        parser_data;

      std::map<std::string, double> constants;


      std::vector<std::string> var_names;

      bool initialized;

      unsigned int n_vars;
    };
  } // namespace FunctionParserCustom
} // namespace internal


#endif // LETHE_MU_PARSER_INTERNAL_H
