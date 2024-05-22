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
#include <deal.II/base/mu_parser_internal.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/thread_local_storage.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>

#include <muParser.h>
#include <muParserFixes.h>

#include <cmath>
#include <ctime>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <random>
#include <string>
#include <vector>

namespace internal
{
  namespace FunctionParserCustom
  {
    /**
     * @brief Rounds the value to the nearest integer
     * @param[in] val Value to round
     * @return Rounded value
     */
    int
    mu_round(double val);

    /**
     * @brief If-condition
     * @param[in] condition Boolean expression
     * @param[in] thenvalue Expression if true
     * @param[in] elsevalue Expression if false
     * @return Resulting value
     */
    double
    mu_if(double condition, double thenvalue, double elsevalue);

    /**
     * @brief Logical or
     * @param[in] left First parameter
     * @param[in] right Second parameter
     * @return 1 if either is true, 0 else
     */
    double
    mu_or(double left, double right);

    /**
     * @brief Logical and
     * @param[in] left First parameter
     * @param[in] right Second parameter
     * @return 1 if both are true, 0 else
     */
    double
    mu_and(double left, double right);

    /**
     * Rounds the value to the nearest integer (same as round)
     * @param[in] value Value to round
     * @return Rounded value
     */
    double
    mu_int(double value);

    /**
     * Rounds the value to the upper integer
     * @param[in] value Value to round
     * @return Rounded value
     */
    double
    mu_ceil(double value);

    /**
     * @brief Rounds the value to the lower integer
     * @param[in] value Value to round
     * @return Rounded value
     */
    double
    mu_floor(double value);

    /**
     * @brief Performs the cotangent operation
     * @param[in] value Value to process
     * @return Cotan(value)
     */
    double
    mu_cot(double value);

    /**
     * @brief Performs the cosecant operation
     * @param[in] value Value to process
     * @return Cosec(value)
     */
    double
    mu_csc(double value);

    /**
     * @brief Performs the secant operation
     * @param[in] value Value to process
     * @return Sec(value)
     */
    double
    mu_sec(double value);

    /**
     * @brief Performs the natural-logarithm operation
     * @param[in] value Value to process
     * @return std::log(value)
     */
    double
    mu_log(double value);

    /**
     * @brief Performs the a to the power b operation
     * @param[in] a Base
     * @param[in] b Exponent
     * @return std::pow(a,b)
     */
    double
    mu_pow(double a, double b);

    /**
     * @brief Performs the erf operation
     * @param[in] value Value to process
     * @return std::erf(value)
     */
    double
    mu_erf(double value);

    /**
     * @brief Performs the erfc operation
     * @param[in] value Value to process
     * @return std::erfc(value)
     */
    double
    mu_erfc(double value);

    /**
     * @brief Generate number in [0,1] with a specific seed
     * @param[in] seed Number to be able to reproduce random behavior between
     * runs
     * @return Random number
     */
    double
    mu_rand_seed(double seed);

    /**
     * @brief Generate number in [0,1]
     * @return Random number
     */
    double
    mu_rand();

    /**
     * @brief Get the array of all function names.
     * @return Available implemented functions
     */
    std::vector<std::string>
    get_function_names();

    DeclException2(ExcParseError,
                   int,
                   std::string,
                   << "Parsing Error at Column " << arg1
                   << ". The parser said: " << arg2);


    /**
     * @brief Base class used for thread safe instances of the parser.
     */
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

    /**
     * @brief Class that interfaces with muParser
     * @tparam n_components
     */
    template <int n_components>
    class ParserImplementation
    {
    public:
      ParserImplementation();

      virtual ~ParserImplementation() = default;

      /**
       * @brief Treat the input and initialize the muParser object
       * @param[in] vars Variable names
       * @param[in] expressions Expressions
       * @param[in] constants Constants
       * @param[in] time_dependent Boolean to define if the function is time
       * dependent
       */
      void
      initialize(const std::string                   &vars,
                 const std::vector<std::string>      &expressions,
                 const std::map<std::string, double> &constants,
                 const bool                           time_dependent = false);

      /**
       * @brief Actual initialization of the muParser object
       */
      void
      init_muparser() const;

      /**
       * @brief Evaluation of the function at a specific point
       * @param[in] p Evalution point
       * @param[in] time Current time
       * @param[in] component Component of interest
       * @return Value
       */
      double
      do_value(const dealii::Tensor<1, n_components> &p,
               const double                           time,
               unsigned int                           component) const;

      /**
       * @brief Evaluation of the function at a specific point for all components
       * @param[in] p Evaluation point
       * @param[in] time Current time
       * @param[out] values Output values
       */
      void
      do_all_values(const dealii::Tensor<1, n_components> &p,
                    const double                           time,
                    dealii::ArrayView<double>             &values) const;

      std::vector<std::string> expressions;

    private:
      /// Thread safe object for using muParser
      mutable dealii::Threads::ThreadLocalStorage<
        ::internal::FunctionParserCustom::ParserData>
        parser_data;

      /// Constants
      std::map<std::string, double> constants;

      /// Variable name
      std::vector<std::string> var_names;

      /// Keeps track if the muParser object is initialized
      bool initialized;

      /// Number of variables
      unsigned int n_vars;
    };
  } // namespace FunctionParserCustom
} // namespace internal

#endif // lethe_mu_parser_internal_h
