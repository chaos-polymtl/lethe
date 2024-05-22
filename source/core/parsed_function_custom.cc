
#include <core/parsed_function_custom.h>


template <int n_components>
ParsedFunctionCustom<n_components>::ParsedFunctionCustom(const double h)
  : h(h)
  , ht()
  , vnames("")
  , expression("")
  , constants_list("")
  , initialized(false)
{
  for (unsigned int i = 0; i < n_components; ++i)
    ht[i][i] = h;
}

template <int n_components>
void
ParsedFunctionCustom<n_components>::declare_parameters(ParameterHandler &prm)
{
  Assert(n_components > 0, ExcZero());

  // The expression of the variables
  vnames = "c_";
  for (unsigned int i = 1; i < n_components; i++)
    vnames = vnames + ",c_" + std::to_string(i);
  std::cout << vnames << std::endl;
  prm.declare_entry(
    "Variable names",
    vnames,
    Patterns::Anything(),
    "The names of the variables as they will be used in the "
    "function, separated by commas. By default, the names of variables "
    "with which the function will be evaluated are the `c_i' (i from 0 to "
    "n_components-1) and `t' for time. You can then use these variable names "
    "in your function expression and they will be "
    "replaced by the values of these variables at which the function is "
    "currently evaluated. However, you can also choose a different set "
    "of names for the variables at which to evaluate your function "
    "expression.");

  // The expression of the function
  std::string expr = "0";
  for (unsigned int i = 1; i < n_components; ++i)
    expr += "; 0";
  prm.declare_entry(
    "Function expression",
    expr,
    Patterns::Anything(),
    "The formula that denotes the function you want to evaluate for "
    "particular values of the variables. This expression "
    "may contain any of the usual operations such as addition or "
    "multiplication, as well as all of the common functions such as "
    "`sin' or `cos'. In addition, it may contain expressions like "
    "`if(x>0, 1, -1)' where the expression evaluates to the second "
    "argument if the first argument is true, and to the third argument "
    "otherwise. "
    "If the function you are describing represents a vector-valued "
    "function with multiple components, then separate the expressions "
    "for individual components by a semicolon.");

  prm.declare_entry(
    "Function constants",
    "",
    Patterns::Anything(),
    "Sometimes it is convenient to use symbolic constants in the "
    "expression that describes the function, rather than having to "
    "use its numeric value everywhere the constant appears. These "
    "values can be defined using this parameter, in the form "
    "`var1=value1, var2=value2, ...'."
    "\n\n"
    "A typical example would be to set this runtime parameter to "
    "`pi=3.1415926536' and then use `pi' in the expression of the "
    "actual formula. (That said, for convenience this class actually "
    "defines both `pi' and `Pi' by default, but you get the idea.)");
}


template <int n_components>
void
ParsedFunctionCustom<n_components>::initialize()
{
  if (!initialized)
    throw ExcMessage(
      "The required fields for muParser are not initialized yet. "
      "This function should be called only after that is done.");
  std::vector<std::string> const_list =
    Utilities::split_string_list(constants_list, ',');
  std::map<std::string, double> constants;
  for (const auto &constant : const_list)
    {
      std::vector<std::string> this_c =
        Utilities::split_string_list(constant, '=');
      AssertThrow(this_c.size() == 2,
                  ExcMessage("The list of constants, <" + constants_list +
                             ">, is not a comma-separated list of "
                             "entries of the form 'name=value'."));
      constants[this_c[0]] = Utilities::string_to_double(this_c[1]);
    }
  // set pi and Pi as synonyms for the corresponding value. note that
  // this overrides any value a user may have given
  constants["pi"] = numbers::PI;
  constants["Pi"] = numbers::PI;

  const unsigned int nn = (Utilities::split_string_list(vnames)).size();
  switch (nn)
    {
      case n_components:
        // Time independent function
        this->::internal::FunctionParserCustom::ParserImplementation<
          n_components>::initialize(vnames,
                                    Utilities::split_string_list(expression,
                                                                 ';'),
                                    constants);
        break;
      case n_components + 1:
        // Time dependent function
        this->::internal::FunctionParserCustom::ParserImplementation<
          n_components>::initialize(vnames,
                                    Utilities::split_string_list(expression,
                                                                 ';'),
                                    constants,
                                    true);
        break;
      default:
        AssertThrow(false,
                    ExcMessage(
                      "The list of variables specified is <" + vnames +
                      "> which is a list of length " +
                      Utilities::int_to_string(nn) +
                      " but it has to be a list of length equal to" +
                      " either n_components (for a time-independent function)" +
                      " or n_components+1 (for a time-dependent function)."));
    }
}

template <int n_components>
void
ParsedFunctionCustom<n_components>::parse_parameters(ParameterHandler &prm)
{
  vnames         = prm.get("Variable names");
  expression     = prm.get("Function expression");
  constants_list = prm.get("Function constants");
  initialized    = true;
  initialize();
}

template <int n_components>
void
ParsedFunctionCustom<n_components>::initialize(const std::string vnames,
                                               const std::string expression,
                                               const std::string constants_list)
{
  this->vnames         = vnames;
  this->expression     = expression;
  this->constants_list = constants_list;
  initialized          = true;
  initialize();
}

template <int n_components>
void
ParsedFunctionCustom<n_components>::vector_value(
  const Tensor<1, n_components> &p,
  Tensor<1, n_components>       &values) const
{
  for (unsigned int i = 0; i < n_components; i++)
    values[i] = this->do_value(p, this->get_time(), i);
}

template <int n_components>
double
ParsedFunctionCustom<n_components>::value(const Tensor<1, n_components> &p,
                                          unsigned int comp) const
{
  return this->do_value(p, this->get_time(), comp);
}

template <int n_components>
Tensor<1, n_components>
ParsedFunctionCustom<n_components>::gradient(const Tensor<1, n_components> &p,
                                             unsigned int comp) const
{
  Tensor<1, n_components> grad;
  // FourthOrder:
  Tensor<1, n_components> q1, q2, q3, q4;
  for (unsigned int i = 0; i < n_components; ++i)
    {
      q2      = p + ht[i];
      q1      = q2 + ht[i];
      q3      = p - ht[i];
      q4      = q3 - ht[i];
      grad[i] = (-this->value(q1, comp) + 8 * this->value(q2, comp) -
                 8 * this->value(q3, comp) + this->value(q4, comp)) /
                (12 * h);
    }
  return grad;
}

template <int n_components>
void
ParsedFunctionCustom<n_components>::vector_gradient(
  const Tensor<1, n_components> &p,
  Tensor<2, n_components>       &gradients) const
{
  Point<n_components>     q1, q2, q3, q4;
  Tensor<1, n_components> v1{}, v2{}, v3{}, v4{};
  const double            h_inv_12 = 1. / (12 * h);
  for (unsigned int i = 0; i < n_components; ++i)
    {
      q2 = p + ht[i];
      q1 = q2 + ht[i];
      q3 = p - ht[i];
      q4 = q3 - ht[i];
      this->vector_value(q1, v1);
      this->vector_value(q2, v2);
      this->vector_value(q3, v3);
      this->vector_value(q4, v4);

      for (unsigned int comp = 0; comp < n_components; ++comp)
        gradients[comp][i] =
          (-v1[comp] + 8 * v2[comp] - 8 * v3[comp] + v4[comp]) * h_inv_12;
    }
}

template <int n_components>
void
ParsedFunctionCustom<n_components>::set_time(const double /*newtime*/)
{}

// Explicit instantiations
template class ParsedFunctionCustom<1>;
template class ParsedFunctionCustom<2>;
template class ParsedFunctionCustom<3>;
template class ParsedFunctionCustom<4>;
template class ParsedFunctionCustom<5>;
template class ParsedFunctionCustom<6>;
template class ParsedFunctionCustom<7>;
template class ParsedFunctionCustom<8>;
template class ParsedFunctionCustom<9>;
template class ParsedFunctionCustom<10>;
template class ParsedFunctionCustom<11>;
template class ParsedFunctionCustom<12>;
template class ParsedFunctionCustom<13>;
template class ParsedFunctionCustom<14>;
template class ParsedFunctionCustom<15>;
template class ParsedFunctionCustom<16>;
template class ParsedFunctionCustom<17>;
template class ParsedFunctionCustom<18>;
template class ParsedFunctionCustom<19>;
template class ParsedFunctionCustom<20>;
