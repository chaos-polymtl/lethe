#include <dem/output_force_torque_calculation.h>

template <int dim>
void
write_forces_torques_output_results(
  std::map<unsigned int, Tensor<1, dim>> force_on_walls,
  std::map<unsigned int, Tensor<1, dim>> torque_on_walls,
  const ConditionalOStream &             pcout)
{
  for (const auto &it : force_on_walls)
    {
      pcout << "Boundary " << it.first << " :\n"
            << "Force = " << it.second
            << "\nMoment = " << torque_on_walls[it.first] << "\n\n";
    }
}

template void
write_forces_torques_output_results(
  std::map<unsigned int, Tensor<1, 2>> force_on_walls,
  std::map<unsigned int, Tensor<1, 2>> torque_on_walls,
  const ConditionalOStream &           pcout);

template void
write_forces_torques_output_results(
  std::map<unsigned int, Tensor<1, 3>> force_on_walls,
  std::map<unsigned int, Tensor<1, 3>> torque_on_walls,
  const ConditionalOStream &           pcout);
