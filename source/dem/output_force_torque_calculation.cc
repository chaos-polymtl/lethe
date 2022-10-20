#include <dem/output_force_torque_calculation.h>

void
write_forces_torques_output_locally(
  std::unordered_map<unsigned int, Tensor<1, 3>> force_on_walls,
  std::unordered_map<unsigned int, Tensor<1, 3>> torque_on_walls)
{
  TableHandler table;

  for (const auto &it : force_on_walls)
    {
      table.add_value("B_id", it.first);
      table.add_value("Fx", force_on_walls[it.first][0]);
      table.add_value("Fy", force_on_walls[it.first][1]);
      table.add_value("Fz", force_on_walls[it.first][2]);
      table.add_value("Tx", torque_on_walls[it.first][0]);
      table.add_value("Ty", torque_on_walls[it.first][1]);
      table.add_value("Tz", torque_on_walls[it.first][2]);
    }
  table.set_precision("Fx", 9);
  table.set_precision("Fy", 9);
  table.set_precision("Fz", 9);
  table.set_precision("Tx", 9);
  table.set_precision("Ty", 9);
  table.set_precision("Tz", 9);

  table.write_text(std::cout);
}
