#include <dem/effective_properties.h>

DEMEffectiveProperties *DEMEffectiveProperties::instance = nullptr;

DEMEffectiveProperties *
DEMEffectiveProperties::get_effective_properties()
{
  if (instance == nullptr)
    instance = new DEMEffectiveProperties();

  return instance;
}
