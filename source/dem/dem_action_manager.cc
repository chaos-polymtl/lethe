#include <dem/dem_action_manager.h>

DEMActionManager *DEMActionManager::instance = nullptr;

DEMActionManager *
DEMActionManager::get_action_manager()
{
  if (instance == nullptr)
    instance = new DEMActionManager();

  return instance;
}
