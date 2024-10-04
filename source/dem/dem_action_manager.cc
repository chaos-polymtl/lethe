// SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/dem_action_manager.h>

DEMActionManager *DEMActionManager::instance = nullptr;

DEMActionManager *
DEMActionManager::get_action_manager()
{
  if (instance == nullptr)
    instance = new DEMActionManager();

  return instance;
}
