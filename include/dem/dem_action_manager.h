

#ifndef LETHE_TRIGGER_MANAGER_H
#define LETHE_TRIGGER_MANAGER_H

class DEMActionManager
{
public:
  // move in .cc file
  static DEMActionManager *
  get_action_manager();

  // Make sure can not be clonable
  DEMActionManager(const DEMActionManager &copy) = delete;

  // Make sure can not be assigned
  DEMActionManager &
  operator=(const DEMActionManager &copy) = delete;

  inline void
  reset_triggers()
  {
    // First mobility status identification of the CFD time step (from the
    // velocity computed at the first DEM time step (counter = 0) of the CFD
    // time step) The contact search is executed to make sure the mobility
    // status of cell match the particles that are in.
    repartition_trigger                  = false;
    contact_search_trigger               = mobility_status_reset_trigger;
    clear_contact_structures_trigger     = false;
    solid_object_search_trigger          = false;
    sparse_contacts_cells_update_trigger = false;
    restart_simulation_trigger           = false;
    mobility_status_reset_trigger        = false;
  }


  inline void
  set_sparse_contacts_enabling()
  {
    this->sparse_contacts_enabled = true;

    // Allowing the full broad search without mobility status at first iteration
    mobility_status_reset_trigger = true;
  }


  inline bool
  check_sparse_contacts_enabling()
  {
    return sparse_contacts_enabled;
  }

  inline void
  set_periodic_boundaries_enabling()
  {
    this->periodic_boundaries_enabled = true;
  }

  inline bool
  check_periodic_boundaries_enabling()
  {
    return periodic_boundaries_enabled;
  }

  inline void
  set_solid_objects_enabling()
  {
    this->solid_objects_enabled = true;

    // Allowing the first contact search of solid objects
    solid_object_search_trigger = true;
  }

  inline bool
  check_solid_objects_enabling()
  {
    return solid_objects_enabled;
  }


  // Need something to say that the previous step was a contact detection
  inline void
  load_balance_step()
  {
    repartition_trigger              = true;
    contact_search_trigger           = true;
    clear_contact_structures_trigger = true;
    solid_object_search_trigger =
      solid_objects_enabled ? true : false; // Not check but should
    sparse_contacts_cells_update_trigger =
      sparse_contacts_enabled ? true : false;
  }

  inline void
  particle_insertion_step()
  {
    contact_search_trigger = true;
  }

  inline void
  contact_search_step()
  {
    contact_search_trigger = true;
  }


  inline void
  restart_simulation()
  {
    restart_simulation_trigger       = true;
    contact_search_trigger           = true;
    clear_contact_structures_trigger = true;
    solid_object_search_trigger      = solid_objects_enabled ? true : false;
    sparse_contacts_cells_update_trigger =
      sparse_contacts_enabled ? true : false;
  }

  inline bool
  check_restart_simulation()
  {
    return restart_simulation_trigger;
  }

  inline bool
  check_update_sparse_contacts_cells()
  {
    return sparse_contacts_cells_update_trigger;
  }

  inline bool
  check_solid_object_search()
  {
    return solid_object_search_trigger;
  }

  inline void
  solid_objects_search_step()
  {
    solid_object_search_trigger = true;
    contact_search_trigger      = true;
  }


  inline bool
  check_contact_search()
  {
    return contact_search_trigger;
  }

  inline bool
  check_repartition() const
  {
    return repartition_trigger;
  }

  inline bool
  check_clear_contact_structures()
  {
    return clear_contact_structures_trigger;
  }

  inline bool
  check_mobility_status_reset()
  {
    return mobility_status_reset_trigger;
  }

  inline void
  first_dem_step()
  {
    // Only triggered if the sparse contacts are enabled
    mobility_status_reset_trigger = sparse_contacts_enabled ? true : false;
  }

  inline void
  last_dem_iteration_with_pbc_step()
  {
    contact_search_trigger = true;
  }



private:
  DEMActionManager()
    : solid_objects_enabled(false)
    , sparse_contacts_enabled(false)
    , periodic_boundaries_enabled(false)
    , grid_motion_enabled(false)
    , repartition_trigger(false)
    , contact_search_trigger(true)
    , solid_object_search_trigger(false)
    , sparse_contacts_cells_update_trigger(false)
    , clear_contact_structures_trigger(false)
    , restart_simulation_trigger(false)
    , mobility_status_reset_trigger(false)
  {}

  static DEMActionManager *instance;



  // Enabled feature parameters
  bool solid_objects_enabled;
  bool sparse_contacts_enabled;
  bool periodic_boundaries_enabled;
  bool grid_motion_enabled;

  // Step action


  // Action trigger
  bool repartition_trigger;
  bool contact_search_trigger;
  bool solid_object_search_trigger;
  bool sparse_contacts_cells_update_trigger;
  bool clear_contact_structures_trigger;
  bool restart_simulation_trigger;
  bool mobility_status_reset_trigger;
};



#endif // LETHE_TRIGGER_MANAGER_H
