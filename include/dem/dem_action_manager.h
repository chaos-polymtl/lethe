#ifndef dem_action_manager_h
#define dem_action_manager_h

class DEMActionManager
{
public:
  /**
   * @brief Copy constructor as a delete function to make sure it can not be
   * copied. Will never be used.
   *
   * @param copy The object to be copied
   */
  DEMActionManager(const DEMActionManager &copy) = delete;

  /**
   * @brief Copy constructor as a delete function to make sure it can not be
   * assigned. Will never be used.
   *
   * @param copy The object to be copied
   */
  DEMActionManager &
  operator=(const DEMActionManager &copy) = delete;

  /**
   * @brief Getter of the unique instance of the DEMActionManager.
   */
  static DEMActionManager *
  get_action_manager();

  /**
   * @brief Reset all triggers to false.
   * TODO
   */
  inline void
  reset_triggers()
  {
    // First mobility status identification of the CFD time step (from the
    // velocity computed at the first DEM time step (counter = 0) of the CFD
    // time step) The contact search is executed to make sure the mobility
    // status of cell match the particles that are in.
    repartition_trigger                  = false;
    contact_search_trigger               = mobility_status_reset_trigger;
    clear_tangential_overlap_trigger     = false;
    solid_object_search_trigger          = false;
    sparse_contacts_cells_update_trigger = false;
    read_checkpoint_trigger              = false;
    mobility_status_reset_trigger        = false;
  }

  /**
   * @brief Flag that there are periodic boundary in the simulation.
   */
  inline void
  set_periodic_boundaries_enabling()
  {
    this->periodic_boundaries_enabled = true;
  }

  /**
   * @brief Check if the periodic boundaries are enabled to perform some actions.
   *
   * @return True if the periodic boundaries are enabled.
   */
  inline bool
  check_periodic_boundaries_enabling()
  {
    return periodic_boundaries_enabled;
  }

  /**
   * @brief Flag that there are solid objects in the simulation. Already trigger
   * a solid object search in the background map.
   */
  inline void
  set_solid_objects_enabling()
  {
    this->solid_objects_enabled = true;

    // Allowing the first contact search of solid objects
    solid_object_search_trigger = true;
  }

  /**
   * @brief Check if the solid objects are enabled to perform some actions.
   *
   * @return True if the solid objects are enabled.
   */
  inline bool
  check_solid_objects_enabling()
  {
    return solid_objects_enabled;
  }

  /**
   * @brief Flag that the sparse contacts are enabled and trigger a
   * mobility status reset to initialize the containers.
   */
  inline void
  set_sparse_contacts_enabling()
  {
    this->sparse_contacts_enabled = true;

    // Allowing the full broad search without mobility status at first iteration
    mobility_status_reset_trigger = true;
  }

  /**
   * @brief Check if the sparse contacts are enabled to perform some actions.
   * Es
   *
   * @return True if the sparse contacts are enabled.
   */
  inline bool
  check_sparse_contacts_enabling()
  {
    return sparse_contacts_enabled;
  }

  /**
   * @brief Flag that the grid will motion in the simulation.
   */
  inline void
  set_grid_motion_enabling()
  {
    this->grid_motion_enabled = true;
  }

  /**
   * @brief Check if the solid objects are enabled to perform some actions.
   *
   * @return True if the solid objects are enabled.
   */
  inline bool
  check_grid_motion_enabling()
  {
    return grid_motion_enabled;
  }

  /**
   * @brief Set triggers for action to be performed in the current time step
   * because the simulation starts from restart files (previously called
   * checkpoint step).
   *
   * It triggers the reading of checkpoints, the contact search,
   * clearing the contact structures, the solid object mapping (if enabled), and
   * the update of the sparse contacts cells (if enabled).
   */
  inline void
  restart_simulation()
  {
    read_checkpoint_trigger          = true;
    contact_search_trigger           = true;
    clear_tangential_overlap_trigger = true;
    solid_object_search_trigger      = solid_objects_enabled ? true : false;
    sparse_contacts_cells_update_trigger =
      sparse_contacts_enabled ? true : false;
  }

  /**
   * @brief Set triggers for action to be performed in the current time step
   * because of a load balancing step.
   *
   * It triggers:
   * - load balancing: the used method has determined that the load needs to be
   *                   balanced
   * - contact search: cells and particles will wont be on the same processor,
   *                   so all the contact search needs to be performed
   * - clearing the contact structures: same as above, the contact structures
   *                                    need to be cleared prior of it
   * - solid object search (if enabled): the solid objects need to be mapped to
   *                                     the rebalanced cells
   * - update the sparse contacts cells (if enabled): the sparse contacts cells
   *                                 need to be updated since they only contain
   *                                 local and ghost cells
   */
  inline void
  load_balance_step()
  {
    repartition_trigger              = true;
    contact_search_trigger           = true;
    clear_tangential_overlap_trigger = true;
    solid_object_search_trigger      = solid_objects_enabled ? true : false;
    sparse_contacts_cells_update_trigger =
      sparse_contacts_enabled ? true : false;
  }

  /**
   * @brief Set trigger for the contact search to be performed in the current
   * time step because of a particle insertion step.
   *
   * It triggers:
   * - contact search: since there are new particles in the simulation, their
   *                   potential contacts need to be found
   */
  inline void
  particle_insertion_step()
  {
    contact_search_trigger = true;
  }

  /**
   * @brief Set trigger for the contact search to be performed in the current
   * time step because of a contact detection step (from the frequent or dynamic
   * contact detection method).
   *
   * It triggers:
   * - contact search: might be from a selected frequency for the contact search
   *                   (frequent method) or from the displacement evaluation
   *                   (dynamic)
   */
  inline void
  contact_detection_step()
  {
    contact_search_trigger = true;
  }

  /**
   * @brief Set triggers for action to be performed in the current time step
   * because of a solid object search step.
   *
   * It triggers the solid object search and the contact search.
   */
  inline void
  solid_objects_search_step()
  {
    solid_object_search_trigger = true;
    contact_search_trigger      = true;
  }

  /**
   * @brief Set trigger that the mobility status need to be reset to mobile at
   * the first DEM time step if the sparse contacts are enabled.
   *
   * It trigger:
   * - mobility status reset: all the mobility status of the cells need to be
   *                          mobile since the velocity of the particles are
   *                          not yet computed and the status can not be
   *                          evaluated. In the case of CFD-DEM, we force the
   *                          full computation of the particles at the first DEM
   *                          time step since hydrodynamic forces, which have
   *                          just been calculated, may change the behavior of
   *                          particles.
   */
  inline void
  first_dem_step()
  {
    // Only triggered if the sparse contacts are enabled
    mobility_status_reset_trigger = sparse_contacts_enabled ? true : false;
  }

  /**
   * @brief Set trigger for the contact search to be performed at the last DEM
   * iteration of the CFD iteration.
   *
   * It triggers:
   * - contact search: the contact search is performed in order to call the
   *                   sorting of particles in cells and for the update of the
   *                   particle reference locations prior of the void fraction
   *                   calculation.
   */
  inline void
  last_dem_of_cfd_iteration_step()
  {
    contact_search_trigger = true;
  }

  /**
   * @brief Check if checkpoint needs to be read.
   */
  inline bool
  check_restart_simulation()
  {
    return read_checkpoint_trigger;
  }

  /**
   * @brief Check if the repartitioning (load balancing) needs to be performed.
   */
  inline bool
  check_repartition() const
  {
    return repartition_trigger;
  }

  /**
   * @brief Check if the contact search needs to be performed.
   */
  inline bool
  check_contact_search()
  {
    return contact_search_trigger;
  }

  /**
   * @brief Check if the sparse contacts cells need to be updated.
   */
  inline bool
  check_update_sparse_contacts_cells()
  {
    return sparse_contacts_cells_update_trigger;
  }

  /**
   * @brief Check if the solid object has to be searched.
   */
  inline bool
  check_solid_object_search()
  {
    return solid_object_search_trigger;
  }

  /**
   * @brief Check if the tangential overlap history needs to be cleared.
   */
  inline bool
  check_clear_tangential_overlap()
  {
    return clear_tangential_overlap_trigger;
  }

  /**
   * @brief Check if the mobility status need to be reset to mobile.
   */
  inline bool
  check_mobility_status_reset()
  {
    return mobility_status_reset_trigger;
  }

private:
  /**
   * @brief Constructor of the DEMActionManager.
   */
  DEMActionManager()
    : periodic_boundaries_enabled(false)
    , solid_objects_enabled(false)
    , sparse_contacts_enabled(false)
    , grid_motion_enabled(false)
    , read_checkpoint_trigger(false)
    , repartition_trigger(false)
    , clear_tangential_overlap_trigger(false)
    , contact_search_trigger(true)
    , solid_object_search_trigger(false)
    , sparse_contacts_cells_update_trigger(false)
    , mobility_status_reset_trigger(false)
  {}

  /**
   * @brief Pointer to the unique of the DEMActionManager.
   */
  static DEMActionManager *instance;

  /**
   * @brief Flag for periodic boundaries in the simulation.
   */
  bool periodic_boundaries_enabled;

  /**
   * @brief Flag for solid objects in the simulation.
   */
  bool solid_objects_enabled;

  /**
   * @brief Flag for the enabling of the sparse contacts.
   */
  bool sparse_contacts_enabled;

  /**
   * @brief Flag for motion of grid in the simulation.
   * TODO add grid motion to action manager
   */
  bool grid_motion_enabled;

  /**
   * @brief Flag of the trigger for reading the checkpoint.
   */
  bool read_checkpoint_trigger;

  /**
   * @brief Flag of the trigger for load balancing execution.
   */
  bool repartition_trigger;

  /**
   * @brief Flag of the trigger for clearing the tangential overlap history.
   */
  bool clear_tangential_overlap_trigger;

  /**
   * @brief Flag of the trigger for the contact search.
   */
  bool contact_search_trigger;

  /**
   * @brief Flag of the trigger for the solid object search, mapping the solid
   * with the background triangulation.
   */
  bool solid_object_search_trigger;

  /**
   * @brief Flag of the trigger for the update of the sparse contacts cells.
   */
  bool sparse_contacts_cells_update_trigger;

  /**
   * @brief Flag of the trigger for the mobility status reset to mobile status.
   */
  bool mobility_status_reset_trigger;
};
#endif
