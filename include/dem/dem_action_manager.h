// SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_dem_action_manager_h
#define lethe_dem_action_manager_h

/**
 * @brief Manage the actions to be performed in the DEM solver and on the DEM
 * part of the coupling CFD-DEM solver.
 *
 * Any event that needs some following actions is handled by the action
 * manager. The event communicates to the action manager that it happened and
 * the action manager will change the trigger flags to the actions to be taken.
 * When this action is about to be performed, it checks with the action manager
 * if it should be triggered.
 * Many events can trigger the same actions, but when the action is about to
 * be performed, it only checks if it has to be performed, not what event have
 * happened.
 *
 * If you implement a new feature that needs resetting or contact search after
 * (or other actions), you should add a new "step" function to the action
 * manager that will set the trigger to the action to be performed.
 * If you implement new actions, you should add a new trigger flag and a new
 * check function to the action manager. Obviously, you should also have to
 * initialize the trigger in the constructor, and reset it in the reset
 * function called after a DEM iteration.
 * If the action depends on a enabling of a feature, you should add a new
 * function to set the enabled flag.
 */
class DEMActionManager
{
public:
  /**
   * @brief Copy constructor as a delete function to make sure it can not be
   * copied. It will never be used.
   *
   * @param copy The object to be copied
   */
  DEMActionManager(const DEMActionManager &copy) = delete;

  /**
   * @brief Copy constructor as a delete function to make sure it can not be
   * assigned. It will never be used.
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
   * TODO: modify the handling of the contact search at counter = 1 because of
   * the resetting of the mobility status.
   */
  inline void
  reset_triggers()
  {
    // First mobility status identification of the CFD time step (from the
    // velocity computed at the first DEM time step (counter = 0) of the CFD
    // time step) The contact search is executed to make sure the mobility
    // status of the cell matches the particles that are in.
    load_balance_trigger                 = false;
    contact_search_trigger               = mobility_status_reset_trigger;
    clear_tangential_displacement_trigger     = false;
    solid_object_search_trigger          = false;
    sparse_contacts_cells_update_trigger = false;
    read_checkpoint_trigger              = false;
    mobility_status_reset_trigger        = false;
  }

  /**
   * @brief Flag that there are periodic boundary in the simulation.
   */
  inline void
  set_periodic_boundaries_enabled()
  {
    this->periodic_boundaries_enabled = true;
  }

  /**
   * @brief Check if the periodic boundaries are enabled to perform some actions
   * that are not handled by triggers of the action manager.
   *
   * @return True if the periodic boundaries are enabled.
   */
  inline bool
  check_periodic_boundaries_enabled()
  {
    return periodic_boundaries_enabled;
  }

  /**
   * @brief Flag that there are solid objects in the simulation. Already trigger
   * a solid object search in the background map.
   */
  inline void
  set_solid_objects_enabled()
  {
    this->solid_objects_enabled = true;

    // Allowing the first contact search of solid objects
    solid_object_search_trigger = true;
  }

  /**
   * @brief Check if the solid objects are enabled to perform some actions
   * that are not handled by triggers of the action manager.
   *
   * @return True if the solid objects are enabled.
   */
  inline bool
  check_solid_objects_enabled()
  {
    return solid_objects_enabled;
  }

  /**
   * @brief Flag that the sparse contacts is enabled and trigger a
   * mobility status reset to initialize the containers.
   */
  inline void
  set_sparse_contacts_enabled()
  {
    this->sparse_contacts_enabled = true;

    // Allowing the full broad search without mobility status at first iteration
    mobility_status_reset_trigger = true;
  }

  /**
   * @brief Check if the sparse contacts is enabled to perform some actions
   * that are not handled by triggers of the action manager.
   *
   * @return True if the sparse contacts is enabled.
   */
  inline bool
  check_sparse_contacts_enabled()
  {
    return sparse_contacts_enabled;
  }

  /**
   * @brief Flag that the grid will motion in the simulation.
   */
  inline void
  set_grid_motion_enabled()
  {
    this->grid_motion_enabled = true;
  }

  /**
   * @brief Check if the solid objects are enabled to perform some actions
   * that are not handled by triggers of the action manager.
   *
   * @return True if the solid objects are enabled.
   */
  inline bool
  check_grid_motion_enabled()
  {
    return grid_motion_enabled;
  }

  /**
   * @brief Set triggers for actions to be performed in the current time step
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
    clear_tangential_displacement_trigger = true;
    solid_object_search_trigger      = solid_objects_enabled ? true : false;
    sparse_contacts_cells_update_trigger =
      sparse_contacts_enabled ? true : false;
  }

  /**
   * @brief Set triggers for actions to be performed in the current iteration
   * because it is a load balancing step.
   *
   * It triggers:
   * - load balancing: the used method has determined that the load needs to be
   *                   balanced
   * - contact search: cells and particles will not be on the same processor,
   *                   so all the contact search needs to be performed.
   * - clearing the tangential displacement: the tangential displacement history of the
   *                                    particles pairs (or particle-wall) are
   *                                    lost when particles are handled by
   *                                    another processor. For consistency
   *                                    reasons, all the tangential displacement
   *                                    history is reset for all particles.
   * - resizing the containers: the containers need to be resized to the new
   *                            number of local particles.
   * - solid object search (if enabled): the solid objects need to be mapped to
   *                                     the rebalanced cells.
   * - update the sparse contacts cells (if enabled): the sparse contacts cells
   *                                    need to be updated since they only
   *                                    contain local and ghost cells.
   */
  inline void
  load_balance_step()
  {
    load_balance_trigger             = true;
    contact_search_trigger           = true;
    clear_tangential_displacement_trigger = true;
    resize_containers_trigger        = true;
    solid_object_search_trigger      = solid_objects_enabled ? true : false;
    sparse_contacts_cells_update_trigger =
      sparse_contacts_enabled ? true : false;
  }

  /**
   * @brief Set trigger for the contact search to be performed in the current
   * time step because of a particle insertion step.
   *
   * It triggers:
   * - contact search: since there are new particles in the simulation, contact
   *                   search needs to be performed for all particles.
   */
  inline void
  particle_insertion_step()
  {
    contact_search_trigger = true;
  }

  /**
   * @brief Set trigger for the contact search to be performed in the current
   * time step because of a contact detection step.
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
   * @brief Set triggers for actions to be performed in the current time step
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
  first_dem_of_cfddem_iteration_step()
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
  last_dem_of_cfddem_iteration_step()
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
  check_load_balance() const
  {
    return load_balance_trigger;
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
   * @brief Check if the tangential displacement history needs to be cleared.
   */
  inline bool
  check_clear_tangential_displacement()
  {
    return clear_tangential_displacement_trigger;
  }

  /**
   * @brief Check if the containers need to be resized.
   */
  inline bool
  check_resize_containers()
  {
    return resize_containers_trigger;
  }

  /**
   * @brief Check if the mobility status need to be reset to mobile.
   */
  inline bool
  check_mobility_status_reset()
  {
    return mobility_status_reset_trigger;
  }

  /**
   * @brief Check if the default broad search functions should be used. It
   * depends if adaptive sparse contacts is enabled, if so, it is a step where
   * the mobility status are reset.
   */
  inline bool
  use_default_broad_search_functions()
  {
    return !sparse_contacts_enabled || mobility_status_reset_trigger;
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
    , load_balance_trigger(false)
    , clear_tangential_displacement_trigger(false)
    , resize_containers_trigger(true)
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
   */
  bool grid_motion_enabled;

  /**
   * @brief Flag of the trigger for reading the checkpoint.
   */
  bool read_checkpoint_trigger;

  /**
   * @brief Flag of the trigger for load balancing execution.
   */
  bool load_balance_trigger;

  /**
   * @brief Flag of the trigger for clearing the tangential displacement history.
   */
  bool clear_tangential_displacement_trigger;

  /**
   * @brief Flag of the trigger for resize the vector dependant of the local
   * particles.
   */
  bool resize_containers_trigger;

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
