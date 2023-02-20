#include "GrainTrackerGG2.h"

registerMooseObject("yinglongApp", GrainTrackerGG2);

InputParameters
GrainTrackerGG2::validParams()
{
  InputParameters params = GrainTrackerGG::validParams();
  params.addClassDescription("Grain Tracker object for running reduced order parameter simulations "
                             "without grain coalescence.");
  params.addParam<bool>("remerge_grains", false, "Grain merge would be considered if true");
  params.addRequiredParam<UserObjectName>("euler_angle_provider",
                                          "Name of Euler angle provider user object");
  return params;
}

GrainTrackerGG2::GrainTrackerGG2(const InputParameters & parameters)
  : GrainTrackerGG(parameters),
    _remerge_grains(getParam<bool>("remerge_grains")),
    _euler(getUserObject<EulerAngleProvider>("euler_angle_provider"))
{
}

GrainTrackerGG2::~GrainTrackerGG2() {}

void
GrainTrackerGG2::trackGrains()
{
  TIME_SECTION("trackGrains", 3, "Tracking Grains");

  mooseAssert(!_first_time, "Track grains may only be called when _tracking_step > _t_step");

  // Used to track indices for which to trigger the new grain callback on (used on all ranks)
  auto _old_max_grain_id = _max_curr_grain_id;

  /**
   * Only the primary rank does tracking, the remaining ranks
   * wait to receive local to global indices from the primary.
   */
  if (_is_primary)
  {
    _inactive_grains_id.clear();
    // Reset Status on active unique grains
    std::vector<unsigned int> map_sizes(_maps_size);
    for (auto & grain : _feature_sets_old)
    {
      if (grain._status != Status::INACTIVE)
      {
        grain._status = Status::CLEAR;
        map_sizes[grain._var_index]++;
      }
    }

    // Print out stats on overall tracking changes per var_index
    if (_verbosity_level > 0)
    {
      _console << "\nGrain Tracker Status:";
      for (const auto map_num : make_range(_maps_size))
      {
        _console << "\nGrains active index " << map_num << ": " << map_sizes[map_num] << " -> "
                 << _feature_counts_per_map[map_num];
        if (map_sizes[map_num] > _feature_counts_per_map[map_num])
          _console << "--";
        else if (map_sizes[map_num] < _feature_counts_per_map[map_num])
          _console << "++";
      }
      _console << '\n' << std::endl;
    }

    // Before we track grains, lets sort them so that we get parallel consistent answers
    std::sort(_feature_sets.begin(), _feature_sets.end());

    /**
     * To track grains across time steps, we will loop over our unique grains and link each one up
     * with one of our new unique grains. The criteria for doing this will be to find the unique
     * grain in the new list with a matching variable index whose centroid is closest to this
     * unique grain.
     */
    std::vector<std::size_t> new_grain_index_to_existing_grain_index(_feature_sets.size(),
                                                                     invalid_size_t);

    for (const auto old_grain_index : index_range(_feature_sets_old))
    {
      auto & old_grain = _feature_sets_old[old_grain_index];

      if (old_grain._status == Status::INACTIVE) // Don't try to find matches for inactive grains
        continue;

      std::size_t closest_match_index = invalid_size_t;
      Real min_centroid_diff = std::numeric_limits<Real>::max();

      /**
       * The _feature_sets vector is constructed by _var_index so we can avoid looping over all
       * indices. We can quickly jump to the first matching index to reduce the number of
       * comparisons and terminate our loop when our variable index stops matching.
       */
      // clang-format off
      auto start_it =
          std::lower_bound(_feature_sets.begin(), _feature_sets.end(), old_grain._var_index,
                           [](const FeatureData & item, std::size_t var_index)
                           {
                             return item._var_index < var_index;
                           });
      // clang-format on

      // We only need to examine grains that have matching variable indices
      bool any_boxes_intersect = false;
      for (MooseIndex(_feature_sets)
               new_grain_index = std::distance(_feature_sets.begin(), start_it);
           new_grain_index < _feature_sets.size() &&
           _feature_sets[new_grain_index]._var_index == old_grain._var_index;
           ++new_grain_index)
      {
        auto & new_grain = _feature_sets[new_grain_index];

        /**
         * Don't try to do any matching unless the bounding boxes at least overlap. This is to avoid
         * the corner case of having a grain split and a grain disappear during the same time step!
         */
        if (new_grain.boundingBoxesIntersect(old_grain))
        {
          any_boxes_intersect = true;
          Real curr_centroid_diff = centroidRegionDistance(old_grain._bboxes, new_grain._bboxes);
          if (curr_centroid_diff <= min_centroid_diff)
          {
            closest_match_index = new_grain_index;
            min_centroid_diff = curr_centroid_diff;
          }
        }
      }

      if (_verbosity_level > 2 && !any_boxes_intersect)
        _console << "\nNo intersecting bounding boxes found while trying to match grain "
                 << old_grain;

      // found a match
      if (closest_match_index != invalid_size_t)
      {
        /**
         * It's possible that multiple existing grains will map to a single new grain (indicated by
         * finding multiple matches when we are building this map). This will happen any time a
         * grain disappears during this time step. We need to figure out the rightful owner in this
         * case and inactivate the old grain.
         */
        auto curr_index = new_grain_index_to_existing_grain_index[closest_match_index];
        if (curr_index != invalid_size_t)
        {
          // The new feature being competed for
          auto & new_grain = _feature_sets[closest_match_index];

          // The other old grain competing to match up to the same new grain
          auto & other_old_grain = _feature_sets_old[curr_index];

          auto centroid_diff1 = centroidRegionDistance(new_grain._bboxes, old_grain._bboxes);
          auto centroid_diff2 = centroidRegionDistance(new_grain._bboxes, other_old_grain._bboxes);

          auto & inactive_grain = (centroid_diff1 < centroid_diff2) ? other_old_grain : old_grain;

          inactive_grain._status = Status::INACTIVE;
          if (_verbosity_level > 0)
          {
            _console << COLOR_GREEN << "Marking Grain " << inactive_grain._id
                     << " as INACTIVE (variable index: " << inactive_grain._var_index << ")\n"
                     << COLOR_DEFAULT;
            if (_verbosity_level > 1)
              _console << inactive_grain;
          }

          /**
           * If the grain we just marked inactive was the one whose index was in the new grain
           * to existing grain map (other_old_grain). Then we need to update the map to point
           * to the new match winner.
           */
          if (&inactive_grain == &other_old_grain)
            new_grain_index_to_existing_grain_index[closest_match_index] = old_grain_index;
        }
        else
          new_grain_index_to_existing_grain_index[closest_match_index] = old_grain_index;
      }
    }

    // Mark all resolved grain matches
    for (const auto new_index : index_range(new_grain_index_to_existing_grain_index))
    {
      auto curr_index = new_grain_index_to_existing_grain_index[new_index];

      // This may be a new grain, we'll handle that case below
      if (curr_index == invalid_size_t)
        continue;

      mooseAssert(_feature_sets_old[curr_index]._id != invalid_id,
                  "Invalid ID in old grain structure");

      _feature_sets[new_index]._id = _feature_sets_old[curr_index]._id; // Transfer ID
      _feature_sets[new_index]._status = Status::MARKED;      // Mark the status in the new set
      _feature_sets_old[curr_index]._status = Status::MARKED; // Mark the status in the old set
    }

    /**
     * At this point we have should have only two cases left to handle:
     * Case 1: A grain in the new set who has an unset status (These are new grains, previously
     *         untracked) This case is easy to understand. Since we are matching up grains by
     *         looking at the old set and finding closest matches in the new set, any grain in
     *         the new set that isn't matched up is simply new since some other grain satisfied
     *         each and every request from the old set.
     *
     * Case 2: A grain in the old set who has an unset status (These are inactive grains that
     *         haven't been marked) We can only fall into this case when the very last grain on
     *         a given variable disappears during the current time step. In that case we never have
     *         a matching _var_index in the comparison loop above so that old grain never competes
     *         for any new grain which means it can't be marked inactive in the loop above.
     */
    // Case 1 (new grains in _feature_sets):
    for (const auto grain_num : index_range(_feature_sets))
    {
      auto & grain = _feature_sets[grain_num];

      // New Grain
      if (grain._status == Status::CLEAR)
      {
        /**
         * Now we need to figure out what kind of "new" grain this is. Is it a nucleating grain that
         * we're just barely seeing for the first time or is it a "splitting" grain. A grain that
         * gets pinched into two or more pieces usually as it is being absorbed by other grains or
         * possibly due to external forces. We have to handle splitting grains this way so as to
         * no confuse them with regular grains that just happen to be in contact in this step.
         *
         * Splitting Grain: An grain that is unmatched by any old grain
         *                  on the same order parameter with touching halos.
         *
         * Nucleating Grain: A completely new grain appearing somewhere in the domain
         *                   not overlapping any other grain's halo.
         *
         * To figure out which case we are dealing with, we have to make another pass over all of
         * the existing grains with matching variable indices to see if any of them have overlapping
         * halos.
         */

        // clang-format off
        auto start_it =
            std::lower_bound(_feature_sets.begin(), _feature_sets.end(), grain._var_index,
                             [](const FeatureData & item, std::size_t var_index)
                             {
                               return item._var_index < var_index;
                             });
        // clang-format on

        // Loop over matching variable indices
        for (MooseIndex(_feature_sets)
                 new_grain_index = std::distance(_feature_sets.begin(), start_it);
             new_grain_index < _feature_sets.size() &&
             _feature_sets[new_grain_index]._var_index == grain._var_index;
             ++new_grain_index)
        {
          auto & other_grain = _feature_sets[new_grain_index];

          // Splitting grain?
          if (grain_num != new_grain_index && // Make sure indices aren't pointing at the same grain
              other_grain._status == Status::MARKED && // and that the other grain is indeed marked
              other_grain.boundingBoxesIntersect(grain) && // and the bboxes intersect
              other_grain.halosIntersect(grain))           // and the halos also intersect
          // TODO: Inspect combined volume and see if it's "close" to the expected value
          {
            grain._id = other_grain._id;    // Set the duplicate ID
            grain._status = Status::MARKED; // Mark it

            if (_verbosity_level > 0)
              _console << COLOR_YELLOW << "Split Grain Detected #" << grain._id
                       << " (variable index: " << grain._var_index << ")\n"
                       << COLOR_DEFAULT;
            if (_verbosity_level > 1)
              _console << grain << other_grain;
          }
        }

        if (grain._var_index < _reserve_op_index)
        {
          /**
           * The "try-harder loop":
           * OK so we still have an extra grain in the new set that isn't matched up against the
           * old set and since the order parameter isn't reserved. We aren't really expecting a new
           * grain. Let's try to make a few more attempts to see if this is a split grain even
           * though it failed to match the criteria above. This might happen if the halo front is
           * advancing too fast!
           *
           * In this loop we'll make an attempt to match up this new grain to the old halos. If
           * adaptivity is happening this could fail as elements in the new set may be at a
           * different level than in the old set. If we get multiple matches, we'll compare the
           * grain volumes (based on elements, not integrated to choose the closest).
           *
           * Future ideas:
           * Look at the volume fraction of the new grain and overlay it over the volume fraction
           * of the old grain (would require more saved information, or an aux field hanging around
           * (subject to projection problems).
           */
          if (_verbosity_level > 1)
            _console << COLOR_YELLOW
                     << "Trying harder to detect a split grain while examining grain on variable "
                        "index "
                     << grain._var_index << '\n'
                     << COLOR_DEFAULT;

          std::vector<std::size_t> old_grain_indices;
          for (const auto old_grain_index : index_range(_feature_sets_old))
          {
            auto & old_grain = _feature_sets_old[old_grain_index];

            if (old_grain._status == Status::INACTIVE)
              continue;

            /**
             * Note that the old grains we are looking at will already be marked from the earlier
             * tracking phase. We are trying to see if this unmatched grain is part of a larger
             * whole. To do that we'll look at the halos across the time step.
             */
            if (grain._var_index == old_grain._var_index &&
                grain.boundingBoxesIntersect(old_grain) && grain.halosIntersect(old_grain))
              old_grain_indices.push_back(old_grain_index);
          }

          if (old_grain_indices.size() == 1)
          {
            grain._id = _feature_sets_old[old_grain_indices[0]]._id;
            grain._status = Status::MARKED;

            if (_verbosity_level > 0)
              _console << COLOR_YELLOW << "Split Grain Detected #" << grain._id
                       << " (variable index: " << grain._var_index << ")\n"
                       << COLOR_DEFAULT;
          }
          else if (old_grain_indices.size() > 1)
            _console
                << COLOR_RED << "Split Grain Likely Detected #" << grain._id
                << " Need more information to find correct candidate - contact a developer!\n\n"
                << COLOR_DEFAULT;
        }

        // Must be a nucleating grain (status is still not set)
        if (grain._status == Status::CLEAR)
        {
          // auto new_index = getNextUniqueID();
          auto new_index = getTopoRelaGrainID(grain); // by weipeng
          grain._id = new_index;          // Set the ID
          grain._status = Status::MARKED; // Mark it

          if (_verbosity_level > 0)
            _console << COLOR_YELLOW << "Nucleating Grain Detected "
                     << " (variable index: " << grain._var_index 
                     << ", grain index: " << grain._id << ")\n"
                     << COLOR_DEFAULT;
          if (_verbosity_level > 1)
            _console << grain;
        }
      }
    }

    createAdjacentIDVector(); // by weipeng

    if (_remerge_grains && _t_step > 2) // by weipeng
      mergeGrainsBasedMisorientation();

    // Case 2 (inactive grains in _feature_sets_old)
    for (auto & grain : _feature_sets_old)
    {
      if (grain._status == Status::CLEAR)
      {
        grain._status = Status::INACTIVE;
        _inactive_grains_id.push_back(grain._id);

        if (_verbosity_level > 0)
        {
          _console << COLOR_GREEN << "Marking Grain " << grain._id
                   << " as INACTIVE (variable index: " << grain._var_index << ")\n"
                   << COLOR_DEFAULT;
          if (_verbosity_level > 1)
            _console << grain;
        }
      }
    }

    if (_inactive_grains_id.size() > 0)
      _invalid_feature_map.clear();
    
    for (auto & id_grain : _inactive_grains_id)
      for (auto & grain_old : _feature_sets_old)
        if (grain_old._id == id_grain && grain_old._adjacent_id.size() > 0)
        {
          std::vector<unsigned int> invalid_adjacentID_sort = grain_old._adjacent_id;
          std::sort(invalid_adjacentID_sort.begin(), invalid_adjacentID_sort.end());

          _invalid_feature_map[invalid_adjacentID_sort].push_back(grain_old._id);
          _invalid_feature_map[invalid_adjacentID_sort].push_back(grain_old._var_index);
          break;
        }
  } // is_primary

  /*************************************************************
   ****************** COLLECTIVE WORK SECTION ******************
   *************************************************************/

  // Make IDs on all non-primary ranks consistent
  scatterAndUpdateRanks();

  // Build up an id to index map
  _communicator.broadcast(_max_curr_grain_id);
  buildFeatureIdToLocalIndices(_max_curr_grain_id);

  /**
   * Trigger callback for new grains
   */
  if (_old_max_grain_id < _max_curr_grain_id)
  {
    for (auto new_id = _old_max_grain_id + 1; new_id <= _max_curr_grain_id; ++new_id)
    {
      // Don't trigger the callback on the reserve IDs
      if (new_id >= _reserve_grain_first_index + _n_reserve_ops)
      {
        // See if we've been instructed to terminate with an error
        if (!_first_time && _error_on_grain_creation)
          mooseError(
              "Error: New grain detected and \"error_on_new_grain_creation\" is set to true");
        else
          newGrainCreated(new_id);
      }
    }
  }
}

void
GrainTrackerGG2::remapGrains()
{
  // Don't remap grains if the current simulation step is before the specified tracking step
  if (_t_step < _tracking_step)
    return;

  TIME_SECTION("remapGrains", 3, "Remapping Grains");

  if (_verbosity_level > 1)
    _console << "Running remap Grains\n" << std::endl;

  /**
   * Map used for communicating remap indices to all ranks
   * This map isn't populated until after the remap loop.
   * It's declared here before we enter the root scope
   * since it's needed by all ranks during the broadcast.
   */
  std::map<unsigned int, std::size_t> grain_id_to_new_var;

  // Items are added to this list when split grains are found
  std::list<std::pair<std::size_t, std::size_t>> split_pairs;

  /**
   * The remapping algorithm is recursive. We will use the status variable in each FeatureData
   * to track which grains are currently being remapped so we don't have runaway recursion.
   * To begin we need to clear all of the active (MARKED) flags (CLEAR).
   *
   * Additionally we need to record each grain's variable index so that we can communicate
   * changes to the non-root ranks later in a single batch.
   */
  if (_is_primary)
  {
    // Build the map to detect difference in _var_index mappings after the remap operation
    std::map<unsigned int, std::size_t> grain_id_to_existing_var_index;
    for (auto & grain : _feature_sets)
    {
      // Unmark the grain so it can be used in the remap loop
      grain._status = Status::CLEAR;

      grain_id_to_existing_var_index[grain._id] = grain._var_index;
    }

    // Make sure that all split pieces of any grain are on the same OP
    for (const auto i : index_range(_feature_sets))
    {
      auto & grain1 = _feature_sets[i];

      for (const auto j : index_range(_feature_sets))
      {
        auto & grain2 = _feature_sets[j];
        if (i == j)
          continue;

        // The first condition below is there to prevent symmetric checks (duplicate values)
        if (i < j && grain1._id == grain2._id)
        {
          split_pairs.push_front(std::make_pair(i, j));
          if (grain1._var_index != grain2._var_index)
          {
            if (_verbosity_level > 0)
              _console << COLOR_YELLOW << "Split Grain (#" << grain1._id
                       << ") detected on unmatched OPs (" << grain1._var_index << ", "
                       << grain2._var_index << ") attempting to remap to " << grain1._var_index
                       << ".\n"
                       << COLOR_DEFAULT;

            /**
             * We're not going to try very hard to look for a suitable remapping. Just set it to
             * what we want and hope it all works out. Make the GrainTrackerGG great again!
             */
            grain1._var_index = grain2._var_index;
            grain1._status |= Status::DIRTY;

            if (_remerge_grains)  // by weipeng
              grain_id_to_new_var.emplace_hint(
                  grain_id_to_new_var.end(),
                  std::pair<unsigned int, std::size_t>(grain1._id, grain1._var_index));
          }
        }
      }
    }

    /**
     * Loop over each grain and see if any grains represented by the same variable are "touching"
     */
    bool any_grains_remapped = false;
    bool grains_remapped;

    std::set<unsigned int> notify_ids;
    do
    {
      grains_remapped = false;
      notify_ids.clear();

      for (auto & grain1 : _feature_sets)
      {
        // We need to remap any grains represented on any variable index above the cuttoff
        if (grain1._var_index >= _reserve_op_index)
        {
          if (_verbosity_level > 0)
            _console << COLOR_YELLOW << "\nGrain #" << grain1._id
                     << " detected on a reserved order parameter #" << grain1._var_index
                     << ", remapping to another variable\n"
                     << COLOR_DEFAULT;

          for (const auto max : make_range(0, _max_remap_recursion_depth + 1))
            if (max < _max_remap_recursion_depth)
            {
              if (attemptGrainRenumber(grain1, 0, max))
                break;
            }
            else if (!attemptGrainRenumber(grain1, 0, max))
            {
              _console << std::flush;
              std::stringstream oss;
              oss << "Unable to find any suitable order parameters for remapping while working "
                  << "with Grain #" << grain1._id << ", which is on a reserve order parameter.\n"
                  << "\n\nPossible Resolutions:\n"
                  << "\t- Add more order parameters to your simulation (8 for 2D, 28 for 3D)\n"
                  << "\t- Increase adaptivity or reduce your grain boundary widths\n"
                  << "\t- Make sure you are not starting with too many grains for the mesh size\n";
              mooseError(oss.str());
            }

          grains_remapped = true;
        }

        for (auto & grain2 : _feature_sets)
        {
          // Don't compare a grain with itself and don't try to remap inactive grains
          if (&grain1 == &grain2)
            continue;

          if (grain1._var_index == grain2._var_index && // grains represented by same variable?
              grain1._id != grain2._id &&               // are they part of different grains?
              grain1.boundingBoxesIntersect(grain2) &&  // do bboxes intersect (coarse level)?
              grain1.halosIntersect(grain2))            // do they actually overlap (fine level)?
          {
            if (_verbosity_level > 0)
              _console << COLOR_YELLOW << "Grain #" << grain1._id << " intersects Grain #"
                       << grain2._id << " (variable index: " << grain1._var_index << ")\n"
                       << COLOR_DEFAULT;

            for (const auto max : make_range(0, _max_remap_recursion_depth + 1))
            {
              if (max < _max_remap_recursion_depth)
              {
                if (attemptGrainRenumber(grain1, 0, max))
                {
                  grains_remapped = true;
                  break;
                }
              }
              else if (!attemptGrainRenumber(grain1, 0, max) &&
                       !attemptGrainRenumber(grain2, 0, max))
              {
                notify_ids.insert(grain1._id);
                notify_ids.insert(grain2._id);
              }
            }
          }
        }
      }
      any_grains_remapped |= grains_remapped;
    } while (grains_remapped);

    if (!notify_ids.empty())
    {
      _console << std::flush;
      std::stringstream oss;
      oss << "Unable to find any suitable order parameters for remapping while working "
          << "with the following grain IDs:\n"
          << Moose::stringify(notify_ids, ", ", "", true) << "\n\nPossible Resolutions:\n"
          << "\t- Add more order parameters to your simulation (8 for 2D, 28 for 3D)\n"
          << "\t- Increase adaptivity or reduce your grain boundary widths\n"
          << "\t- Make sure you are not starting with too many grains for the mesh size\n";

      if (_tolerate_failure)
        mooseWarning(oss.str());
      else
        mooseError(oss.str());
    }

    // Verify that split grains are still intact
    for (auto & split_pair : split_pairs)
      if (_feature_sets[split_pair.first]._var_index != _feature_sets[split_pair.first]._var_index)
        mooseError("Split grain remapped - This case is currently not handled");

    /**
     * The remapping loop is complete but only on the primary process.
     * Now we need to build the remap map and communicate it to the
     * remaining processors.
     */
    for (auto & grain : _feature_sets)
    {
      mooseAssert(grain_id_to_existing_var_index.find(grain._id) !=
                      grain_id_to_existing_var_index.end(),
                  "Missing unique ID");

      auto old_var_index = grain_id_to_existing_var_index[grain._id];

      if (old_var_index != grain._var_index)
      {
        mooseAssert(static_cast<bool>(grain._status & Status::DIRTY), "grain status is incorrect");

        grain_id_to_new_var.emplace_hint(
            grain_id_to_new_var.end(),
            std::pair<unsigned int, std::size_t>(grain._id, grain._var_index));

        /**
         * Since the remapping algorithm only runs on the root process,
         * the variable index on the primary's grains is inconsistent from
         * the rest of the ranks. These are the grains with a status of
         * DIRTY. As we build this map we will temporarily switch these
         * variable indices back to the correct value so that all
         * processors use the same algorithm to remap.
         */
        grain._var_index = old_var_index;
        // Clear the DIRTY status as well for consistency
        grain._status &= ~Status::DIRTY;
      }
    }

    if (!grain_id_to_new_var.empty())
    {
      if (_verbosity_level > 1)
      {
        _console << "Final remapping tally:\n";
        for (const auto & remap_pair : grain_id_to_new_var)
          _console << "Grain #" << remap_pair.first << " var_index "
                   << grain_id_to_existing_var_index[remap_pair.first] << " -> "
                   << remap_pair.second << '\n';
        _console << "Communicating swaps with remaining processors..." << std::endl;
      }
    }
  } // root processor

  // Communicate the std::map to all ranks
  _communicator.broadcast(grain_id_to_new_var);

  // Perform swaps if any occurred
  if (!grain_id_to_new_var.empty())
  {
    // Cache for holding values during swaps
    std::vector<std::map<Node *, CacheValues>> cache(_n_vars);

    // Perform the actual swaps on all processors
    for (auto & grain : _feature_sets)
    {
      // See if this grain was remapped
      auto new_var_it = grain_id_to_new_var.find(grain._id);
      if (new_var_it != grain_id_to_new_var.end())
        swapSolutionValues(grain, new_var_it->second, cache, RemapCacheMode::FILL);
    }

    for (auto & grain : _feature_sets)
    {
      // See if this grain was remapped
      auto new_var_it = grain_id_to_new_var.find(grain._id);
      if (new_var_it != grain_id_to_new_var.end())
        swapSolutionValues(grain, new_var_it->second, cache, RemapCacheMode::USE);
    }

    _nl.solution().close();
    _nl.solutionOld().close();
    _nl.solutionOlder().close();

    _fe_problem.getNonlinearSystemBase().system().update();

    if (_verbosity_level > 1)
      _console << "Swaps complete" << std::endl;
  }
}

unsigned int
GrainTrackerGG2::getTopoRelaGrainID(const FeatureData & grain_i)
{
  std::vector<unsigned int> nucleat_adj_ID;
  for (auto & grain_j : _feature_sets)
    if (grain_j._status == Status::MARKED && grain_i.boundingBoxesIntersect(grain_j) && grain_i.halosIntersect(grain_j))
      nucleat_adj_ID.push_back(grain_j._id);

  std::sort(nucleat_adj_ID.begin(), nucleat_adj_ID.end());

  for (auto & grain_i : nucleat_adj_ID)
      _console << COLOR_RED << "getTopoRelaGrainID::nucleat_adj_ID: " << grain_i << std::endl;

  std::cout << "***********" << std::endl;

  auto iter = _invalid_feature_map.find(nucleat_adj_ID);
  if ( iter != _invalid_feature_map.end())
    return iter->second[0];
  else if (_feature_sets_old.size() > 0)
  {
    for (auto & grain_i : _feature_sets_old)
      if (nucleat_adj_ID == grain_i._adjacent_id)
        return grain_i._id;
  }
  else if (nucleat_adj_ID.size() > 0)
  {
    _console << COLOR_RED << "impart nucleation grains to adjacent grains " << nucleat_adj_ID[0] << std::endl;
    return nucleat_adj_ID[0];
  }

  // mooseError("Detecting incorrect grain nucleation");
  return 0;
}

void 
GrainTrackerGG2::mergeGrainsBasedMisorientation()
{
  Real misor_angle = 0;
  const Real & threshold_merge = 0.8; 
  
  for (const auto grain_num_i : index_range(_feature_sets))
  {
    if (_feature_sets[grain_num_i]._status == Status::INACTIVE)
      continue;

    auto & grain_i = _feature_sets[grain_num_i];
    EulerAngles angles_i = _euler.getEulerAngles(grain_i._id);
    for (const auto grain_num_j : index_range(grain_i._adjacent_id))
    {
      auto & grain_j = _feature_sets[grain_i._adjacent_id[grain_num_j]];

      if (grain_j._status == Status::INACTIVE || grain_i._id >= grain_j._id)
        continue;

      EulerAngles angles_j = _euler.getEulerAngles(grain_j._id);

      misor_angle = CalculateMisorientationAngle::calculateMisorientaion(angles_i, angles_j, _s_misoriTwin, "hcp").misor;

      if (misor_angle < threshold_merge)
      {
        _console << COLOR_YELLOW << "Grain #" << grain_i._id << " and Grain #" << grain_j._id
                 << " was merged (misor: " << misor_angle << ").\n"
                 << COLOR_DEFAULT;

        grain_j._id = grain_i._id;
      }
    }
  }
}
