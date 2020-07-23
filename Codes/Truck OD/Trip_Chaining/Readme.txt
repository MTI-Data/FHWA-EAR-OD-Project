The Python script evaluates the unlinked truck trips identified from the GPS data, and chains truck trips according to their trip ends and a set of additional rules. With trip chaining, the long-haul truck trips that are originally separated into several shorter trips could be recovered.

User-defined function list:
- prepare_trip_ends(trip_roster_file, colname_file):
	The function takes the unheaded trip roster file (trip_roster_file) and the table heading file (colname_file) as inputs,
and converts the trip start and trip end into geometry object for future processing.

- mp_filter(full_trip, stop_shp_buffer):
	The function takes the trip ends and stop buffer data (including truck parking lots and gas stations) 
and applies the multiprocessing module to filter trips ending near the defined stops.

- filter_trip(trip_roster_file, parking_shp_nad83, gas_shp_nad83, output_dir, colname_file):
	The main function takes the unheaded trip roster file (trip_roster_file), the table heading file (colname_file), the shapefiles 
of truck parking lots and gas stations (parking_shp_nad83, gas_shp_nad83), and the output directory (output_dir) as inputs. 
And the function filters all trips ending near the defined stops.

- mp_match(overlap_stop, input_list):
	The function takes the gas stations with both in and out flows (overlap_stop), and a list of items (input_list) as inputs, 
including trips starting from the filtered gas stations (start_clean) or ending at the filtered gas stations (end_clean), 
the minimum dwell time of the stop (min_stay), and the maximum dwell time of the stop (max_stay) considered for trip linking.
And applies the multiprocessing module to match trips ending at and starting from the same gas stations at a first-in-first-out basis.

- fifo(stop_type = 'gas', min_stay = 1, max_stay= 180):
	The function takes the stop type considered (gas stations), the minimum dwell time of the stop (min_stay), 
and the maximum dwell time of the stop (max_stay) considered for trip linking. Then the function matches trips ending at 
and starting from the same gas stations at a first-in-first-out basis.

- pre_process(trip_start_park, trip_end_park):
	The function takes trips starting from and ending at the truck parking locations and processes the date and time.

- truck_trip_matching_main(trip_start_park, trip_end_park):
	The function takes processed trip file, and matches trips ending at and starting from the same parking locations 
at a first-in-first-out basis.

- trip_roster_merged(trip_roster_file, colname_file, trip_chain, park_pair_file, gas_pair_file):
	The function takes the unheaded trip roster file (trip_roster_file), the table heading file (colname_file), a dummy variable for
trip chaining (trip_chain), the trip pair matched based on parking locations (park_pair_file) and gas stations (gas_pair_file) as inputs.
Then it finds out the unmatched trips and chained trip OD as output.