# -*- coding: utf-8 -*-


import pandas as pd
import geopandas as gpd
from shapely import geometry as sg
from functools import partial
import time
import gc
from multiprocessing import Pool
from datetime import timedelta
import warnings
import dateutil.parser
import uuid
warnings.filterwarnings('ignore')

import sys
sys.path.append('path/unionfind-master')
from unionfind import UnionFind # License: MIT. https://github.com/deehzee/unionfind

# =============================================================================
# Truck Trip Filtering
# For a particular truck trip, it will be considered if its trip ends fall into one mile from the centroid
# of a truck parking lot and/or fall into quarter mile from the centroid of a gas station.
# =============================================================================

def prepare_trip_ends(trip_roster_file, colname_file):
    # Prepare point file
    start_time = time.time()
    col_names = pd.read_csv(colname_file)
    trip_raw = pd.read_csv(trip_roster_file, header = None, names = col_names.columns)
    print(trip_raw.isna().sum())
    # trip_start
    trip_start = trip_raw[['TripId', 'StartLocLat', 'StartLocLon']]
    trip_start_gdf = gpd.GeoDataFrame(trip_start, crs={'init' :'epsg:4326'}, 
                                      geometry=[sg.Point(xy) for xy in zip(trip_start.StartLocLon, trip_start.StartLocLat)])
    trip_start_gdf = trip_start_gdf.to_crs({'init': 'epsg:6350'})
    trip_end = trip_raw[['TripId', 'EndLocLat', 'EndLocLon']]
    trip_end_gdf = gpd.GeoDataFrame(trip_end, crs={'init' :'epsg:4326'}, 
                                    geometry=[sg.Point(xy) for xy in zip(trip_end.EndLocLon, trip_end.EndLocLat)])
    trip_end_gdf = trip_end_gdf.to_crs({'init': 'epsg:6350'})
    print('Preparing trip ends takes %s secs.' % (time.time() - start_time))
    try:
        trip_short = trip_raw[['TripId', 'DeviceId','ProviderId', 'StartDate', 'EndDate',
                               'StartLocLat', 'StartLocLon', 'EndLocLat','EndLocLon',
                               'GeospatialType','VehicleWeightClass','EndpointType',
                               'TripMeanSpeedKph','TripMaxSpeedKph', 'TripDistanceMeters','MovementType']]
    except:
        trip_short = trip_raw[['TripId', 'DeviceId','ProviderId', 'StartDate', 'EndDate',
                               'StartLocLat', 'StartLocLon', 'EndLocLat','EndLocLon',
                               'GeospatialType','VehicleWeightClass']]
    return trip_start_gdf, trip_end_gdf, trip_short



def mp_filter(full_trip, park_shp_buffer):
    start_time = time.time()
    processes = 36
    chunksize = int(len(full_trip)/processes)
    print("Start parallel processing with chunksize %s" % chunksize)
    l = [full_trip[i:i+chunksize] for i in range(0, len(full_trip), chunksize)]
    
    def worker(park_shp_buffer, sub_trip_gdf):
        start_time = time.time()
        park_trip = gpd.sjoin(park_shp_buffer, sub_trip_gdf[['TripId','geometry']], how = 'inner', op = 'contains')
        print('%s trip takes %s sec' % (len(sub_trip_gdf), time.time()- start_time))
        gc.collect()
        return park_trip[['OBJECTID', 'TripId']]
    
    p = Pool(processes)
    partial_worker = partial(worker, park_shp_buffer)
    result = p.map(partial_worker, l)
    print("---Parallel (%s) processing the spatial join takes %s seconds for %s trips and %s stops---" 
              % (processes,(time.time() - start_time),len(full_trip),len(park_shp_buffer)))
    try:
        p.terminate()
    except:
        print('Pool failed to terminate itself!')
    gc.collect()
    return result

def filter_trip(trip_roster_file, parking_shp_nad83, gas_shp_nad83, 
                output_dir, colname_file):
    trip_start_gdf, trip_end_gdf, trip_short = prepare_trip_ends(trip_roster_file, colname_file)
    # prepare output shapefile
    trip_file_name = trip_roster_file.split("\\")[-1]
    park_trip_start_file = '/'.join([output_dir, trip_file_name.replace('.csv', '_start_park.csv')])
    park_trip_end_file = '/'.join([output_dir, trip_file_name.replace('.csv', '_end_park.csv')])
    gas_trip_start_file = '/'.join([output_dir, trip_file_name.replace('.csv', '_start_gas.csv')])
    gas_trip_end_file = '/'.join([output_dir, trip_file_name.replace('.csv', '_end_gas.csv')])
    
    # Truck parking stops
    def park_buffer(point_geo):
        return point_geo.buffer(1609)
    
    park_shp_buffer = gpd.GeoDataFrame(parking_shp_nad83.geometry.apply(park_buffer), 
                                       crs = {'init': 'epsg:6350'})
    park_shp_buffer['OBJECTID'] = parking_shp_nad83['OBJECTID']
    
    def flatten_list_df(list_df):
        flatten_df = pd.DataFrame()
        for df in list_df:
            flatten_df = flatten_df.append(df, ignore_index = True)
        return flatten_df
    
    park_trip_start = mp_filter(trip_start_gdf, park_shp_buffer)
    park_trip_start = flatten_list_df(park_trip_start)
    park_trip_start = park_trip_start.drop_duplicates(subset = 'TripId', keep = 'first')
    park_trip_end = mp_filter(trip_end_gdf, park_shp_buffer)
    park_trip_end = flatten_list_df(park_trip_end)
    park_trip_end = park_trip_end.drop_duplicates(subset = 'TripId', keep = 'first')
    
    # Attach additional info of trip
    park_trip_start = park_trip_start.merge(trip_short, how = 'left', sort = False)
    park_trip_end = park_trip_end.merge(trip_short, how = 'left', sort = False)
    park_trip_start.to_csv(park_trip_start_file, index = False)
    park_trip_end.to_csv(park_trip_end_file, index = False)
    
    # Gas stations
    gas_shp_nad83.loc[:, 'OBJECTID'] = list(range(1, len(gas_shp_nad83)+1, 1))
    def gas_buffer(point_geo):
        return point_geo.buffer(402)
    
    gas_shp_buffer = gpd.GeoDataFrame(gas_shp_nad83.geometry.apply(gas_buffer), 
                                      crs = {'init': 'epsg:6350'})
    gas_shp_buffer['OBJECTID'] = gas_shp_nad83['OBJECTID']
    gas_trip_start = mp_filter(trip_start_gdf, gas_shp_buffer)
    gas_trip_start = flatten_list_df(gas_trip_start)
    gas_trip_start = gas_trip_start.drop_duplicates(subset = 'TripId', keep = 'first')
    gas_trip_end = mp_filter(trip_end_gdf, gas_shp_buffer)
    gas_trip_end = flatten_list_df(gas_trip_end)
    gas_trip_end = gas_trip_end.drop_duplicates(subset = 'TripId', keep = 'first')
    # Attach additional info of trip
    gas_trip_start = gas_trip_start.merge(trip_short, how = 'left', sort = False)
    gas_trip_end = gas_trip_end.merge(trip_short, how = 'left', sort = False)
    gas_trip_start.to_csv(gas_trip_start_file, index = False)
    gas_trip_end.to_csv(gas_trip_end_file, index = False)

# =============================================================================
# Truck Trip Matching via FIFO with minimum and maximum dwell time at gas stations.
# For a particular gas station, the first trip finished will be the first trip left after 1 min 
# but within 180 min.
# =============================================================================

def mp_match(overlap_stop, input_list):
    overlap_stop = overlap_stop.tolist()
    start_time = time.time()
    processes = 26
    chunksize = int(len(overlap_stop)/processes)
    print("Start parallel processing with chunksize %s" % chunksize)
    l = [overlap_stop[i:i+chunksize] for i in range(0, len(overlap_stop), chunksize)]
    
    def worker(lst_df, lst_stop):
        start_clean, end_clean, min_stay, max_stay = lst_df
    
        def single_stop(end_trip, start_trip, min_stay = 1, max_stay= 180):
            end_trip_sorted = end_trip.sort_values(by = ['DeviceId', 'EndDate'])
            start_trip_sorted = start_trip.sort_values(by = ['DeviceId', 'StartDate'])
            match_pair = pd.DataFrame()
            for index, row in end_trip_sorted.iterrows():
                if len(match_pair) > 0:
                    candidate = start_trip_sorted.loc[(~start_trip_sorted['TripId'].isin(match_pair['StartTripId']))
                    &(start_trip_sorted['StartDate'] > row['EndDate'])]
                else:
                    candidate = start_trip_sorted.loc[start_trip_sorted['StartDate'] > row['EndDate']]
                candidate = candidate.loc[(candidate['ProviderId'] == row['ProviderId'])&(candidate['VehicleWeightClass'] == row['VehicleWeightClass'])]
                # First pair with DeviceId
                device_match = candidate.loc[candidate['DeviceId'] == row['DeviceId']]
                if len(device_match)>0:
                    device_match['sec_delta'] = device_match['StartDate'] - row['EndDate']
                    if device_match['sec_delta'].min().seconds < 2*3600:
                        tmp_pair = [row['TripId'],device_match.loc[device_match['sec_delta'].idxmin(), 'TripId']]
                        match_pair = match_pair.append(pd.DataFrame([tmp_pair], columns = ['EndTripId', 'StartTripId']))
                else:
                    candidate['sec_delta'] = candidate['StartDate'] - row['EndDate']
                    narrow_candidate = candidate.loc[(candidate['sec_delta'] >= timedelta(minutes = min_stay))
                    &(candidate['sec_delta'] < timedelta(minutes = max_stay))]
                    if len(narrow_candidate) > 0:
                        narrow_candidate = narrow_candidate.iloc[0]
                        tmp_pair = [row['TripId'], narrow_candidate['TripId']]
                        match_pair = match_pair.append(pd.DataFrame([tmp_pair], columns = ['EndTripId', 'StartTripId']))
            match_pair['StopId'] = end_trip['OBJECTID'].drop_duplicates().iloc[0]
            return match_pair
        
        match_pair = pd.DataFrame()
        for stop in lst_stop:
            end_trip = end_clean.loc[end_clean['OBJECTID'] == stop]
            start_trip = start_clean.loc[start_clean['OBJECTID']==stop]
            tmp_match_pair = single_stop(end_trip, start_trip, min_stay, max_stay)
            match_pair = match_pair.append(tmp_match_pair, sort = False)
            gc.collect()
        return match_pair
    
    def flatten_list_df(list_df):
        flatten_df = pd.DataFrame()
        for df in list_df:
            flatten_df = flatten_df.append(df, ignore_index = True, sort = False)
        return flatten_df
    
    p = Pool(processes)
    partial_worker = partial(worker, input_list)
    result = p.map(partial_worker, l)
    flatten_match = flatten_list_df(result)
    print("---Parallel (%s) finds %s matches for %s trips and %s stops in %s seconds---" 
              % (processes,len(flatten_match),len(input_list[1]),len(overlap_stop),(time.time() - start_time)))
    p.terminate()
    gc.collect()
    return flatten_match

def fifo(stop_type = 'gas', min_stay = 1, max_stay= 180):
    start_file = 'trip_start_' + stop_type +'.csv'
    end_file = 'trip_end_' + stop_type +'.csv'
    start_here = pd.read_csv(start_file)
    end_here = pd.read_csv(end_file)
    
    def clean_df(df, col = 'StartDate'):
        df.loc[:, col] = df[col].apply(lambda timestring: dateutil.parser.parse(timestring))
        df_sorted = df.sort_values(by = ['OBJECTID',col])
        return df_sorted
    
    start_here_sorted = clean_df(start_here)
    end_here_sorted = clean_df(end_here, 'EndDate')
    # Only keep gas stations with in and out volumes
    unique_stop = start_here_sorted['OBJECTID'].drop_duplicates()
    overlap_stop = unique_stop.loc[unique_stop.isin(end_here_sorted['OBJECTID'].drop_duplicates())]
    start_clean = start_here_sorted.loc[start_here_sorted['OBJECTID'].isin(overlap_stop)].reset_index(drop = True)
    end_clean = end_here_sorted.loc[end_here_sorted['OBJECTID'].isin(overlap_stop)].reset_index(drop = True)
    # Multi-processing
    input_list = [start_clean, end_clean, min_stay, max_stay]
    all_match = mp_match(overlap_stop, input_list)
    all_match = all_match.drop_duplicates('StartTripId')
    all_match = all_match.drop_duplicates('EndTripId')
    return all_match

# =============================================================================
# Truck Trip Matching via FIFO with resting time constraint at truck parking stations.
# For a particular truck stop parking, the first trip finished will be the first trip left after 10 hours.
# =============================================================================
def pre_process(trip_start_park, trip_end_park):
    
    lambda_date = lambda timestring: dateutil.parser.parse(timestring)
    trip_end_park['s_dt'] = list(map(lambda_date,trip_end_park['StartDate']))
    trip_end_park['e_dt'] = list(map(lambda_date,trip_end_park['EndDate']))
    trip_start_park['s_dt'] = list(map(lambda_date,trip_start_park['StartDate']))
    trip_start_park['e_dt'] = list(map(lambda_date,trip_start_park['EndDate']))
    
    lambda_time = lambda date:time.mktime(date.timetuple())
    trip_end_park['s_ts'] = list(map(lambda_time,trip_end_park['s_dt']))
    trip_end_park['e_ts'] = list(map(lambda_time,trip_end_park['e_dt']))
    trip_start_park['s_ts'] = list(map(lambda_time,trip_start_park['s_dt']))
    trip_start_park['e_ts'] = list(map(lambda_time,trip_start_park['e_dt']))
    
    trip_end_park['tt'] =  trip_end_park['e_ts'] - trip_end_park['s_ts']
    trip_start_park['tt'] =  trip_start_park['e_ts'] - trip_start_park['s_ts']
    
    trip_end_park = trip_end_park.sort_values(by=['e_ts']).reset_index().drop(['index'],axis=1)
    trip_start_park = trip_start_park.sort_values(by=['s_ts']).reset_index().drop(['index'],axis=1)
    return trip_start_park,trip_end_park 

def truck_trip_matching_main(trip_start_park, trip_end_park):
    park_lst = trip_end_park['OBJECTID'].unique().tolist()
    
    def Provider_ID_match(ProviderId,candidates):
        provider_candidates = candidates[candidates['ProviderId']==ProviderId]
        return provider_candidates
    def Device_ID_match(DeviceId,candidates):
        device_candidates = candidates[candidates['DeviceId']==DeviceId]
        return device_candidates
    def Weight_Class_match(weight_class,candidates):
        weight_candidates = candidates[candidates['VehicleWeightClass']==weight_class]
        return weight_candidates
    def Resting_Time_match(tt, arrival_ts, candidates):
        # maximum driving 11 hours, have to rest 10 hours
        if tt >= 11*3600:
            resting_time = 10*3600
        else:
            resting_time = 300
        estimate_departure = arrival_ts + resting_time
        upper_bound = arrival_ts + 24*3600
        time_candidates = candidates[(candidates['s_ts']>=estimate_departure) & (candidates['s_ts']<=upper_bound)]
        return time_candidates
    
    Pair_Results = []
    Trip_Rosters = []
    
    for park_id in park_lst:
        end_trip = trip_end_park[trip_end_park['OBJECTID']==park_id]
        start_trip = trip_start_park[trip_start_park['OBJECTID']==park_id]
        for i in range(len(end_trip)):
                        
            one_trip_info = end_trip.iloc[i]
            tt = one_trip_info.tt
            arrival_ts = one_trip_info.e_ts
            weight_class = one_trip_info.VehicleWeightClass
            ProviderId = one_trip_info.ProviderId
            DeviceId = one_trip_info.DeviceId
            TripId = one_trip_info.TripId
            
            # Hard Constraints
            provider_candidates = Provider_ID_match(ProviderId,start_trip)
            
            device_candidates = Device_ID_match(DeviceId,provider_candidates)
            
            weight_candidates = Weight_Class_match(weight_class,device_candidates)
            
            time_candidates = Resting_Time_match(tt,arrival_ts,weight_candidates)
            
            if time_candidates.empty:
                time_candidates = Resting_Time_match(tt,arrival_ts,start_trip)
                
                weight_candidates = Weight_Class_match(weight_class,time_candidates)
                
                provider_candidates = Provider_ID_match(ProviderId,weight_candidates)
                
                device_candidates = Device_ID_match(DeviceId,provider_candidates)
                
                if device_candidates.empty:
                    final_candidates = provider_candidates
                else:
                    final_candidates = device_candidates
            else:
                final_candidates = time_candidates
            
            if final_candidates.empty:
                continue

            paired_trip = final_candidates.iloc[0]
            paired_trip_id = paired_trip.TripId
            
            new_trip_id = str(uuid.uuid4())
            new_trip_provider = ProviderId
            new_trip_olat = one_trip_info.StartLocLat
            new_trip_olon = one_trip_info.StartLocLon
            new_trip_dlat = paired_trip.EndLocLat
            new_trip_dlon = paired_trip.EndLocLon
            new_trip_start_time = one_trip_info.StartDate
            new_trip_end_time = paired_trip.EndDate
            new_trip_midlat = one_trip_info.EndLocLat
            new_trip_midlon = one_trip_info.EndLocLon
            new_trip_mid_start_time = one_trip_info.EndDate
            new_trip_mid_end_time = paired_trip.StartDate
            tt_first = tt
            tt_second = paired_trip.tt
            device_first = DeviceId
            device_second = paired_trip.DeviceId
            new_trip_stay_time = paired_trip.s_ts - one_trip_info.e_ts
             
            Pair_Results.append([new_trip_id,park_id,TripId,paired_trip_id])
            Trip_Rosters.append([new_trip_id,new_trip_provider,
                                 new_trip_olat,new_trip_olon,new_trip_dlat,new_trip_dlon,
                                 new_trip_start_time,new_trip_end_time,
                                 new_trip_midlat,new_trip_midlon,new_trip_mid_start_time,new_trip_mid_end_time,
                                 tt_first,tt_second,new_trip_stay_time,
                                 device_first,device_second])
            
            start_trip = start_trip[start_trip['TripId']!=paired_trip_id]
    
    Pair_Results = pd.DataFrame(Pair_Results,columns=['TripId','OBJECTID','end','start'])
    Trip_Rosters = pd.DataFrame(Trip_Rosters,columns=['TripId','ProviderId',
                                                      'StartLocLat','StartLocLon','EndLocLat','EndLocLon',
                                                      'StartDate','EndDate',
                                                      'MidLocLat','MidLocLon','MidStartDate','MidEndDate',
                                                      'travel_time_1','travel_time_2','stay_time',
                                                      'DeviceId_1','DeviceId_2'])
    Pair_Results.to_csv('Pair_Results.csv',index=False)
    Trip_Rosters.to_csv('Trip_Rosters_New.csv',index=False)
    return Pair_Results,Trip_Rosters

def trip_roster_merged(trip_roster_file, colname_file, trip_chain,
                       park_pair_file, gas_pair_file):
    col_names = pd.read_csv(colname_file)
    trip_roster = pd.read_csv(trip_roster_file, header = None, names = col_names.columns)
    if trip_chain == False:
        matched_trip_pair = pd.read_csv(park_pair_file)
        matched_trip_pair.columns = ['TripId', 'StopId', 'EndTripId', 'StartTripId']
        matched_trip_pair = matched_trip_pair[['EndTripId', 'StartTripId', 'StopId']].append(pd.read_csv(gas_pair_file))
        trip_pair_id = ['EndTripId', 'StopId', 'StartTripId']
        new_pair_df = matched_trip_pair
    else:
        matched_trip_pair = pd.read_csv(park_pair_file, usecols = ['end', 'start'])
        matched_trip_pair.columns = ['EndTripId', 'StartTripId']
        matched_trip_pair = matched_trip_pair.append(pd.read_csv(gas_pair_file, usecols = ['EndTripId', 'StartTripId']))
        trip_pair_id = ['EndTripId', 'StartTripId']
        # trip chaining
        start_time = time.time()
        uf = UnionFind(list(set(matched_trip_pair.values.flatten())))
        for index, row in matched_trip_pair.iterrows():
            uf.union(row['EndTripId'], row['StartTripId'])
        result = uf.components()
        print('Trip chaining takes %s secs for %s trip pairs.' % (time.time() - start_time, len(matched_trip_pair)))
        
        def set_first_last(input_set):
            tmp = list(input_set)
            return [tmp[0], tmp[-1]]
        
        new_pair = map(set_first_last, result)
        new_pair_df = pd.DataFrame(new_pair, columns = ['EndTripId', 'StartTripId'])
    trip_unmatched = trip_roster.loc[~trip_roster['TripId'].isin(
        matched_trip_pair[['EndTripId', 'StartTripId']].values.flatten())]
    # create od file for matched trips/ trip chain
    trip_od = new_pair_df.merge(trip_roster[['TripId', 'StartLocLat', 'StartLocLon']].rename(
        columns = {'TripId': 'EndTripId'}), how = 'left', sort = False)
    trip_od = trip_od.merge(trip_roster[['TripId', 'EndLocLat', 'EndLocLon']].rename(
        columns = {'TripId': 'StartTripId'}), how = 'left', sort = False)
    trip_od['TripId'] = trip_od[trip_pair_id].apply(lambda row: '_'.join(row.tolist()), axis = 1)
    return trip_unmatched, trip_od

if __name__ == '__main__':
    # Filter trips
    filter_trip(trip_roster_file, parking_shp_nad83, gas_shp_nad83, output_dir,
                colname_file)
    
    # Pair Trips ending in gas stations
    all_match = fifo()
    all_match.to_csv('gas_Pair_Results.csv', index = False)
    
    # Pair Trips ending in truck parking stations
    trip_start_park = pd.read_csv('trip_start_park.csv')
    trip_end_park = pd.read_csv('trip_end_park.csv')
    
    trip_start_park, trip_end_park = pre_process(trip_start_park,trip_end_park)
    
    Pair_Results, Trip_Rosters= truck_trip_matching_main(trip_start_park,trip_end_park)
    Pair_Results.to_csv('park_Pair_Results.csv',index=False)
    
    # Chaining trips
    trip_unmatched, trip_od = trip_roster_merged(trip_roster_file, colname_file, True,
                                                 'park_Pair_Results.csv',
                                                 'gas_Pair_Results.csv')
    trip_unmatched.to_csv('Unmatched_Trips.csv', index = False)
    trip_od.to_csv('Chained_Trip_OD.csv', index = False)
    
    
