import fiona
import shapely.geometry as sg
from shapely.geometry import asMultiLineString
import time
import pandas as pd
import geopandas
import csv
import os
import numpy as np


def calc_dist(in_pts_file, in_net_file, out_path, out_name):
    
    # Create shapefile from csv file
    start_time = time.time()
    shapeout = in_pts_file.replace(".csv","_fiona.shp")
    yourschema =  {'geometry': 'Point',
                   'properties': {'point_id': 'int','long': 'float','lat': 'float'}}
    try:
        with fiona.open(shapeout, 'w',crs=fiona.crs.from_epsg(4326),driver='ESRI Shapefile', schema=yourschema) as output:
            reader = pd.read_csv(in_pts_file)
            count = 0
            for index, row in reader.iterrows():
                # geometry       
                tmp_point = sg.Point(float(row['long']), float(row['lat']))
                # attributes
                prop = {'point_id': int(count),'long': float(row['long']),'lat': float(row['lat'])}
                # write the row (geometry + attributes in GeoJSON format)
                output.write({'geometry': sg.mapping(tmp_point), 'properties':prop})
                count += 1
            del row, reader
            output.close()
            print(output.closed)
    except:
        print(output.closed)
    print("---Creating a shapefile takes %s seconds for %s points---" % ((time.time() - start_time),count))
    
    # Convert crs
    start_time = time.time()
    points = geopandas.read_file(shapeout)
    # change CRS to epsg 6350: NAD83
    points = points.to_crs({'init': 'epsg:6350'})
    points.to_file(shapeout.replace(".shp", "_NAD83.shp"))
    print("---Projecting a shapefile takes %s seconds for %sMB file---" % ((time.time() - start_time),
                                                                          int(os.path.getsize(shapeout)*1e-6)))
    
    # Read network file
    start_time = time.time()
    network = geopandas.read_file(in_net_file)
    # change CRS to epsg 6350: NAD83
    network = network.to_crs({'init': 'epsg:6350'})
    network.to_file(in_net_file.replace(".shp", "_NAD83.shp"))
    print("---Projecting a shapefile takes %s seconds for %sMB file---" % ((time.time() - start_time),
                                                                          int(os.path.getsize(in_net_file)*1e-6)))
    net_lines = network['geometry']
    line_output = asMultiLineString(net_lines)
    net_shply = sg.MultiLineString(line_output)
    print("Network converted to Shapely geometry object.")
    
    # Calculate distance
    start_time = time.time()
    distance = []
    count = 0
    with fiona.open(shapeout.replace(".shp", "_NAD83.shp")) as coords:
        for feature in coords:
            geom = sg.shape(feature["geometry"])
            distance.append( geom.distance(net_shply) )
            count += 1
    print("---Distance calculation takes %s seconds for %s points and %sMB network file---" 
          % ((time.time() - start_time),count,int(os.path.getsize(in_net_file)*1e-6)))
    
    # Store distance to csv
    with open(out_path + '\\' + out_name, 'w', newline='') as myfile: 
        w = csv.writer(myfile, delimiter=',')
        w.writerows(zip(distance)) 
    myfile.close()
    print("Distance calculated and stored!")

calc_dist(in_pts_file = 'Point_Directory',
         in_net_file ='Network_Directory' ,
         out_path = 'Output_Directory',
         out_name = 'Output_Name')

