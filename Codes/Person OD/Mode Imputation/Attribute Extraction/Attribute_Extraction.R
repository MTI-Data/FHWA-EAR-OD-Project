#Install the following libraries if they are not already installed (Use install.packages())
library(data.table)
library(plyr)
library(geosphere)
library(foreach)
library(doSNOW)
library(chron)

#function to calculate timedifference
timedifference=function(d,o){
  a=o
  b=d
  date_a=unlist(strsplit(a,"T",fixed=T))[1]
  date_b=unlist(strsplit(b,"T",fixed=T))[1]
  temp_a=unlist(strsplit(a,"T",fixed=T))[2]
  temp_b=unlist(strsplit(b,"T",fixed=T))[2]
  time_a=unlist(strsplit(temp_a,"Z",fixed=T))[1]
  time_b=unlist(strsplit(temp_b,"Z",fixed=T))[1]
  chron_a=chron(dates=date_a,times=time_a,format=c(dates="y-m-d",times="h:m:s"))
  chron_b=chron(dates=date_b,times=time_b,format=c(dates="y-m-d",times="h:m:s"))
  as.numeric(difftime(chron_b,chron_a,unit="secs"))
}
timedifference_vec=Vectorize(timedifference,c("o","d"))

#function to calculate distance, time, and speed attribute for each trip
calc_attr = function(trip_index,file=point_file) {
  
  # temporary point array for for each trip
  tmp_df = file[trip_index[1]:trip_index[2],]
  tmp_df = as.data.frame(tmp_df)
  npoint = nrow(tmp_df)
  # OD distance and time
  od_dist = distGeo(tmp_df[1,c("longitude","latitude")],tmp_df[npoint,c("longitude","latitude")])
  od_time = timedifference(tmp_df[npoint,"start_time"],tmp_df[1,"end_time"])/60
  # Average number of record per hour
  avg_record = round(npoint/od_time*60,2)
  # Average speed of trip meters/sec
  avg_spd = round(od_dist/od_time/60,3)
  #point to point speed distribution
  time_vec = timedifference_vec(tmp_df[2:npoint,"start_time"],tmp_df[1:(npoint-1),"end_time"])
  nonzero_id = which(time_vec>0)
  if(length(nonzero_id)>=3) {
    time_vec = time_vec[nonzero_id]
    dist_vec = distGeo(tmp_df[1:(npoint-1),c("longitude","latitude")],tmp_df[2:npoint,c("longitude","latitude")])
    dist = sum(dist_vec)
    dist_vec= dist_vec[nonzero_id]
    spd_vec = round(dist_vec/time_vec,3)
    max_spd = max(spd_vec,na.rm = T)
    min_spd = min(spd_vec,na.rm = T)
    med_spd = median(spd_vec,na.rm = T)
    Q_spd = quantile(spd_vec,c(.05,.75,.95),na.rm = T)
  } else {
    dist = od_dist
    max_spd = avg_spd
    min_spd = avg_spd
    med_spd = avg_spd
    Q_spd = rep(avg_spd,3)
  }
  #cretae a vetor of attributes and return it as the output
  attr_vec = c(npoint,od_dist, dist, od_time, avg_record, avg_spd, max_spd, min_spd, med_spd, Q_spd)
  attr_vec
}

#function to add distance to network attributes to trip files
append_dist_attr = function(trip_file, point_file,mode) {
  #seperated by mode, calculate average distance, max and min distances, and quantiles of distance
  if(mode=="Rail") {
    dist_sum = ddply(point_file,.(trip_id),summarise,Rail_Avg = mean(rail_dist),
                     Rail_Max = max(rail_dist),Rail_Min = min(rail_dist),
                     Rail_Med = median(rail_dist), Rail_Q5=quantile(rail_dist, 0.05),
                     Rail_Q75 = quantile(rail_dist, 0.75),
                     Rail_Q95 = quantile(rail_dist, 0.95))
  } else if (mode=="Bus") {
    dist_sum = ddply(point_file,.(trip_id),summarise,Bus_Avg = mean(bus_dist),
                     Bus_Max = max(bus_dist),Bus_Min = min(bus_dist),
                     Bus_Med = median(bus_dist), Bus_Q5=quantile(bus_dist, 0.05),
                     Bus_Q75 = quantile(bus_dist, 0.75),
                     Bus_Q95 = quantile(bus_dist, 0.95))
  } else if (mode=="Road") {
    dist_sum = ddply(point_file,.(trip_id),summarise,Road_Avg = mean(road_dist),
                     Road_Max = max(road_dist),Road_Min = min(road_dist),
                     Road_Med = median(road_dist), Road_Q5=quantile(road_dist, 0.05),
                     Road_Q75 = quantile(road_dist, 0.75),
                     Road_Q95 = quantile(road_dist, 0.95))
  } else if (mode=="Air") {
    dist_sum = ddply(point_file,.(trip_id),summarise,Air_Avg = mean(air_dist),
                     Air_Max = max(air_dist),Air_Min = min(air_dist),
                     Air_Med = median(air_dist), Air_Q5=quantile(air_dist, 0.05),
                     Air_Q75 = quantile(air_dist, 0.75),
                     Air_Q95 = quantile(air_dist, 0.95))
  } else {
    return(paste(mode,"doesn't exist!",sep=" "))
  }
  #set column names
  colnames(dist_sum)<-c("trip_id",paste(mode,c("Avg","Max","Min",
                                               "Med","Q5","Q75","Q95"),sep="_"))
  #merge the results with the current trp dataset and return the merged dataset
  new_attr = merge(trip_file,dist_sum,by="trip_id",sort=F,all.x=T,all.y=F)
  new_attr
}

attr_extr=function(core){
    #Read the trip file
    point_file = fread(paste0("Point_File_Name_",core,".csv"),stringsAsFactors = F)

    # prepare trip start/end index array
    trip_start_id=which(!duplicated(point_file$trip_id))
    trip_end_id=which(!duplicated(point_file$trip_id,fromLast=TRUE))
    trip_ind_array = cbind(trip_start_id,trip_end_id)
    dim(trip_ind_array)

    # calculate the distance,time, and speed attributes
    trip_file = as.data.frame(t(apply(trip_ind_array, 1, calc_attr,file = point_file)))
    colnames(trip_file)<-c("No_point","OD_dist", "trip_dist","trip_time","Avg_record", "speed_avg", 
                               "speed_max","speed_min","speed_med",
                               "speed_Q5","speed_Q75","speed_Q95")
    trip_file$trip_id = unique(point_file$trip_id)


    ## Add trip distance-to-network attributes to trip file
      trip_file = append_dist_attr(trip_file, point_file, mode="Rail")
      trip_file = append_dist_attr(trip_file, point_file, mode="Bus")
      trip_file = append_dist_attr(trip_file, point_file, mode="Road")
      trip_file = append_dist_attr(trip_file, point_file, mode="Air")

    # Export full attribute table
    write.csv(trip_file,paste0("Trip_File_Name",core,".csv"),row.names = F)
}

#paralel application
cl = makeCluster(32)
registerDoSNOW(cl)
getDoParWorkers()
clusterExport(cl,c("timedifference","calc_attr","append_dist_attr"),envir=.GlobalEnv)
foreach(core=0:31,.packages=c("chron","data.table","plyr","geosphere")) %dopar%  attr_extr(core)
stopCluster(cl)
rm(cl)
closeAllConnections()

