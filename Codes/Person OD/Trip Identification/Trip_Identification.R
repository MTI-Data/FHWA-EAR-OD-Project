library(hashids)
library(data.table)
library(geosphere)
library(foreach)
library(doSNOW)

#trip identification function
trip_identification_function=function(core,ind,device_ind_array){
  
  #set the parameters
  distance_threshold=300 #meters
  time_threshold=300 #seconds
  speed_threshold= 1.4 #meters/second
  
  first_row=floor(nrow(device_ind_array)*core/20)+1
  last_row=floor(nrow(device_ind_array)*(core+1)/20)
  rows_to_skip=as.numeric(device_ind_array[first_row,"device_start_id"])-1
  nrows_to_read=as.numeric(device_ind_array[last_row,"device_end_id"]-device_ind_array[first_row,"device_start_id"])+1
  #read column names
  col_names=names(fread(paste0("../user_point_file_",ind,"_cluster.csv"),nrow=1))
  #read the data
  if (core==0){
    gps_data=fread(paste0("../user_point_file_",ind,"_cluster.csv"),skip=rows_to_skip,nrow=nrows_to_read,stringsAsFactors = F,header = T)
    
  }else{
    gps_data=fread(paste0("../user_point_file_",ind,"_cluster.csv"),skip=rows_to_skip,nrow=nrows_to_read,stringsAsFactors = F,header = F)
    names(gps_data)<-col_names
  }
  gps_data=gps_data[order(userid,utc_timestamp)]
  

  #removs static points
  gps_data$prev_cl=c(0,gps_data$final_clust[-1])
  gps_data$next_cl=c(gps_data$final_clust[-nrow(gps_data)],0)
  gps_data=gps_data[!((final_clust==prev_cl)&(final_clust==next_cl)&(clust_flag=="static")),1:9]

  #remove devices with insufficient points
  obs_by_dev=gps_data[,.(.N),userid]
  dev_to_remove=obs_by_dev[N<2,userid]
  gps_data=gps_data[!userid%in%dev_to_remove,]
  
  
  #calculate indices for the observations of each device
  device_start_id=which(!duplicated(gps_data$userid))
  device_end_id=which(!duplicated(gps_data$userid,fromLast=TRUE))
  device_ind_array = cbind(device_start_id,device_end_id)
  
  #add_time_from
  gps_data$time_from=mat.or.vec(nrow(gps_data),1)
  add_time=function(i){
    tmp_df = gps_data[device_ind_array[i,1]:device_ind_array[i,2],]
    tmp_df = as.data.frame(tmp_df)
    npoint = nrow(tmp_df)
    time_vec = c(as.numeric(tmp_df[2:npoint,"utc_timestamp"])-as.numeric(tmp_df[1:(npoint-1),"utc_timestamp"]),0)
    gps_data[device_ind_array[i,1]:device_ind_array[i,2],time_from:=time_vec]
  }
  add_time_vec=Vectorize(add_time)
  add_time_vec(1:nrow(device_ind_array))
  
  #add_distance_from
  gps_data$dist_from=mat.or.vec(nrow(gps_data),1)
  add_dist=function(i){
    tmp_df = gps_data[device_ind_array[i,1]:device_ind_array[i,2],c("lat","long")]
    npoint = nrow(tmp_df)
    dist_vec = c(distGeo(tmp_df[1:(nrow(tmp_df)-1),c("long","lat")],tmp_df[2:nrow(tmp_df),c("long","lat")]),0)
    gps_data[device_ind_array[i,1]:device_ind_array[i,2],dist_from:=dist_vec]
  }
  add_dist_vec=Vectorize(add_dist)
  add_dist_vec(1:nrow(device_ind_array))
  
  #add speed from
  gps_data$speed_from=mat.or.vec(nrow(gps_data),1)
  gps_data[time_from!=0,speed_from:=gps_data[time_from!=0,dist_from]/gps_data[time_from!=0,time_from]]
  
  #function to identify trips
  trip_identifier=function(row,datatable){
    
    hash_id=datatable$trip_id[row-1]

    if (datatable[row-1,trip_id]=="0"){
      if(datatable[row,speed_from]<speed_threshold){
        datatable[row,trip_id:="0"]
      }else{
        datatable[row,trip_id:=hash_generator()]
      }
    }else{
      if (datatable[row-1,speed_from]>speed_threshold){
        datatable[row,trip_id:=hash_id]
      }else{
        if (datatable[row-1,dist_from]>distance_threshold){
          if (datatable[row,speed_from]<speed_threshold){
            datatable[row,trip_id:="0"]
          }else{
            datatable[row,trip_id:=hash_generator()]
          }
        }else{
          if (datatable[row-1,time_from]>time_threshold){
            if (datatable[row,speed_from]<speed_threshold){
              datatable[row,trip_id:="0"]
            }else{
              datatable[row,trip_id:=hash_generator()]
            }
          }else{
            datatable[row,trip_id:=hash_id]
            datatable[row,time_from:=(datatable[row,time_from]+datatable[row-1,time_from])]
          }
        }
      }
    }
  }
  
  flagged_trip_identifier=function(row,datatable){

    previous_speed_from=datatable$speed_from[row-1]
    previous_cluster_type=datatable$clust_flag[row-1]
    current_cluster=datatable$final_clust[row]
    current_cluster_type=datatable$clust_flag[row]
    
    if ((current_cluster_type==previous_cluster_type)&(current_cluster_type=="static")&(previous_speed_from<speed_threshold)){
      flagged_id=hash_generator()
      datatable[row-1,flagged_trip_id:=flagged_id]
      datatable[row,flagged_trip_id:=flagged_id]
    }
  }
  
  #set all trip IDs
  gps_data$trip_id=as.character(mat.or.vec(nrow(gps_data),1))
  gps_data$flagged_trip_id=as.character(mat.or.vec(nrow(gps_data),1))
  
  set_trip=function(device){
    first_row= device_ind_array[device,1]
    last_row=device_ind_array[device,2]
    if (gps_data[first_row,speed_from]>speed_threshold){
      gps_data[first_row,trip_id:=hash_generator()]
    }else{
      gps_data[first_row,trip_id:="0"]
    }
    for (obs in (first_row+1):last_row){
      trip_identifier(obs,gps_data)
    }
  }
  
  set_flagged_trip=function(device){
    first_row= device_ind_array[device,1]
    last_row=device_ind_array[device,2]
    for (obs in (first_row+1):last_row){
      flagged_trip_identifier(obs,gps_data)
    }
  }
  
  lapply(1:nrow(device_ind_array),set_trip)  
  lapply(1:nrow(device_ind_array),set_flagged_trip)

  #add device id to trip id
  gps_data$trip_id=paste(gps_data$device_id,gps_data$trip_id,sep='')
  
  #create trip files
  trips=gps_data[trip_id!="0",.(userid,utc_timestamp,tz_offset,lat,long,flag,accuracy,final_clust,clust_flag),trip_id]

  
  ##post processing to remove very short trips
  #calculate indices for the observations of each trip
  trip_start_index=which(!duplicated(trips$trip_id))
  trip_end_index=which(!duplicated(trips$trip_id,fromLast=TRUE))
  trip_ind_array = cbind(trip_start_index,trip_end_index)
  #function to calculate if all trip points belong to the same static cluster
  static_trips=c()
  check_static_trips=function(row){
    if((trip_start_index[row]+1==trip_end_index[row])&(trips[trip_start_index[row],clust_flag]=="static")){
      static_trips<<-c(static_trips,trips[trip_start_index[row],trip_id])
    }
  }
  check_static_trips_vec=Vectorize(check_static_trips)
  check_static_trips_vec(1:nrow(trip_ind_array))
  trips=trips[!(trip_id%in%static_trips)]
  
  #calculate indices for the observations of each trip
  trip_start_index=which(!duplicated(trips$trip_id))
  trip_end_index=which(!duplicated(trips$trip_id,fromLast=TRUE))
  trip_ind_array = cbind(trip_start_index,trip_end_index)
  
  #function to calculate trip time
  calculate_trip_time=function(i){
    start_time=as.numeric(trips$utc_timestamp[trip_start_index[i]])
    end_time=as.numeric(trips$utc_timestamp[trip_end_index[i]])
    trip_time=end_time-start_time
    return(trip_time)
  }
  #vectorize the function
  calculate_trip_time_vec=Vectorize(calculate_trip_time)
  #calculate vector of all trip times
  trip_time_vector=calculate_trip_time_vec(1:nrow(trip_ind_array))
  
  #function to calculate trip distance
  calculate_trip_distance=function(i){
    start_point=trips[trip_start_index[i],c("long","lat")]
    end_point=trips[trip_end_index[i],c("long","lat")]
    trip_distance=distGeo(start_point,end_point)
    return(trip_distance)
  }
  #vectorize the function
  calculate_trip_distance_vec=Vectorize(calculate_trip_distance)
  #calculate vector of all trip times
  trip_distance_vector=calculate_trip_distance_vec(1:nrow(trip_ind_array))
  #produce trip summary
  trip_summary=data.table("trip_id"=unique(trips$trip_id),"time"=trip_time_vector,"distance"=trip_distance_vector)
  #identify short trips
  short_trips=trip_summary[(distance<300),trip_id]
  #remove short trips
  trips=trips[!(trip_id%in%short_trips),]
  
  #summarize flagged_trips and bind to the trips
  flagged_trips=gps_data[flagged_trip_id!="0",.(userid,utc_timestamp,tz_offset,lat,long,flag,accuracy,final_clust,clust_flag),flagged_trip_id]
  all_trips=rbindlist(list(trips,flagged_trips))
  
  #remove files
  rm(gps_data)
  rm(trips)
  rm(flagged_trips)
  
  #return all trips
  return(all_trips)
}

#function to generate hashid for trip ID
h = hashid_settings(salt = 'this is my salt', min_length = 10)
hash_generator=function(){
  random=runif(1,100000000,900000000)
  encode(as.integer(random),h)
}

#Run trip identification function in parallel
indices=sprintf("%02d",c(0:27))
for (ind in indices){
  #read the data
  gps_data=fread(paste0("../user_point_file_",ind,"_cluster.csv"),stringsAsFactors = F)
  gps_data=gps_data[order(userid,utc_timestamp)]
  
  #calculate indices for the observations of each device
  device_start_id=which(!duplicated(gps_data$userid))
  device_end_id=which(!duplicated(gps_data$userid,fromLast=TRUE))
  device_ind_array = cbind(device_start_id,device_end_id)
  rm(gps_data)
  
  cl<- makeCluster(20)
  registerDoSNOW(cl)
  getDoParWorkers()
  trips=foreach (core=c(0:19),.packages = c('hashids','geosphere','data.table'),.export = c("ind","device_ind_array"),.combine="rbind")%dopar%{
    trip_identification_function(core,ind,device_ind_array)
  }
  stopCluster(cl)
  rm(cl)
  closeAllConnections()
  gc()
  fwrite(trips,paste0("trips_",ind,".csv"))
}
  

