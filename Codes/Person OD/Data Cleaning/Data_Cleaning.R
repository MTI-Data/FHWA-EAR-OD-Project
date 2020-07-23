#Data cleaning

#Reading the data and extracting active device list for the specific month
Raw_Data_Path="Raw Data directoy"
#Getting the list of raw data files
Raw_Data_list=list.files(path=Raw_Data_Path,full.names = TRUE)
#Column name in the data
Column_name=c("Timestamp","Device_ID","Latitude","Longitude","Accuracy","Time_zone_offset") #based on raw data schema
#Integrating one month of data
monthly_raw_data=data.table()
for( i in 1:length(Raw_Data_list)){
  daily_raw_data=fread(Raw_Data_list[i], fill = TRUE)
  colnames(daily_raw_data)=Column_name
  monthly_raw_data=rbind(monthly_raw_data,daily_raw_data)
}
fwrite(monthly_raw_data,"Directory to device list ID")

#Cleaning based on Accuracy and duplication
#Removing data based on accuracy
accuracy_threshold=XXX
monthly_raw_data=monthly_raw_data[Accuracy<=accuracy_threshold]
#Remove duplicate observations
monthly_raw_data=unique(monthly_raw_data,by=c("Device_ID","Timestamp"))

#Sorting data based on device ID and timestamp
monthly_raw_data=monthly_raw_data[order(Device_ID,Timestamp),]
#adding ID to each observation (for binding projected coordinates in future)
monthly_raw_data$ID=seq.int(nrow(monthly_raw_data))

#Projecting data to utm
longlattoutm=function (x,y,ID,zone){
  xy=data.frame(ID = ID, X=x, Y=y)
  coordinates(xy)=c("X","Y")
  proj4string(xy)=CRS("+proj=longlat +datum=WGS84")
  res = spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
  return(as.data.table(res))
}
Projection=longlattoutm(monthly_raw_data$Longitude,monthly_raw_data$Latitude, ID=monthly_raw_data$ID,zone)
monthly_raw_data=merge(monthly_raw_data,Projection,by="ID")

#find the starting and ending index for each device ID
device_start_id=which(!duplicated(monthly_raw_data$Device_ID))
device_end_id=which(!duplicated(monthly_raw_data$Device_ID,fromLast = TRUE))
device_index_array=cbind(device_start_id,device_end_id)


#Nighttime activity cluster

#Defining DBSCAN hyper parameter
Epsilon<-XXX
minimum_number_of_points<-XXX

#activity cluster identification function
activity_cluster_identification=function(trip_index,file){
  tmp_df=file[trip_index[1]:trip_index[2],]
  tmp_df=as.data.table(tmp_df)
  npoint=nrow(tmp_df)
  Device_ID=tmp_df$Device_ID[1]
  coordination=tmp_df[,c("X","Y")]
  cluster_result=dbscan(coordination,eps=Epsilon,minPts = minimum_number_of_points)
  no_cluster=length(unique(cluster_result$cluster))
  attr_vec=c(npoint,Device_ID,no_cluster)
  attr_vec
}

activity_cluster_file=as.data.frame(t(apply(device_ind_array, 1, activity_cluster_identification,file=monthly_raw_data)))
colnames(activity_cluster_file)=c("npoint","Device_ID","no_cluster")
activity_cluster_file=as.data.table(activity_cluster_file)
active_device_list=activity_cluster_file[no_cluster>=2]
cleaned_data=active_device_list[device_ID%in%active_device_list$Device_ID]
fwrite(cleaned_data,"Cleaned Data Directory")

