# -*- coding: utf-8 -*-
"""

"""
#Loading the libraries
import pandas as pd
import os
from os import sys
import pickle

#setting the directory
os.chdir(sys.path[0])


#loading the data:
data = pd.read_csv('Trip_roster_Directory')
#adding mode attributes to the data
data['mode']=0
#Predicting air trips
data.loc[data.loc[(data['trip_dist']>=50000) & (data['speed_Q75']>=100)].index.values,'mode']=4

#separating air trips from other trips
airtrips=data.loc[data['mode']==4]
df=data.loc[data['mode']==0] 
#Loading data scaler model
datascaler=pickle.load(open('data_scaler.sav','rb'))

#Scaling test data
test_data=df[df.columns[2:34]]
test_data_scaled = datascaler.transform(test_data)


#loading the Random Forest model
RandomForest=pickle.load(open('Random_Forest.sav','rb'))
#Predicting other Modes
prediction=RandomForest.predict(test_data_scaled)

#adding the prediction results to the data
df.mode=prediction

#Combining all trips and saving
alltrips=df.append(airtrips)
alltrips=pd.DataFrame.sort_index(alltrips)
alltrips.to_csv('Output_Directory')
