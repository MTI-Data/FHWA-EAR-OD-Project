This code assigns an age, gender, and income category to each device ID, based on imputed home Census Block Group.

The method can be summarized as below:
	1. Gather and clean the ACS data that has been downloaded from NHGIS.org
	2. Calculate the age, income, and gender distribution of each zone.
	3. Assign an age, gender and income category to each device ID, based on the distribution calculated.

NOTE:

Income Categories:	
		* Below 20K
		* Between 20K and 50K
		* Between 50K and 100K
		* Above 100K

Gender Categories:
		* MALE
		* FEMALE

Age Categories:
		* Below 35 years old
		* Between 35 and 65 years old
		* Above 65 years old
