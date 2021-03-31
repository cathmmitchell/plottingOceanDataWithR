### Damariscotta River cruise data

















### Challenge A
##
## 1.a. Create a scatter plot of fluorescence_mg_m3 (x-axis) by temperature_degC 
## (y-axis) and color the points by station.
##
## 1.b. Do the same as 1.a. but convert the station values to factors. What is 
## the difference between the two plots?
##
## 2. Create boxplots looking at the distribution of temperature_degC by 
## (tip: change the station values to factors).
##
## 3. Plot temperature by depth for samples for the cruise from September 8th 
## 2016, coloring the points by station (tip: filter the data frame on the date column)















### Visualizing relationships between multiple variables


#### Data manipulation






### Challenge B
##
## What are the steps we need to take to manipulate our data into a data frame 
## with three columns: cruise, station and fluorescence averaged over the top
## 2 m? Describe the step and say the associated dplyr data frame manipulation 
## functions you would use to do it.







#### Visualizing using contour plots






#### Visualizing using bubble plots





### Interpolating and visualizing data





### Challenge C
## 1. Create a geom_tile plot for temperature_degC from the cruise that took place on Sept 12th 2017.
##
## 2. Are there any examples from your own data that you could plot in this way?




#### Interpolation








### Challenge D
## Create a plot like the above that shows temperature interpolated by depth and latitude
## for the cruise that took place on Sept 12th 2017.



### Functions and for loops for creating multiple versions of figures


#### Functions






## Challenge E
## Create a plot like the above that shows fluorescence interpolated by depth and 
## latitude for every cruise in 2016 (one plot per cruise). To get a list of all 
## the dates of cruises in 2016 do unique(data2016$date) and look at the object
## printed to the console





#### For loops










## Challenge F
##
## Adapt the cruiseInterpolationPlot function such that the x, y and color bar
## limits are determined from input arguments and re-plot all the fluorescence 
## data for each cruise in 2016. (Hint: use the xlim and ylim functions with 
## ggplot and include the limits argument in the scale_fill_distiller function.)






