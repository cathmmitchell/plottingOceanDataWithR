Plotting Ocean Data With R: Challenge Solutions
================

Here are the solutions to the Challenges in the main tutorial
`Plotting Ocean Data With R`

## Initialize session

We first need to load the `tidyverse` library:

``` r
library(tidyverse)
```

    ## -- Attaching packages ------------------------------------------------------------------------------------------------------------------------- tidyverse 1.3.0 --

    ## v ggplot2 3.3.2     v purrr   0.3.4
    ## v tibble  3.0.3     v dplyr   1.0.2
    ## v tidyr   1.1.2     v stringr 1.4.0
    ## v readr   1.3.1     v forcats 0.5.0

    ## -- Conflicts ---------------------------------------------------------------------------------------------------------------------------- tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

and the data:

``` r
fieldData <- read_csv('DamariscottaRiverData.csv')
```

    ## Parsed with column specification:
    ## cols(
    ##   date = col_double(),
    ##   station = col_double(),
    ##   depth_m = col_double(),
    ##   year = col_double(),
    ##   month = col_double(),
    ##   day = col_double(),
    ##   temperature_degC = col_double(),
    ##   salinity_psu = col_double(),
    ##   density_kg_m3 = col_double(),
    ##   PAR = col_double(),
    ##   fluorescence_mg_m3 = col_double(),
    ##   oxygenConc_umol_kg = col_double(),
    ##   oxygenSaturation_percent = col_double(),
    ##   latitude = col_double()
    ## )

## Challenge A

1.a. Create a scatter plot of ‘temperature\_degC’ by
‘fluorescence\_mg\_m3’ and color the points by ‘station’.

``` r
ggplot(data = fieldData, mapping = aes(x = temperature_degC, y = fluorescence_mg_m3)) + 
  geom_point(aes(color=station)) 
```

![](plottingOceanDataWithR-challengeSolutions_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->
1.b. Do the same as 1.a. but convert the station values to factors.
What’s the difference between the two plots?

``` r
ggplot(data = fieldData, mapping = aes(x = temperature_degC, y = fluorescence_mg_m3)) + 
  geom_point(aes(color=factor(station)))  
```

![](plottingOceanDataWithR-challengeSolutions_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

In the first plot, the colors are on a continuous colorscale. In the
second plot, the colors are discrete, separate colors. When the data was
read into R, the station values were read in as numeric values, so R
plotted them on a continuous color scale. But the station number can
really be though of as discrete data - they didn’t *need* to be numbers,
they could have been letters, or names. The station numbers are separate
from each other - not part of a continuous scale. So, to plot them as
four separate entities, we need to tell R to treat the station column as
factors (as it would have automatically done if the stations were named
“A”, “B”, “C”, and “D”).

2.  Create boxplots looking at the distribution of temperature\_degC by
    station (tip: change the station values to factors)

``` r
ggplot(data = fieldData, mapping = aes(x = factor(station), y = temperature_degC)) + 
  geom_boxplot() 
```

![](plottingOceanDataWithR-challengeSolutions_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

3.  Plot temperature by depth for samples from 2016, coloring the points
    by station

``` r
datasubset <- fieldData %>% filter(date==20160908)

ggplot(data = datasubset, mapping = aes(x = temperature_degC, y = depth_m)) + 
  geom_point(aes(color=factor(station))) 
```

![](plottingOceanDataWithR-challengeSolutions_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
\#\# Challenge B

What are the steps we need to take to manipulate our data into a data
frame with three columns: cruise, station and fluorescence averaged over
the top 2 m? Describe the step and say the associated `dplyr` data frame
manipulation functions you would use to do it.

1.  Select the rows where the depth is 2 m or less (use `filter`)
2.  Separate the data into date and station groups (use `group_by`)
3.  Take the average for each group of data (use `summarize`)

Use the pipes to connect all three steps together.

## Challenge C

1.  Create a `geom_tile` plot for temperature\_degC from the cruise that
    took place on Septh 12th 2017.

``` r
cruiseData <- filter(fieldData, date==20160908)
ggplot(cruiseData,aes(x=station,y=depth_m)) +
  geom_tile(aes(fill=temperature_degC)) +
  #scale_fill_continuous() +
  labs(fill='temperature (degC') +
  scale_y_reverse()
```

![](plottingOceanDataWithR-challengeSolutions_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

2.  Are there any examples from your own data that you could plot in
    this way?

[Here’s an example with genomic
data](https://science.sciencemag.org/content/sci/358/6366/1046/F2.large.jpg)
that was created using `geom_tile`.
