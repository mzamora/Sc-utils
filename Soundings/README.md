# Soundings code
## Raw files
There are two sources of raw files, located in the `raw` folder.
The csv files are obtained from the [UWyoming website](http://weather.uwyo.edu/upperair/sounding.html), with an automated script to include all the data for a month.
The netCDF files are obtained manually from the [NOAA Raobs website](https://ruc.noaa.gov/raobs/).

## Reading soundings
To read the csv files, `Get_sounding_Var.m` will obtain all the variables for a selected timerange

To read the netCDF files, `getSounding_raobs.m` can read the profiles and derive other variables

## Processing sounding data
`process_all_NKX_soundings.m` is a script that goes through all needed dates reading the csv soundings and obtaining the most important parameters.

To do so, it also detects temperature inversions in the lower 3km of the atmosphere (where Sc is likely to occur), with the `TMP_Inversion_Strength_Cal.m` function.

___
Mónica Zamora Z., 2017. SRAF at UCSD solar.ucsd.edu
