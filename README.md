matlab_utilities
================

Short MATLAB functions for reading, analysing and deriving quantities from atmospheric data. 

Inventory:

WRF UTILITIES
read_wrf_full.m
Reads WRF NetCDF output at a given time interval and extracts the full 3D field given a file path and name, 
variable name, starting time index, and time interval. Detects if a field is staggered and unstaggers it 
to mass points. For use with "new" MATLAB NetCDF libraries, tested on version R2014a.

read_wrf_height.m
Reads WRF NetCDF output at a at a given time index and extracts a 2D field interpolated to 
the height (in m) of your choosing given a file path and name, variable name, time index, and height. 
Detects if a field is staggered and unstaggers it to mass points. For use with "new" MATLAB NetCDF libraries, 
tested on version R2014a.

RADAR UTILITIES

get_cfradial_filenames.m
cfrad NetCDF files have names based on unevenly spaced time intervals. This reads those names from a text file created by a short python script of my creation. For opening cfrad files in MATLAB.

radar_pdf.m
Create a joint frequency distribution of horizontal reflectivity and differential reflectivity over a time and height region. Reads radar or WRF output. Normalization options. This hasn't been cleaned up for user-friendliness.

cfad.m
Create a contoured frequency by altitude diagram (CFAD) of chosen variable. Reads radar or WRF output. Normalization options. This hasn't been cleaned up for user-friendliness.

radar_mosaic.m
Combines data from more than one radar - takes the data with the highest reflectivity at each pixel to create a composite over area covered by at least one radar. Some hard-coded file reading sections will have to be modified.

conv_stratiform_partition.m
Convective stratiform partition based on reflectivity. For looping over multiple times, creates a separate netCDF file for output.

METEOROLGY UTILITIES

equiv_theta.m
calculate equivalent potential temperature from potential temperature, pressure, water vapor mixing ratio.

CM1 UTILITIES

read_cm1_traj.m
read trajectories from CM1

read_cm1r17.m

read standard output from CM1 release 17

