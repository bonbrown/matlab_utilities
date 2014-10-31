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
