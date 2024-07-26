General Information: 
--------------------

GPR_SEAL stands for Ground Penetrating Radar Simple Easily Accessible Library.

This library contains a bunch of Python functions to process 
3D Ground Penetrating Radar data. The focus lies on simplicity. 
The goal is to provide a package that is easy to use for students 
new to the topic, such that they can be acquainted with 3D-GPR 
processing without frustration and get good results quickly without 
needing to write a lot of code.

Please use the script GPR_SEAL_example.py to get an example of the most 
simple processing flow that contains the absolute minimum of processing steps. 
Using the functions provided, the example flow can be extended with additional steps.

Check the readme file of the Matlab version for detailed information about the functions. 
Please note, that the Matlab version of GPR_SEAL is at the moment more extensive
than the Python version. More functions will be implemented in the future. 

Which functions are implemented?
--------------------------------

- plot_profile(yfrac,fs): Plot a profile
  * yfrac: Fractional position of measurement line in crossline direction, 
    where 0 is the first line, 0.5 the middle line and 1 the last line. 
  * fs: Fontsize to be used in plots. 

- time_zero_correction(tshift): Correct for time zero
  * tshift: The time shift in nanoseconds. A negative time shift removes
    the samples at the top, a positive time shift adds samples at the top. 

- subtract_mean_trace(xwin,tstart,tend): Subtract mean trace
  * xwin: The size of the moving window in meters. Traces inside this
    window are taken into account to calculate the mean trace. 
  * tstart: Start time in nanoseconds from which onwards the mean trace is 
    subtracted. 
  * tend: End time in nanoseconds before which the mean trace is
    subtracted.

- gain(startgain,lingain,expgain,doplotgain,fs): Apply gain function
  * startgain: Time in nanoseconds after which arrivals should be gained. 
  * lingain: Increase of linear gain (put to 0 to avoid any linear gain)
  * expgain: Increase of exponential gain (put to 0 to avoid any
    exponential gain)
  * doplotgain: Plot the gain function (1) or not (0)
  * fs: Fontsize for plot of gain function

- set_caxis(yfrac): Set color axis
  * yfrac: Fractional position of measurement line in crossline direction, 
    where 0 is the first line, 0.5 the middle line and 1 the last line. 

- velocity_analysis(yfrac,ant_sep): Velocity analysis
  * yfrac: Fractional position of measurement line in crossline direction, 
    where 0 is the first line, 0.5 the middle line and 1 the last line. 
  * ant_sep: Antenna separation (distance between source and receiver) in
    meters.

- browse_profiles(yfrac): 2D-Browser
  * yfrac: Fractional position of measurement line in crossline direction, 
    where 0 is the first line, 0.5 the middle line and 1 the last line. 

- browse(): 3D-Browser
  no input arguments

How to run the code: 
--------------------

Run the file GPR_SEAL_example.py in Python. In Linux, you can simply type in a terminal: 
python GPR_SEAL_example.py

Which Python packages are required?
-----------------------------------

time, numpy, matplotlib, sys, tkinter

Known bugs:
-----------
- When multiple interactive windows are open (also inactive ones), adjusting the 
  font size only works in the first one. 

License: 
--------

This file is part of the GPR_SEAL python package. GPR_SEAL is free software: 
you can redistribute it and/or modify it under the terms of the GNU General 
Public License as published by the Free Software Foundation, either 
version 3 of the License, or (at your option) any later version. 
GPR_SEAL is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
more details. You should have received a copy of the GNU General Public 
License along with GPR_SEAL. If not, see <https://www.gnu.org/licenses/>. 
