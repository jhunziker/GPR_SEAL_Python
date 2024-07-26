#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of the GPR_SEAL python package. GPR_SEAL is free software:
you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version.
GPR_SEAL is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
more details. You should have received a copy of the GNU General Public
License along with GPR_SEAL. If not, see <https://www.gnu.org/licenses/>.
"""

import time
from lib.GPR_SEAL import GPR_SEAL

# Start a timer
tstart = time.time()

# Specify the folder where the GPR data are stored
foldername = './data/'

# Create a list of files that contain the profiles of the GPR data cube
filenamelist = []
for ifile in range(40):
    if ifile<9:
        filenamelist.append('beach_000' + str(ifile+1) + '_0.iprh')
    else:
        filenamelist.append('beach_00' + str(ifile+1) + '_0.iprh')

# Crossline distance between lines in meters
dy = 0.2

# Create a data object of class GPR_SEAL and import the data.
data = GPR_SEAL(foldername,filenamelist,dy)

# Correct time zero, such that the direct wave arrives at the correct time.
data.time_zero_correction(-8)

# Subtract a mean trace in a moving window to suppress the direct wave. 
data.subtract_mean_trace(2.0,0.0,8.5)

# Increase the amplitude of late arrivals. 
data.gain(8.5,0.25,0.02,0,40)
#data.browse_profiles(0.5)

# Determine the velocity of the background medium. 
data.velocity_analysis(0.5,0.23)

# Check if a velocity has been stored and a depth vector created
if hasattr(data,'zvec'):
    print('Stored velocity: %0.3f m/ns' % (data.vel))
    print('Maximum depth: %0.4f m' % (data.zvec[-1]))

# Plot a horizontal slice
data.browse()

# Stop timer and print runtime
tend = time.time()
print("Runtime: %.1f s" %(tend-tstart))
