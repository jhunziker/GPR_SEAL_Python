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

import numpy as np
from matplotlib import pyplot as plt
import sys

import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# from mpl_toolkits import mplot3d

class GPR_SEAL:
    ###########################################################################
    # Initialize data cube and load impulseradar data                         #
    ###########################################################################
    def __init__(self,foldername,filenamelist,dy):
        # Load the header of the first file
        header = []
        with open(foldername + filenamelist[0], 'r') as fid:
            fid.seek(0,2) # Set the file cursor at the end of the file
            filesize = fid.tell()
            fid.seek(0,0) # Set the file cursor at the beginning of the file
            while fid.tell()<filesize:
                header.append(fid.readline())
                
        # Read the length of the time window from the header
        for iline in range(len(header)):
            if len(header[iline])>10:
                if header[iline][0:10]=='TIMEWINDOW':
                    timewindow = float(header[iline][11:-1])
                    break

        # Read the amount of samples in time from the header
        for iline in range(len(header)):
            if len(header[iline])>7:
                if header[iline][0:7]=='SAMPLES':
                    nsamples = int(header[iline][8:-1])
                    break

        # Read the inline spacing of traces from the header
        for iline in range(len(header)):
            if len(header[iline])>17:
                if header[iline][0:17]=='DISTANCE INTERVAL':
                    dx = float(header[iline][18:-1])
                    break

        # Determine the size of the longest profile
        for ifile in range(len(filenamelist)):
            # Load the header of each file
            tempheader = []
            with open(foldername + filenamelist[ifile], 'r') as fid:
                fid.seek(0,2) # Set the file cursor at the end of the file
                filesize = fid.tell()
                fid.seek(0,0) # Set the file cursor at the beginning of the file
                while fid.tell()<filesize:
                    tempheader.append(fid.readline())
            # Read the amount of traces in each profile and keep the longest
            for iline in range(len(tempheader)):
                if len(tempheader[iline])>10:
                    if tempheader[iline][0:10]=='LAST TRACE':
                        if ifile==0:
                            ntraces = int(tempheader[iline][11:-1])
                        else: 
                            temp_ntraces = int(tempheader[iline][11:-1])
                            if temp_ntraces>ntraces:
                                ntraces=temp_ntraces
                        break
        del temp_ntraces, tempheader

        # Setup vectors of data cube
        self.tvec = np.linspace(0,timewindow,nsamples)
        self.xvec = np.linspace(0,ntraces-1,ntraces)*dx
        self.yvec = np.linspace(0,len(filenamelist)-1,len(filenamelist))*dy

        # Initialize data cube
        self.cube = np.zeros((nsamples,ntraces,len(filenamelist)))

        # Load binary data and fill data cube
        for ifile in range(len(filenamelist)):
            # Load the header of each file
            tempheader = []
            with open(foldername + filenamelist[ifile], 'r') as fid:
                fid.seek(0,2) # Set the file cursor at the end of the file
                filesize = fid.tell()
                fid.seek(0,0) # Set the file cursor at the beginning of the file
                while fid.tell()<filesize:
                    tempheader.append(fid.readline())
            for iline in range(len(tempheader)):
                if len(tempheader[iline])>10:
                    if tempheader[iline][0:10]=='LAST TRACE':
                        thisfile_ntraces = int(tempheader[iline][11:-1])
            with open(foldername + filenamelist[ifile][0:-1] + 'b', 'rb') as fid:
                self.cube[:,0:thisfile_ntraces,ifile] = np.fromfile(fid, dtype=np.int32).reshape((thisfile_ntraces,nsamples)).transpose()
    
    ###########################################################################
    # Plot a profile                                                          #
    ###########################################################################
    def plot_profile(self,yfrac,fs):
        cm = 1/2.54  # convert centimeters to inches for figure size
        
        yel = round(yfrac*self.cube.shape[2]);
        if yel<0:
            yel=0
        elif yel>self.cube.shape[2]-1:
            yel = self.cube.shape[2]-1
            
        fig, ax = plt.subplots(figsize=(100*cm, 50*cm))
        if hasattr(self,'zvec'):
            imag=ax.imshow(self.cube[:,:,yel],
                           extent=(self.xvec[0],self.xvec[-1],self.zvec[-1],self.zvec[0]),
                           cmap=plt.colormaps['binary_r'])
            ax.set_ylabel('Depth (m)')
        else:
            imag=ax.imshow(self.cube[:,:,yel],
                           extent=(self.xvec[0],self.xvec[-1],self.tvec[-1],self.tvec[0]),
                           cmap=plt.colormaps['binary_r'])
            ax.set_ylabel('Time (ns)')
        ax.set_xlabel('Horizontal distance (m)')
        ax.set_title('Profile at %0.1f m crossline distance' % (self.yvec[yel]))
        ax.set_aspect('auto')
        if hasattr(self,'cax'):
            imag.set_clim(self.cax[0],self.cax[1])
        plt.colorbar(imag)
        plt.rc('font', size=fs) 
        plt.show()

    ###########################################################################
    # Correct for time zero                                                   #
    ###########################################################################
    def time_zero_correction(self,tshift):
        # Check if the supplied value for tshift is reasonable
        if abs(tshift)>self.tvec[-1]:
            print('Error: Supplied time shift is larger than length of time vector.')
            sys.exit()
            
        # Time sampling
        dt = self.tvec[1]-self.tvec[0]
        
        # Determine necessary shift in amount of elements
        elshift = round(abs(tshift)/dt)            
        
        # Carry out the shift by either deleting rows or adding rows of zeros
        if tshift<0:
            self.cube=self.cube[elshift:,:,:]
        else:
            self.cube=np.concatenate((np.zeros((elshift,self.cube.shape[1],self.cube.shape[2])), self.cube))
        
        # Adjusting the time vector correspondingly
        self.tvec=np.linspace(0,self.cube.shape[0]-1,self.cube.shape[0])*dt
        
        # Adjusting the depth vector if it exists
        if hasattr(self,'zvec'):
            self.zvec = self.tvec/2*self.v
            
    ###########################################################################
    # Subtract mean trace                                                     #
    ###########################################################################
    def subtract_mean_trace(self,xwin,tstart,tend):
        # Check if the size of the spatial window is larger than 0. 
        if xwin<=0.0:
            print('Error: The size of the moving window has to be a positive value.')
            sys.exit()
            
        # If tstart is larger than tend, switch them.
        if tstart>tend:
            temp = tend
            tend = tstart
            tstart = temp
            del temp

        # Determine dx, dt and amount of elements for input parameters.
        dx = self.xvec[1]-self.xvec[0]
        xwinel = round(xwin/dx)
        dt = self.tvec[1]-self.tvec[0]
        tstartel = round(tstart/dt)
        tendel = round(tend/dt)
            
        # Check that the start element is not negative. 
        if tstartel<0:
            tstartel=0
            
        # Check that the last element is not larger than the size of the time dimension.
        if tendel>self.cube.shape[0]:
            tendel=self.cube.shape[0]
            
        # Calculate and subtract the mean trace.
        for iy in range(self.cube.shape[2]):
            for ix in range(self.cube.shape[1]):
                if ix-round(xwinel/2)<0 and ix+round(xwinel/2)>self.cube.shape[1]:
                    meantrace = np.mean(self.cube[:,:,iy],axis=1)
                elif ix-round(xwinel/2)<0:
                    meantrace = np.mean(self.cube[:,0:ix+round(xwinel/2),iy],axis=1)
                elif ix+round(xwinel/2)>self.cube.shape[1]:
                    meantrace = np.mean(self.cube[:,ix-round(xwinel/2):,iy],axis=1)
                else:
                    meantrace = np.mean(self.cube[:,ix-round(xwinel/2):ix+round(xwinel/2),iy],axis=1)
                self.cube[tstartel:tendel,ix,iy] = self.cube[tstartel:tendel,ix,iy]-meantrace[tstartel:tendel]
                
    ###########################################################################
    # Apply gain function                                                     #
    ###########################################################################
    def gain(self,startgain,lingain,expgain,doplotgain,fs):
        # Check if the start time of the gain is reasonable. 
        if startgain<0:
            print('Error: The starting point of the gain function has to be larger than zero.')
            sys.exit()
        elif startgain>self.tvec[-1]:
            print('Error: The starting point of the gain function cannot be later than the last sample of the time vector.')
            sys.exit()
        
        # Check if the value for the linear gain is reasonable. 
        if lingain<0: 
            print('Error: The linear gain has to be larger or equal to zero.')
            sys.exit()
        
        # Check if the value for the exponential gain is reasonable. 
        if expgain<0:
            print('Error: The exponential gain has to be larger or equal to zero.')
            sys.exit()
            
        # Calculate gain function. 
        dt = self.tvec[1] - self.tvec[0]
        startgainel = round(startgain/dt)
        gainfun = np.ones(self.tvec.shape)
        gaincoords = np.linspace(1,self.tvec.shape[0]-startgainel,self.tvec.shape[0]-startgainel)
        gainfun[startgainel:] = gainfun[startgainel:] + gaincoords*lingain
        gainfun[startgainel:] = gainfun[startgainel:] + np.exp(gaincoords*expgain) - 1.0
        
        # Apply gain function. 
        for iy in range(self.cube.shape[2]):
            for ix in range(self.cube.shape[1]):
                self.cube[:,ix,iy] = self.cube[:,ix,iy]*gainfun
                
        # Plot gain function.
        if doplotgain==1:
            cm = 1/2.54  # convert centimeters to inches for figure size
            fig, ax = plt.subplots(figsize=(50*cm, 50*cm))
            ax.plot(gainfun,self.tvec)
            ax.set_ylim(self.tvec[0],self.tvec[-1])
            ax.set_xlabel('Gain (-)')
            ax.set_xlabel('Time (ns)')
            plt.gca().invert_yaxis()
            plt.rc('font', size=fs)
                
    ###########################################################################
    # Set color axis                                                          #
    ###########################################################################
    def set_caxis(self,yfrac):
        yel = round(yfrac*self.cube.shape[2]);
        if yel<0:
            yel=0
        elif yel>self.cube.shape[2]-1:
            yel = self.cube.shape[2]-1
            
        root = tk.Tk()
        app = setcax(self,yel,root)
        root.mainloop()
                
    ###########################################################################
    # Velocity analysis                                                       #
    ###########################################################################
    def velocity_analysis(self,yfrac,ant_sep):
        yel = round(yfrac*self.cube.shape[2]);
        if yel<0:
            yel=0
        elif yel>self.cube.shape[2]-1:
            yel = self.cube.shape[2]-1
            
        root = tk.Tk()
        app = vel_ana(self,yel,ant_sep,root)
        root.mainloop()
                
    ###########################################################################
    # 2D-Browser                                                              #
    ###########################################################################
    def browse_profiles(self,yfrac):
        yel = round(yfrac*self.cube.shape[2]);
        if yel<0:
            yel=0
        elif yel>self.cube.shape[2]-1:
            yel = self.cube.shape[2]-1
            
        root = tk.Tk()
        app = browse2D(self,yel,root)
        root.mainloop()

    ###########################################################################
    # 3D-Browser                                                              #
    ###########################################################################
    def browse(self):
        root = tk.Tk()
        app = browse3D(self,root)
        root.mainloop()

#-----------------------------------------------------------------------------#

class setcax:
    def __init__(self,data,yel,root):
        # Copy input parameters into data members
        self.data=data
        self.yel=yel
        self.root=root
        if hasattr(self.data,'cax')==0:
            self.data.cax = np.zeros((2)) 
            self.data.cax[0] = self.data.cube.min()
            self.data.cax[1] = self.data.cube.max()

        # Window properties
        root.title("Set color axis")

        # Buttons in window
        frm_nav = tk.Frame(master=self.root)
        btn_broad = tk.Button(master=frm_nav,text="Broaden",command=self.broaden)
        btn_broad.pack(fill=tk.X, side=tk.TOP)
        btn_narro = tk.Button(master=frm_nav,text="Narrow",command=self.narrow)
        btn_narro.pack(fill=tk.X, side=tk.TOP)
        btn_quit = tk.Button(master=frm_nav,text="quit",command=self.root.quit)
        btn_quit.pack(fill=tk.X, side=tk.TOP)
        frm_nav.pack(fill=tk.Y, side=tk.LEFT)

        # Figure space in window
        frm_fig = tk.Frame(master=self.root)
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=frm_fig)
        self.widget = self.canvas.get_tk_widget()
        self.widget.pack(fill=tk.BOTH, side=tk.TOP, expand=1)
        # Note: expand=1 is necessary to let the plot be scaled
        frm_fig.pack(fill=tk.BOTH, side=tk.LEFT, expand=1)

        # Update plot
        self.update_plot()

    def broaden(self):
        self.data.cax[0] = self.data.cax[0]*2.0
        self.data.cax[1] = self.data.cax[1]*2.0
        self.update_plot()

    def narrow(self):
        self.data.cax[0] = self.data.cax[0]/2.0
        self.data.cax[1] = self.data.cax[1]/2.0
        self.update_plot()

    def update_plot(self):
        if hasattr(self.data,'zvec'):
            depthvec = self.data.zvec
        else:
            depthvec = self.data.tvec
        # Inline profile
        self.ax.clear()
        self.imag=self.ax.imshow(self.data.cube[:,:,self.yel],
                  extent=(self.data.xvec[0],self.data.xvec[-1],
                          depthvec[-1],depthvec[0]),
                       cmap=plt.colormaps['binary_r'],
                       interpolation='none')
        self.imag.set_clim(self.data.cax[0],self.data.cax[1])
        self.ax.set_xlabel('Inline distance (m)')
        if hasattr(self.data,'zvec'):
            self.ax.set_ylabel('Depth (m)')
        else:
            self.ax.set_ylabel('Time (ns)')
        self.ax.set_title('Inline profile %d' % (self.yel))
        self.ax.set_aspect('auto')
        self.canvas.draw()

#-----------------------------------------------------------------------------#

class vel_ana:
    def __init__(self,data,yel,ant_sep,root):
        # Copy input parameters into data members
        self.data=data
        self.yel=yel
        self.ant_sep = ant_sep
        self.fs = 20
        self.font = tk.font.Font(family="Helvetica", size=self.fs, weight="bold")
        self.root=root
        if hasattr(self.data,'cax')==0:
            self.data.cax = np.zeros((2)) 
            self.data.cax[0] = self.data.cube.min()
            self.data.cax[1] = self.data.cube.max()
        self.halfxwinel = 80
        self.vel = 0.1
        self.xpickel = round(data.cube.shape[1]/2)
        self.tpickel = round(data.cube.shape[0]/8)

        # Window properties
        root.title("Velocity Analysis")

        # Buttons in window
        frm_nav = tk.Frame(master=self.root)
        btn_next = tk.Button(master=frm_nav,text="Next", font=self.font,
                             command=self.next_profile)
        btn_next.pack(fill=tk.X, side=tk.TOP)
        btn_prev = tk.Button(master=frm_nav,text="Previous", font=self.font,
                             command=self.prev_profile)
        btn_prev.pack(fill=tk.X, side=tk.TOP)
        btn_broad = tk.Button(master=frm_nav,text="Broaden", font=self.font,
                              command=self.broaden)
        btn_broad.pack(fill=tk.X, side=tk.TOP)
        btn_narro = tk.Button(master=frm_nav,text="Narrow", font=self.font,
                              command=self.narrow)
        btn_narro.pack(fill=tk.X, side=tk.TOP)
        btn_left = tk.Button(master=frm_nav,text="Left", font=self.font,
                             command=self.left)
        btn_left.pack(fill=tk.X, side=tk.TOP)
        btn_right = tk.Button(master=frm_nav,text="Right", font=self.font,
                              command=self.right)
        btn_right.pack(fill=tk.X, side=tk.TOP)
        btn_up = tk.Button(master=frm_nav,text="Up", font=self.font,
                           command=self.up)
        btn_up.pack(fill=tk.X, side=tk.TOP)
        btn_down = tk.Button(master=frm_nav,text="Down", font=self.font,
                             command=self.down)
        btn_down.pack(fill=tk.X, side=tk.TOP)
        btn_bigleft = tk.Button(master=frm_nav,text="Big left", font=self.font,
                                command=self.bigleft)
        btn_bigleft.pack(fill=tk.X, side=tk.TOP)
        btn_bigright = tk.Button(master=frm_nav,text="Big right", font=self.font,
                                 command=self.bigright)
        btn_bigright.pack(fill=tk.X, side=tk.TOP)
        btn_bigup = tk.Button(master=frm_nav,text="Big up", font=self.font,
                              command=self.bigup)
        btn_bigup.pack(fill=tk.X, side=tk.TOP)
        btn_bigdown = tk.Button(master=frm_nav,text="Big down", font=self.font,
                                command=self.bigdown)
        btn_bigdown.pack(fill=tk.X, side=tk.TOP)
        btn_winplus = tk.Button(master=frm_nav,text="Window +", font=self.font,
                                command=self.winplus)
        btn_winplus.pack(fill=tk.X, side=tk.TOP)
        btn_winmin = tk.Button(master=frm_nav,text="Window -", font=self.font,
                               command=self.winmin)
        btn_winmin.pack(fill=tk.X, side=tk.TOP)
        btn_velplus = tk.Button(master=frm_nav,text="Velocity +", font=self.font,
                                command=self.velplus)
        btn_velplus.pack(fill=tk.X, side=tk.TOP)
        btn_velmin = tk.Button(master=frm_nav,text="Velocity -", font=self.font,
                               command=self.velmin)
        btn_velmin.pack(fill=tk.X, side=tk.TOP)
        btn_setvel = tk.Button(master=frm_nav,text="Set velocity", font=self.font,
                               command=self.setvel)
        btn_setvel.pack(fill=tk.X, side=tk.TOP)
        btn_fontplus = tk.Button(master=frm_nav,text="font +", font=self.font, 
                                 command=self.fontplus)
        btn_fontplus.pack(fill=tk.X, side=tk.TOP)
        btn_fontmin = tk.Button(master=frm_nav,text="font -", font=self.font, 
                                command=self.fontmin)
        btn_fontmin.pack(fill=tk.X, side=tk.TOP)
        btn_quit = tk.Button(master=frm_nav,text="quit", font=self.font,
                             command=self.root.quit)
        btn_quit.pack(fill=tk.X, side=tk.TOP)
        frm_nav.pack(fill=tk.Y, side=tk.LEFT)

        # Figure space in window
        frm_fig = tk.Frame(master=self.root)
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=frm_fig)
        self.widget = self.canvas.get_tk_widget()
        self.widget.pack(fill=tk.BOTH, side=tk.TOP, expand=1)
        # Note: expand=1 is necessary to let the plot be scaled
        frm_fig.pack(fill=tk.BOTH, side=tk.LEFT, expand=1)

        # Update plot
        self.update_plot()

    def next_profile(self):
        self.yel = self.yel+1
        if self.yel==self.data.cube.shape[2]:
            self.yel=0
        self.update_plot()

    def prev_profile(self):
        self.yel = self.yel-1
        if self.yel<0:
            self.yel=self.data.cube.shape[2]-1
        self.update_plot()

    def broaden(self):
        self.data.cax[0] = self.data.cax[0]*2.0
        self.data.cax[1] = self.data.cax[1]*2.0
        self.update_plot()

    def narrow(self):
        self.data.cax[0] = self.data.cax[0]/2.0
        self.data.cax[1] = self.data.cax[1]/2.0
        self.update_plot()
    
    def left(self):
        self.xpickel = self.xpickel-1
        if self.xpickel<0:
            self.xpickel = self.data.cube.shape[1]-1
        self.update_plot()
    
    def right(self):
        self.xpickel = self.xpickel+1
        if self.xpickel==self.data.cube.shape[1]:
            self.xpickel = 0
        self.update_plot()
    
    def up(self):
        self.tpickel = self.tpickel-1
        if self.tpickel<0:
            self.tpickel = self.data.cube.shape[0]-1
        self.update_plot()
    
    def down(self):
        self.tpickel = self.tpickel+1
        if self.tpickel==self.data.cube.shape[0]:
            self.tpickel = 0
        self.update_plot()
    
    def bigleft(self):
        self.xpickel = self.xpickel-10
        if self.xpickel<0:
            self.xpickel = self.data.cube.shape[1]-1
        self.update_plot()
    
    def bigright(self):
        self.xpickel = self.xpickel+10
        if self.xpickel>=self.data.cube.shape[1]:
            self.xpickel = 0
        self.update_plot()
    
    def bigup(self):
        self.tpickel = self.tpickel-10
        if self.tpickel<0:
            self.tpickel = self.data.cube.shape[0]-1
        self.update_plot()
    
    def bigdown(self):
        self.tpickel = self.tpickel+10
        if self.tpickel>=self.data.cube.shape[0]:
            self.tpickel = 0
        self.update_plot()

    def winplus(self):
        self.halfxwinel = self.halfxwinel+10
        if self.halfxwinel>round(self.data.cube.shape[1]/2)-1:
            self.halfxwinel = round(self.data.cube.shape[1]/2)-1
        self.update_plot()
        
    def winmin(self):
        self.halfxwinel = self.halfxwinel-10
        if self.halfxwinel<10:
            self.halfxwinel=10
        self.update_plot()

    def velplus(self):
        self.vel = self.vel+0.005
        if self.vel>0.3:
            self.vel=0.3
        self.update_plot()

    def velmin(self):
        self.vel = self.vel-0.005
        if self.vel<0.01:
            self.vel=0.01
        self.update_plot()

    def setvel(self):
        self.data.vel = self.vel
        # When considering the antenna separation when calculating the
        # depth vector, the depth-increment becomes non-constant and values
        # of time equal to zero or close to zero become complex. Therefore,
        # the antenna separation is here ignored. 
        # self.data.zvec = np.sqrt((self.data.tvec/2.0*self.vel)**2-(self.ant_sep/2.0)**2)
        self.data.zvec = self.data.tvec/2.0*self.vel
        self.update_plot()

    def fontplus(self):
        self.fs = self.fs+2
        self.font.config(size=self.fs)
        self.update_plot()
        
    def fontmin(self):
        self.fs = self.fs-2
        if self.fs<8:
            self.fs=8
        self.font.config(size=self.fs)
        self.update_plot()

    def update_plot(self):
        # Calculate diffraction hyperbola
        # Convert time pick into depth
        # For picks close to the surface (close to t = 0) the argument of the sqrt
        # becomes negative when the antenna separation is taken into account. 
        # To avoid this problem, the depth is set to 0.0 in that case. 
        # Anyway, there should not be any reflections at smaller times,  
        # because a signal needs a minimum time to travel from the source to the receiver. 
        tpick = self.data.tvec[self.tpickel]
        temp = (tpick/2.0*self.vel)**2-(self.ant_sep/2.0)**2
        if temp<0:
            zpick = 0.0
        else: 
            zpick = np.sqrt(temp)
        # Set up horizontal vector
        xstartel = self.xpickel-self.halfxwinel
        xendel = self.xpickel+self.halfxwinel
        if xstartel<0:
            xstartel=0
        if xendel>=self.data.cube.shape[1]:
            xendel=self.data.cube.shape[1]-1
        x_hyper = np.linspace(self.data.xvec[xstartel],self.data.xvec[xendel],21)
        # Calculate hyperbola
        xpick = self.data.xvec[self.xpickel]
        t_hyper = np.sqrt(zpick**2+(abs(x_hyper-xpick)+self.ant_sep/2.0)**2)/self.vel
        t_hyper = t_hyper + np.sqrt(zpick**2+(abs(x_hyper-xpick)-self.ant_sep/2.0)**2)/self.vel
        # Plot profile
        self.ax.clear()
        self.imag=self.ax.imshow(self.data.cube[:,:,self.yel],
                  extent=(self.data.xvec[0],self.data.xvec[-1],
                          self.data.tvec[-1],self.data.tvec[0]),
                       cmap=plt.colormaps['binary_r'],
                       interpolation='none')
        self.imag.set_clim(self.data.cax[0],self.data.cax[1])
        # self.ax.plot(self.data.xvec[self.xpickel],self.data.tvec[self.tpickel], 
        #              marker='+', markersize=20, color="red") 
        self.ax.plot(x_hyper,t_hyper, linewidth=3, color="red") 
        self.ax.set_xlabel('Inline distance (m)',fontsize=self.fs)
        self.ax.set_ylabel('Depth (ns)',fontsize=self.fs)
        # self.ax.set_title('Inline profile %d' % (self.yel))
        self.ax.set_title('Current velocity: %0.3f m/ns' % (self.vel),fontsize=self.fs)
        self.ax.tick_params(axis='both', labelsize=self.fs)
        self.ax.set_aspect('auto')
        self.canvas.draw()

#-----------------------------------------------------------------------------#

class browse2D:
    def __init__(self,data,yel,root):
        # Copy input parameters into data members
        self.data=data
        self.yel=yel
        self.fs = 20
        self.font = tk.font.Font(family="Helvetica", size=self.fs, weight="bold")
        self.root=root
        if hasattr(self.data,'cax')==0:
            self.data.cax = np.zeros((2)) 
            self.data.cax[0] = self.data.cube.min()
            self.data.cax[1] = self.data.cube.max()

        # Window properties
        root.title("2D Browser")

        # Buttons in window
        frm_nav = tk.Frame(master=self.root)
        btn_next = tk.Button(master=frm_nav,text="Next", font=self.font,
                             command=self.next_profile)
        btn_next.pack(fill=tk.X, side=tk.TOP)
        btn_prev = tk.Button(master=frm_nav,text="Previous", font=self.font,
                             command=self.prev_profile)
        btn_prev.pack(fill=tk.X, side=tk.TOP)
        btn_broad = tk.Button(master=frm_nav,text="Broaden", font=self.font,
                              command=self.broaden)
        btn_broad.pack(fill=tk.X, side=tk.TOP)
        btn_narro = tk.Button(master=frm_nav,text="Narrow", font=self.font,
                              command=self.narrow)
        btn_narro.pack(fill=tk.X, side=tk.TOP)
        btn_fontplus = tk.Button(master=frm_nav,text="font +", font=self.font, 
                                 command=self.fontplus)
        btn_fontplus.pack(fill=tk.X, side=tk.TOP)
        btn_fontmin = tk.Button(master=frm_nav,text="font -", font=self.font, 
                                command=self.fontmin)
        btn_fontmin.pack(fill=tk.X, side=tk.TOP)
        btn_quit = tk.Button(master=frm_nav,text="quit", font=self.font,
                             command=self.root.quit)
        btn_quit.pack(fill=tk.X, side=tk.TOP)
        frm_nav.pack(fill=tk.Y, side=tk.LEFT)

        # Figure space in window
        frm_fig = tk.Frame(master=self.root)
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=frm_fig)
        self.widget = self.canvas.get_tk_widget()
        self.widget.pack(fill=tk.BOTH, side=tk.TOP, expand=1)
        # Note: expand=1 is necessary to let the plot be scaled
        frm_fig.pack(fill=tk.BOTH, side=tk.LEFT, expand=1)

        # Update plot
        self.update_plot()

    def next_profile(self):
        self.yel = self.yel+1
        if self.yel==self.data.cube.shape[2]:
            self.yel=0
        self.update_plot()

    def prev_profile(self):
        self.yel = self.yel-1
        if self.yel<0:
            self.yel=self.data.cube.shape[2]-1
        self.update_plot()

    def broaden(self):
        self.data.cax[0] = self.data.cax[0]*2.0
        self.data.cax[1] = self.data.cax[1]*2.0
        self.update_plot()

    def narrow(self):
        self.data.cax[0] = self.data.cax[0]/2.0
        self.data.cax[1] = self.data.cax[1]/2.0
        self.update_plot()

    def fontplus(self):
        self.fs = self.fs+2
        self.font.config(size=self.fs)
        self.update_plot()
        
    def fontmin(self):
        self.fs = self.fs-2
        if self.fs<8:
            self.fs=8
        self.font.config(size=self.fs)
        self.update_plot()

    def update_plot(self):
        if hasattr(self.data,'zvec'):
            depthvec = self.data.zvec
        else:
            depthvec = self.data.tvec
        # Inline profile
        self.ax.clear()
        self.imag=self.ax.imshow(self.data.cube[:,:,self.yel],
                  extent=(self.data.xvec[0],self.data.xvec[-1],
                          depthvec[-1],depthvec[0]),
                       cmap=plt.colormaps['binary_r'],
                       interpolation='none')
        self.imag.set_clim(self.data.cax[0],self.data.cax[1])
        self.ax.set_xlabel('Inline distance (m)',fontsize=self.fs)
        if hasattr(self.data,'zvec'):
            self.ax.set_ylabel('Depth (m)',fontsize=self.fs)
        else:
            self.ax.set_ylabel('Time (ns)',fontsize=self.fs)
        self.ax.set_title('Inline profile %d' % (self.yel),fontsize=self.fs)
        self.ax.tick_params(axis='both', labelsize=self.fs)
        self.ax.set_aspect('auto')
        self.canvas.draw()

#-----------------------------------------------------------------------------#

class browse3D:
    def __init__(self,data,root):
        # Copy input parameters into data members
        self.data=data
        self.root=root
        self.timel=round(data.cube.shape[0]/8)
        self.crossel=round(data.cube.shape[2]/2)
        self.inel=round(data.cube.shape[1]/2)
        self.doinv = 1
        self.fs = 20
        self.font = tk.font.Font(family="Helvetica", size=self.fs, weight="bold")

        if hasattr(self.data,'cax')==0:
            self.data.cax = np.zeros((2)) 
            self.data.cax[0] = self.data.cube.min()
            self.data.cax[1] = self.data.cube.max()

        # Window properties
        root.title("3D Browser")

        # First column of buttons
        frm_nav = tk.Frame(master=self.root)
        lbl11 = tk.Label(master=frm_nav, text="Navigation:", font=self.font)
        lbl11.pack(fill=tk.X, side=tk.TOP)
        btn_up = tk.Button(master=frm_nav,text="up", font=self.font,
                           command=self.moveup)
        btn_up.pack(fill=tk.X, side=tk.TOP)
        btn_do = tk.Button(master=frm_nav,text="down", font=self.font,
                           command=self.movedo)
        btn_do.pack(fill=tk.X, side=tk.TOP)
        btn_back = tk.Button(master=frm_nav,text="back", font=self.font,
                             command=self.moveback)
        btn_back.pack(fill=tk.X, side=tk.TOP)
        btn_front = tk.Button(master=frm_nav,text="front", font=self.font,
                              command=self.movefront)
        btn_front.pack(fill=tk.X, side=tk.TOP)
        btn_left = tk.Button(master=frm_nav,text="left", font=self.font,
                             command=self.moveleft)
        btn_left.pack(fill=tk.X, side=tk.TOP)
        btn_right = tk.Button(master=frm_nav,text="right", font=self.font,
                              command=self.moveright)
        btn_right.pack(fill=tk.X, side=tk.TOP)
        btn_bigleft = tk.Button(master=frm_nav,text="big left", font=self.font,
                                command=self.movebigleft)
        btn_bigleft.pack(fill=tk.X, side=tk.TOP)
        btn_bigright = tk.Button(master=frm_nav,text="big right", font=self.font,
                                 command=self.movebigright)
        btn_bigright.pack(fill=tk.X, side=tk.TOP)
        lbl12 = tk.Label(master=frm_nav, text="Color:", font=self.font)
        lbl12.pack(fill=tk.X, side=tk.TOP)
        btn_broad = tk.Button(master=frm_nav,text="broaden", font=self.font,
                              command=self.broaden)
        btn_broad.pack(fill=tk.X, side=tk.TOP)
        btn_narro = tk.Button(master=frm_nav,text="narrow", font=self.font,
                              command=self.narrow)
        btn_narro.pack(fill=tk.X, side=tk.TOP)
        lbl13 = tk.Label(master=frm_nav, text="Other:", font=self.font)
        lbl13.pack(fill=tk.X, side=tk.TOP)
        btn_dir = tk.Button(master=frm_nav,text="direction", font=self.font, command=self.direction)
        btn_dir.pack(fill=tk.X, side=tk.TOP)
        btn_fontplus = tk.Button(master=frm_nav,text="font +", font=self.font, command=self.fontplus)
        btn_fontplus.pack(fill=tk.X, side=tk.TOP)
        btn_fontmin = tk.Button(master=frm_nav,text="font -", font=self.font, command=self.fontmin)
        btn_fontmin.pack(fill=tk.X, side=tk.TOP)
        btn_quit = tk.Button(master=frm_nav,text="quit", font=self.font,
                             command=self.root.quit)
        btn_quit.pack(fill=tk.X, side=tk.TOP)
        frm_nav.pack(fill=tk.Y, side=tk.LEFT)

        # Figure space
        frm_fig = tk.Frame(master=self.root)
        self.fig, self.ax = plt.subplots(2,2, gridspec_kw={'width_ratios': [3, 1],
                                                           'height_ratios': [3, 1]})
        self.fig.delaxes(self.ax[1,1])
        self.canvas = FigureCanvasTkAgg(self.fig, master=frm_fig)
        self.widget = self.canvas.get_tk_widget()
        self.widget.pack(fill=tk.BOTH, side=tk.TOP, expand=1)
        # Note: expand=1 is necessary to let the plot be scaled
        frm_fig.pack(fill=tk.BOTH, side=tk.LEFT, expand=1)

        # Second column of buttons
        frm_nav2 = tk.Frame(master=self.root)
        lbl21 = tk.Label(master=frm_nav2, text="Jump:", font=self.font)
        lbl21.pack(fill=tk.X, side=tk.TOP)
        btn_j0 = tk.Button(master=frm_nav2,text="0/8", font=self.font, command=self.jump0)
        btn_j0.pack(fill=tk.X, side=tk.TOP)
        btn_j1 = tk.Button(master=frm_nav2,text="1/8", font=self.font, command=self.jump1)
        btn_j1.pack(fill=tk.X, side=tk.TOP)
        btn_j2 = tk.Button(master=frm_nav2,text="2/8", font=self.font, command=self.jump2)
        btn_j2.pack(fill=tk.X, side=tk.TOP)
        btn_j3 = tk.Button(master=frm_nav2,text="3/8", font=self.font, command=self.jump3)
        btn_j3.pack(fill=tk.X, side=tk.TOP)
        btn_j4 = tk.Button(master=frm_nav2,text="4/8", font=self.font, command=self.jump4)
        btn_j4.pack(fill=tk.X, side=tk.TOP)
        btn_j5 = tk.Button(master=frm_nav2,text="5/8", font=self.font, command=self.jump5)
        btn_j5.pack(fill=tk.X, side=tk.TOP)
        btn_j6 = tk.Button(master=frm_nav2,text="6/8", font=self.font, command=self.jump6)
        btn_j6.pack(fill=tk.X, side=tk.TOP)
        btn_j7 = tk.Button(master=frm_nav2,text="7/8", font=self.font, command=self.jump7)
        btn_j7.pack(fill=tk.X, side=tk.TOP)
        btn_j8 = tk.Button(master=frm_nav2,text="8/8", font=self.font, command=self.jump8)
        btn_j8.pack(fill=tk.X, side=tk.TOP)
        frm_nav2.pack(fill=tk.Y, side=tk.LEFT)

        # Update plot
        self.update_plot()

    def moveup(self):
        self.timel = self.timel-1
        if self.timel<0:
            self.timel=self.data.cube.shape[0]-1
        self.update_plot()

    def movedo(self):
        self.timel = self.timel+1
        if self.timel==self.data.cube.shape[0]:
            self.timel=0
        self.update_plot()

    def moveback(self):
        self.crossel = self.crossel-1
        if self.crossel<0:
            self.crossel=self.data.cube.shape[2]-1
        self.update_plot()

    def movefront(self):
        self.crossel = self.crossel+1
        if self.crossel==self.data.cube.shape[2]:
            self.crossel=0
        self.update_plot()

    def moveleft(self):
        self.inel = self.inel-1
        if self.inel<0:
            self.inel=self.data.cube.shape[1]-1
        self.update_plot()

    def moveright(self):
        self.inel = self.inel+1
        if self.inel==self.data.cube.shape[1]:
            self.inel=0
        self.update_plot()

    def movebigleft(self):
        self.inel = self.inel-10
        if self.inel<0:
            self.inel=self.data.cube.shape[1]-1
        self.update_plot()

    def movebigright(self):
        self.inel = self.inel+10
        if self.inel>=self.data.cube.shape[1]:
            self.inel=0
        self.update_plot()

    def broaden(self):
        self.data.cax[0] = self.data.cax[0]*2.0
        self.data.cax[1] = self.data.cax[1]*2.0
        self.update_plot()

    def narrow(self):
        self.data.cax[0] = self.data.cax[0]/2.0
        self.data.cax[1] = self.data.cax[1]/2.0
        self.update_plot()

    def jump0(self):
        self.timel = 0
        self.update_plot()

    def jump1(self):
        self.timel=round(self.data.cube.shape[0]/8*1)
        self.update_plot()

    def jump2(self):
        self.timel=round(self.data.cube.shape[0]/8*2)
        self.update_plot()

    def jump3(self):
        self.timel=round(self.data.cube.shape[0]/8*3)
        self.update_plot()

    def jump4(self):
        self.timel=round(self.data.cube.shape[0]/8*4)
        self.update_plot()

    def jump5(self):
        self.timel=round(self.data.cube.shape[0]/8*5)
        self.update_plot()

    def jump6(self):
        self.timel=round(self.data.cube.shape[0]/8*6)
        self.update_plot()

    def jump7(self):
        self.timel=round(self.data.cube.shape[0]/8*7)
        self.update_plot()

    def jump8(self):
        self.timel=round(self.data.cube.shape[0]-1)
        self.update_plot()
        
    def direction(self):
        if self.doinv==0:
            self.doinv=1
        else:
            self.doinv=0
        self.update_plot()
        
    def fontplus(self):
        self.fs = self.fs+2
        self.font.config(size=self.fs)
        self.update_plot()
        
    def fontmin(self):
        self.fs = self.fs-2
        if self.fs<8:
            self.fs=8
        self.font.config(size=self.fs)
        self.update_plot()

    def update_plot(self):
        dx = self.data.xvec[1]-self.data.xvec[0]
        dy = self.data.yvec[1]-self.data.yvec[0]
        dt = self.data.tvec[1]-self.data.tvec[0]
        xip = [self.data.xvec[0],self.data.xvec[-1]]
        xcp = [self.data.yvec[0],self.data.yvec[-1]]
        yip = [self.data.yvec[self.crossel],self.data.yvec[self.crossel]]
        ycp = [self.data.xvec[self.inel],self.data.xvec[self.inel]]
        if hasattr(self.data,'zvec'):
            depthvec = self.data.zvec
            xhs = [self.data.zvec[0],self.data.zvec[-1]]
            yhs = [self.data.zvec[self.timel],self.data.zvec[self.timel]]
        else:
            depthvec = self.data.tvec
            xhs = [self.data.tvec[0],self.data.tvec[-1]]
            yhs = [self.data.tvec[self.timel],self.data.tvec[self.timel]]
        # Horizontal slice
        self.ax[0,0].clear()
        if self.doinv==1:
            self.imag_horz=self.ax[0,0].imshow(self.data.cube[self.timel,:,:].transpose(),
                           extent=(self.data.xvec[0]-dx/2,self.data.xvec[-1]+dx/2,
                                   self.data.yvec[-1]+dy/2,self.data.yvec[0]-dy/2),
                           cmap=plt.colormaps['binary_r'],
                           interpolation='none')
        else: 
            self.imag_horz=self.ax[0,0].imshow(np.flip(self.data.cube[self.timel,:,:]
                                                 .transpose(),axis=0),
                           extent=(self.data.xvec[0]-dx/2,self.data.xvec[-1]+dx/2,
                                   self.data.yvec[0]-dy/2,self.data.yvec[-1]+dy/2),
                           cmap=plt.colormaps['binary_r'],
                           interpolation='none')
        if hasattr(self.data,'cax'):
            self.imag_horz.set_clim(self.data.cax[0],self.data.cax[1])
        self.ax[0,0].plot(xip, yip, linewidth=2, color="red") 
        self.ax[0,0].plot(ycp, xcp, linewidth=2, color="red") 
        self.ax[0,0].set_xlabel('Inline distance (m)',fontsize=self.fs)
        self.ax[0,0].set_ylabel('Crossline distance (m)',fontsize=self.fs)
        #self.ax[0,0].set_title('Horizontal section at time element %d' % (self.timel),fontsize=self.fs)
        self.ax[0,0].tick_params(axis='both', labelsize=self.fs)
        self.ax[0,0].set_aspect('auto')
        # Inline profile
        self.ax[1,0].clear()
        self.imag_profin=self.ax[1,0].imshow(self.data.cube[:,:,self.crossel],
                              extent=(self.data.xvec[0],self.data.xvec[-1],
                                      depthvec[-1],depthvec[0]),
                       cmap=plt.colormaps['binary_r'],
                       interpolation='none')
        if hasattr(self.data,'cax'):
            self.imag_profin.set_clim(self.data.cax[0],self.data.cax[1])
        self.ax[1,0].plot(xip, yhs, linewidth=2, color="red") 
        self.ax[1,0].plot(ycp, xhs, linewidth=2, color="red") 
        self.ax[1,0].set_xlabel('Inline distance (m)',fontsize=self.fs)
        if hasattr(self.data,'zvec'):
            self.ax[1,0].set_ylabel('Depth (m)',fontsize=self.fs)
        else:
            self.ax[1,0].set_ylabel('Time (ns)',fontsize=self.fs)
        #self.ax[1,0].set_title('Inline profile %d' % (self.crossel),fontsize=self.fs)
        self.ax[1,0].tick_params(axis='both', labelsize=self.fs)
        self.ax[1,0].set_aspect('auto')
        # Crossline profile
        self.ax[0,1].clear()
        if self.doinv==1:
            self.imag_profcross=self.ax[0,1].imshow(self.data.cube[:,self.inel,:].transpose(),
                                     extent=(depthvec[0],depthvec[-1],
                                             self.data.yvec[-1],self.data.yvec[0]),
                           cmap=plt.colormaps['binary_r'],
                           interpolation='none')
        else:
            self.imag_profcross=self.ax[0,1].imshow(np.flip(self.data.cube[:,self.inel,:]
                                                      .transpose(),axis=0),
                                     extent=(depthvec[0],depthvec[-1],
                                             self.data.yvec[0],self.data.yvec[-1]),
                           cmap=plt.colormaps['binary_r'],
                           interpolation='none')

        if hasattr(self.data,'cax'):
            self.imag_profcross.set_clim(self.data.cax[0],self.data.cax[1])
        self.ax[0,1].plot(yhs, xcp, linewidth=2, color="red") 
        self.ax[0,1].plot(xhs, yip, linewidth=2, color="red") 
        if hasattr(self.data,'zvec'):
            self.ax[0,1].set_xlabel('Depth (m)',fontsize=self.fs)
        else:
            self.ax[0,1].set_xlabel('Time (ns)',fontsize=self.fs)
        self.ax[0,1].set_ylabel('Crossline distance (m)',fontsize=self.fs)
        #self.ax[0,1].set_title('Crossline profile %d' % (self.inel),fontsize=self.fs)
        self.ax[0,1].tick_params(axis='both', labelsize=self.fs)
        self.ax[0,1].set_aspect('auto')
        self.canvas.draw()
