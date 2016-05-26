# -*- coding: utf-8 -*-
"""
Created on Fri May 06 16:05:50 2016

@author: Cristina Gallego
"""

import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.image import AxesImage
import SimpleITK as sitk
import re
from dictionaries import mha_data_loc
import fnmatch

#!/usr/bin/env python

class Imgpick(object):
    """Imgage picking pixel coordinates and display functions """
    
    def __init__(self):
        self.x = None
        self.y = None
        self.imgArray = None
        self.diag = []
        self.annot_ax = False
        self.annot_axVoln = None
        self.sel_annots = []

    def __call__(self):       
        """ Turn Class into a callable object """
        Imgpick()
    
    def filen_patt_match(self, mha_data_loc, patternstr):
        selfilenm = []
        for filenm in os.listdir(mha_data_loc):
            if fnmatch.fnmatch(filenm, patternstr):
                selfilenm = filenm      
            
        return selfilenm
    

    def onpick(self, event):
        artist = event.artist
        if isinstance(artist, AxesImage):
            mouseevent = event.mouseevent
            self.x = mouseevent.xdata
            self.y = mouseevent.ydata
            print('image', [self.y,self.x])
            self.diag.append([self.y,self.x])
    
        
    def extract_annot(self, list_annots):
        '''Parse list of annotations, put markers according to notes and color code according to sequence order'''
        
        annots_dict_list=[]
        for one_annot in list_annots:
            annots_dict = {}
              
            # iterate throuhg attributes of annotation
            annots_dict['AccessionNumber'] = one_annot[one_annot.find("':")+4:one_annot.find("',")] 
            one_annot = one_annot[one_annot.find("',")+2:]
            
            annots_dict['SeriesDate'] = one_annot[one_annot.find("':")+4:one_annot.find("',")] 
            one_annot = one_annot[one_annot.find("',")+2:]
            
            annots_dict['SeriesNumber'] = one_annot[one_annot.find("':")+4:one_annot.find("',")] 
            one_annot = one_annot[one_annot.find("',")+2:]
            
            annots_dict['SliceLocation'] = one_annot[one_annot.find("':")+4:one_annot.find("',")] 
            one_annot = one_annot[one_annot.find("',")+2:]
            
            annots_dict['SeriesDescription'] = one_annot[one_annot.find("':")+4:one_annot.find("',")] 
            one_annot = one_annot[one_annot.find("',")+2:]
            
            annots_dict['PatientID'] = one_annot[one_annot.find("':")+4:one_annot.find("',")] 
            one_annot = one_annot[one_annot.find("',")+2:]
            
            annots_dict['StudyID'] = one_annot[one_annot.find("':")+4:one_annot.find("',")] 
            one_annot = one_annot[one_annot.find("',")+2:]
                            
            # get the type of annotation: e.g CALLIPER, ELLIPSE, ARROW
            annots_dict['note'] = one_annot[one_annot.find("':")+4:one_annot.find("\\")]
            
            # extract annotation coordinate location
            coords_str = one_annot[one_annot.find("\\"):one_annot.find("',")]
            non_dec = re.compile(r'[^\d.]+')
            coords = non_dec.sub(',', coords_str).split(',')
            if coords != ['']:
                annots_dict['xi']=float(coords[1])
                annots_dict['yi']=float(coords[2])
                annots_dict['xf']=float(coords[3])
                annots_dict['yf']=float(coords[4])
            
            # finish last attribute of annotations
            one_annot = one_annot[one_annot.find("',")+2:]
            annots_dict['SeriesInstanceUID'] = one_annot[one_annot.find("':")+4:]
            annots_dict_list.append(annots_dict)
                
        return annots_dict_list
        
        
    def selectDiagonal(self, imgArray, sno, annots_dict_list):
        """
        selectDiagonal takes to picked pixels for the diagonal of the patch 
        
        Usage:
        # plot slice location over time
        slicestime = [ Vol1[:,:,s_img_no-1], 
                       Vol2[:,:,s_img_no-1],
                       Vol3[:,:,s_img_no-1], 
                       Vol4[:,:,s_img_no-1]]          
        # display to pick pixels
        selcPatch.plotimgs(slicestime)
        if annotationsfound not [] dispplay mark on slice
        
        """
        self.fig = plt.figure()
        self.ax1 = self.fig.add_subplot(1, 4, 1)
        self.ax2 = self.fig.add_subplot(1, 4, 2)
        self.ax3 = self.fig.add_subplot(1, 4, 3)
        self.ax4 = self.fig.add_subplot(1, 4, 4)
        
        t1 = self.ax1.imshow(sitk.GetArrayFromImage(imgArray[0]), interpolation=None, picker=True)
        t1.set_cmap("gray")
        t2 = self.ax2.imshow(sitk.GetArrayFromImage(imgArray[1]), interpolation=None, picker=True)
        t2.set_cmap("gray")
        t3 = self.ax3.imshow(sitk.GetArrayFromImage(imgArray[2]), interpolation=None, picker=True)
        t3.set_cmap("gray")
        t4 = self.ax4.imshow(sitk.GetArrayFromImage(imgArray[3]), interpolation=None, picker=True)
        t4.set_cmap("gray")
        self.fig.canvas.mpl_connect('pick_event', self.onpick)

        # Go through list of annots, if found any, mathc sno and display
        self.sel_annots = []
        if(annots_dict_list):
            for annots_dict in annots_dict_list:
                StudyID = annots_dict['StudyID']
                AccessionN = annots_dict['AccessionNumber']
                annotSeries = annots_dict['SeriesNumber'] # '603'
                SliceLocation = float(annots_dict['SliceLocation']) # e.g -73.34339142
                
                # read Volume
                strtomatch = '*_'+AccessionN+'_'+annotSeries
                filen_Vol = self.filen_patt_match(mha_data_loc, strtomatch+'*.mha')
                Vol = sitk.ReadImage(mha_data_loc+os.sep+filen_Vol)
                
                # find annot slice No
                Sliceind = Vol.TransformPhysicalPointToIndex([SliceLocation, Vol.GetOrigin()[1], Vol.GetOrigin()[2]])
                print('Sliceind: ', Sliceind)
                print("====================")
                
                if(sno == int(Sliceind[2])): # does it match the slice of the enhan?
                    self.annot_ax = True
                    self.sel_annots = self.sel_annots.append(annots_dict)
                    self.annot_axVoln = annotSeries[-1] # last digit: 0,1,2,3,4
                
                if(self.annot_ax):
                    if(annots_dict['note'] == 'LINEARROW'):
                        if(self.annot_axVoln=='1'):
                            # xy (arrow tip) and xytext locations (text location) are in data coordinates.
                            self.ax1.annotate(' ', xy=(annots_dict['xf'], annots_dict['yf']), xytext=(annots_dict['xi'], annots_dict['yi']),
                                              arrowprops=dict(arrowstyle="->",linewidth = 5.,color = 'red'))
                        if(self.annot_axVoln=='2'):
                            # xy (arrow tip) and xytext locations (text location) are in data coordinates.
                            self.ax2.annotate('', xy=(annots_dict['xf'], annots_dict['yf']), xytext=(annots_dict['xi'], annots_dict['yi']),
                                              arrowprops=dict(arrowstyle="->",linewidth = 5.,color = 'red'))
                        if(self.annot_axVoln=='3'):
                            # xy (arrow tip) and xytext locations (text location) are in data coordinates.
                            self.ax3.annotate('', xy=(annots_dict['xf'], annots_dict['yf']), xytext=(annots_dict['xi'], annots_dict['yi']),
                                              arrowprops=dict(arrowstyle="->",linewidth = 5.,color = 'red'))
                        if(self.annot_axVoln=='4'):
                            # xy (arrow tip) and xytext locations (text location) are in data coordinates.
                            self.ax4.annotate('', xy=(annots_dict['xf'], annots_dict['yf']), xytext=(annots_dict['xi'], annots_dict['yi']),
                                              arrowprops=dict(arrowstyle="->",linewidth = 5.,color = 'red'))                    
                    
                    if(annots_dict['note'] == 'CALIPER'):
                        if(self.annot_axVoln=='1'):
                            # xy (arrow tip) and xytext locations (text location) are in data coordinates.
                            self.ax1.annotate(' ', xy=(annots_dict['xf'], annots_dict['yf']), xytext=(annots_dict['xi'], annots_dict['yi']),
                                              arrowprops=dict(arrowstyle="-",linewidth = 5.,color = 'red'))
                        if(self.annot_axVoln=='2'):
                            # xy (arrow tip) and xytext locations (text location) are in data coordinates.
                            self.ax2.annotate('', xy=(annots_dict['xf'], annots_dict['yf']), xytext=(annots_dict['xi'], annots_dict['yi']),
                                              arrowprops=dict(arrowstyle="-",linewidth = 5.,color = 'red'))
                        if(self.annot_axVoln=='3'):
                            # xy (arrow tip) and xytext locations (text location) are in data coordinates.
                            self.ax3.annotate('', xy=(annots_dict['xf'], annots_dict['yf']), xytext=(annots_dict['xi'], annots_dict['yi']),
                                              arrowprops=dict(arrowstyle="-",linewidth = 5.,color = 'red'))
                        if(self.annot_axVoln=='4'):
                            # xy (arrow tip) and xytext locations (text location) are in data coordinates.
                            self.ax4.annotate('', xy=(annots_dict['xf'], annots_dict['yf']), xytext=(annots_dict['xi'], annots_dict['yi']),
                                              arrowprops=dict(arrowstyle="-",linewidth = 5.,color = 'red'))                         
                        
        
        mng = plt.get_current_fig_manager()
        mng.window.showMaximized()
        plt.show(block=True)
        
        
    def myshow(self, Vol, sn, title=None):
        """support 3d images
        include a title
        use physical pixel size for axis labels (mm)
        show the image as gray values
        
        Usage:
        # show asingle slice
        Vol1 = sitk.Cast(Vol1,sitk.sitkFloat32)
        selcPatch.myshow(Vol1, s_img_no) 
        """
    
        spacing = Vol.GetSpacing()   
        nda = sitk.GetArrayFromImage(Vol)             
        ysize = nda.shape[1]
        xsize = nda.shape[2]
        # set up margin and dpi of display
        margin=0.05
        dpi=80
        
        # Make a figure big enough to accomodate an axis of xpixels by ypixels
        # as well as the ticklabels, etc...
        figsize = (1 + margin) * ysize / dpi, (1 + margin) * xsize / dpi
    
        fig = plt.figure(figsize=figsize, dpi=dpi)
        
        # Make the axis the right size...
        # Add an axes at position rect [left, bottom, width, height] 
        ax = fig.add_axes([margin, margin, 1 - 3*margin, 1 - 2*margin])
        extent = (0, xsize*spacing[1], ysize*spacing[0], 0)
        
        if(sn):
            t = ax.imshow(nda[sn-1,:,:], extent=extent, interpolation=None)
        else:
            t = ax.imshow(nda, extent=extent, interpolation=None)
        
        t.set_cmap("gray")
        if(title):
            plt.title(title)
        plt.show()
        
    
    def myshowSlices(self, slices, title=None, margin=0.01, dpi=80):
        """
        Usage:
        # plot slice location over time
        slicestime = [ Vol1[:,:,s_img_no-1], 
                       Vol2[:,:,s_img_no-1],
                       Vol3[:,:,s_img_no-1], 
                       Vol4[:,:,s_img_no-1]]  
        selcPatch.myshowSlices(sitk.Tile(slicestime, [4,1]), dpi=120)
        """
        nda = sitk.GetArrayFromImage(slices)
        spacing = slices.GetSpacing()
                    
        ysize = nda.shape[0]
        xsize = nda.shape[1]
       
        # Make a figure big enough to accomodate an axis of xpixels by ypixels
        # as well as the ticklabels, etc...
        figsize = (1 + margin) * ysize / dpi, (1 + margin) * xsize / dpi
    
        fig = plt.figure(figsize=figsize, dpi=dpi)
        # Make the axis the right size...
        ax = fig.add_axes([margin, margin, 1 - 2*margin, 1 - 2*margin])
        
        extent = (0, xsize*spacing[1], ysize*spacing[0], 0)
        
        t = ax.imshow(nda, extent=extent, interpolation=None)
        
        if nda.ndim == 2:
            t.set_cmap("gray")
        
        if(title):
            plt.title(title)
            
        mng = plt.get_current_fig_manager()
        mng.window.showMaximized()
        plt.show()

