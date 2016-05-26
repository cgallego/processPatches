# -*- coding: utf-8 -*-
"""
Created on Mon May 09 12:22:14 2016

@author: Cristina Gallego
"""

import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import SimpleITK as sitk
import math

from skimage import data, img_as_float
from skimage import exposure
import matplotlib.patches as mpatches
import matplotlib.ticker as plticker

from dictionaries import snapshot_loc

#!/usr/bin/env python
class Patchify(object):
    """Imgage Patchify coordinates and extract Patches functions """
    
    def __init__(self):
        self.ptcsize = None
        self.img_size = None

    def __call__(self):       
        """ Turn Class into a callable object """
        Patchify()
        
    def contrast_stretch(self, img):  
        ## Floating point images are between 0 and 1 (unsigned images) or -1 and 1 (signed images), 
        # while 8-bit images are expected to have values in {0,255}.
        #imgorig = img.view(np.uint8)
        p2, p999 = np.percentile(img, (0, 99.9))
        # Contrast stretching
        img_rescale = exposure.rescale_intensity(img, in_range=(p2, p999))
        # rescale to values again between o and 1 (unsigned images) 
        #img_rescale = imgres.view(np.float32)
        
        return img_rescale
        
            
    def extractPatches(self, imgArray, x1, y1, x2, y2, name):
        """
        Currently only supports rectangular patches
        extractPatch takes pixel coordinates of 2 points defining patch diagonals
        e.g: x1, y1, x2, y2
        
        Usage:
        # plot slice location over time
        extrPatch = Patchify()
        extrPatch.extractPatch(imgArray, x1, y1, x2, y2)
        
        imgArray:: contains preimg, img1 to img4 (total 5 slices)
        """
        # get img size
        self.img_size = imgArray[0].shape
        self.ptcsize = imgArray[0][x1:x2,y1:y2].shape # patch1.flatten().reshape(64L, 52L) == patch1
        
        patches = []
        # show       
        fig, axes = plt.subplots(ncols=2, nrows=5, figsize=(4, 4))
        a = axes.flat 
        
        # finally append all elemnet arrays
        allpatches = []
        
        ## for the post-contract imgs
        for k in range(0,len(imgArray)):
            # extract patch inside the rectangular ROI
            patchk = imgArray[k][x1:x2,y1:y2]
            patches.append( patchk.flatten() ) 
            a[2*k].imshow(imgArray[k], cmap=plt.cm.gray) 
            
            allpatches = np.insert(allpatches, len(allpatches), patches[k])
            # reshape and extract
            eimg = allpatches.reshape(k+1, self.ptcsize[1], self.ptcsize[0])
            img = eimg[k,:,:]
            eslice = img.reshape(self.ptcsize[0],self.ptcsize[1])
            #plot patch
            a[2*k+1].imshow(eslice, cmap=plt.cm.gray) 
            
        # display
        mng = plt.get_current_fig_manager()
        mng.window.showMaximized()
        plt.savefig(snapshot_loc+os.sep+name+'.pdf')
        plt.show(block=True) 
        
        return allpatches
        
        
    def extractUnlabeledPatches(self, pimg1, pimg2, pimg3, pimg4, selectedsl, selectednx, selectedny, ha, wa, name):
        """
        Currently only supports rectangular patches
        extractPatch takes pimg0, pimg1, pimg2, pimg3, pimg4 vol arrays and a list of selected slices, 
        selected x,y init patch locations
        
        Usage:
        # plots patches location over time
        extrPatch = Patchify()
        extrPatch.extractUnlabeledPatches(pimg0, pimg1, pimg2, pimg3, pimg4, selectedsl, selectednx, selectedny, name):
        
        preimg, img1 to img4 (total 5 slices)
        self.ptcsize = imgslicestime[0][x1:x2,y1:y2].shape # patch1.flatten().reshape(64L, 52L) == patch1
        
        """        
        self.ptcsize = [ha,wa]
        allUpatches = []
        ncols = len(selectednx)
        nrows = len(selectedsl)
 
        # show   
        fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=(4, 4), sharex=True, sharey=True)
        # add a big axes, hide frame
        # hide tick and tick label of the big axes
        ax = axes.flat 
            
        for ssl,i in zip(selectedsl,range(len(selectedsl))):
            # extract slice
            imgslicestime = [ pimg1[ssl,:,:], pimg2[ssl,:,:], pimg3[ssl,:,:], pimg4[ssl,:,:]]
    
            for snx,j in zip(selectednx,range(len(selectednx))):
                # select nx, ny for patches
                x1 = snx
                x2 = x1+ha
                y1 = selectedny[j]
                y2 = y1+wa
                
                # finally append all elemnet arrays
                allpatches = []
                
                ## for the post-contract imgs
                for k in range(len(imgslicestime)): 
                    # extract patch inside the rectangular ROI
                    patchk = imgslicestime[k][x1:x2,y1:y2]
                    allpatches = np.insert(allpatches, len(allpatches), patchk.flatten())
                    
                    # reshape and extract
                    eimg = allpatches.reshape(k+1, self.ptcsize[1], self.ptcsize[0])
                    img = eimg[k,:,:]
                    eslice = img.reshape(self.ptcsize[0],self.ptcsize[1])
                
                # append extracted patch of size 3600L (30x30x4) and append a total of nrowsxcols = 50 
                allUpatches.append(allpatches)
                #plot patch of time point 4
                ax[(ncols)*i+j].imshow(eslice, cmap=plt.cm.gray) 
              
                
        # display
        mng = plt.get_current_fig_manager()
        mng.window.showMaximized()
        plt.savefig(snapshot_loc+os.sep+name+'.pdf')
        #plt.show(block=True) 
        plt.close()
        
        # TODO: put your data in a format that Theano can work with.
        # train_set_x = theano.shared(numpy.array(my_x, dtype='float64'))
        # train_set_y = theano.shared(numpy.array(my_y, dtype='float64'))
        # train_set_x = theano.shared(np.array(allUpatches, dtype='float64'))
        
        return allUpatches
        
        
    
    def extractROIPatches(self, pimg1, pimg2, pimg3, pimg4, centroid, patch_size, 
                          patch_diag1, patch_diag2, ha, wa, name):
        """
        Currently only supports rectangular patches
        extractPatch takes pixel coordinates of 2 points defining patch diagonals
        e.g: x1, y1, x2, y2
        
        Usage:
        # plot slice location over time
        extrPatch = Patchify()
        extrPatch.extractPatch(imgArray, x1, y1, x2, y2)
        
        imgArray:: contains preimg, img1 to img4 (total 5 slices)
        """
                
        # get patches sizes and quatities
        sno = int(centroid[2])-1
        x0 = float(centroid[0])
        y0 = float(centroid[1])
        hp = float(long(patch_size[0]))
        wp = float(long(patch_size[1]))
        self.optcsize = [hp-1, wp-1]
        self.ptcsize = [ha, wa]
        
        x1 = float(patch_diag1[0])
        y1 = float(patch_diag1[1])
        x2 = float(patch_diag2[0])
        y2 = float(patch_diag2[1])
        if(x1<x2):
            x1 = x1-10
            x2 = x2+10
        if(x1>x2):
            xt = x1
            x1 = x2-10
            x2 = xt+10
        if(y1<y2):
            y1 = y1-10
            y2 = y2+10
        if(y1>y2):
            yt = y1
            y1 = y2-10
            y2 = yt+10
        
        print '[x1:x2][y1:y2] [%d,%d][%d,%d]' % (x1,x2,y1,y2)
        x0 = x1 + (x2-x1)/2
        y0 = y1 + (y2-y1)/2
        z0 = sno
        print 'indexed centroid:', x0,y0,z0
        
        # derive radius of proportional samplling
        ra = np.sqrt(ha**2 + wa**2)/2
        rp = np.sqrt((hp-1)**2 + (wp-1)**2)/2
        rationp = int(round(rp/ra))
        skewp = wp/hp
        
        print("====================")
        print 'ra, rp, rationp, skewp:', ra,rp,rationp, skewp
        print("====================")
        
        # will sample np times in each directions
        # from centroid x0,y0 to a radius distance of  ra (hypothenuse of triangle)
        directions = [0, math.radians(45), math.radians(90), math.radians(135), math.radians(180), math.radians(225), math.radians(270), math.radians(315)]
        #[0, math.radians(45), math.radians(90), math.radians(135), math.radians(180), math.radians(225), math.radians(270), math.radians(315)]
        npatches = len(directions)*rationp+1
        print 'sampling %d times, for a total of %d patches.' % (rationp, npatches)
        xs=[]
        ys=[]
        allLpatches = []

        # show       
        fig, ax = plt.subplots(nrows=4, ncols=rationp*len(directions)+2, figsize=(4, 4))
        
        ## for the post-contract imgs
        # samples from image at xs, ys locations
        imgslicestime = [ pimg1[sno,:,:], pimg2[sno,:,:], pimg3[sno,:,:], pimg4[sno,:,:]]
        
        # init with a patch size centered in the orig patch
        xinit = (x0-ra/2)  
        yinit = (y0-ra/2)  
        centeredpatches = []
        minO = 6
        
        # extract patch inside the rectangular ROI
        for kimg in range(0,len(imgslicestime)):   
            patchk = imgslicestime[kimg][xinit:xinit+self.ptcsize[0],yinit:yinit+self.ptcsize[1]]  
            # add centered patch only once, the first patches
            centeredpatches = np.insert(centeredpatches, len(centeredpatches), patchk.flatten())
            
        # append extracted patch of size npratio * 8directions (30x30x4) and append 
        allLpatches.append(centeredpatches)        
      
        kp=2  
        for k in range(1,rationp+1):
            # per each k repeat a sample per each direction
            print 'sampling %d times......' % (k)
            for angle,j in zip(directions,range(len(directions))): 
                
                x = (x0-ra/2) - k*wa*rationp/minO*np.cos(angle)*wp/hp 
                y = (y0-ra/2) - k*ha*rationp/minO*np.sin(angle)*hp/wp
                
                # append and continue
                print 'sampling %f angle, (cos=%f, sin=%f), patch Origin: %s' % (int(math.degrees(angle)), np.cos(angle), np.sin(angle), str([x,y]) )
                xs.append(x)
                ys.append(y)
                
                # append all elemnet arrays
                allsubpatches = []
                
                for kimg in range(0,len(imgslicestime)):                    
                    print 'sampling image post-contrast time-point %d...... [-x %f,-y %f]' % ((kimg+1), k*wa*rationp/minO*np.cos(angle)*wp/hp , k*ha*rationp/minO*np.sin(angle)*hp/wp )
                    # extract patch inside the rectangular ROI
                    patchk = imgslicestime[kimg][x1:x2,y1:y2]                
                    #plot original patch
                    ax[kimg,0].imshow(patchk, cmap=plt.cm.gray) 
                    ax[kimg,0].set_title('imgslicestime,\ntime='+str(kimg+1), fontsize=10)
                    xmajor_ticks = [int(xi) for xi in [x1,x2]]
                    ymajor_ticks = [int(yi) for yi in [y1,y2]]
                    ax[kimg,0].set_xlabel(str(xmajor_ticks)) 
                    ax[kimg,0].set_ylabel(str(ymajor_ticks))
                                                    
                    #plot original patch
                    ax[kimg,1].imshow(patchk, cmap=plt.cm.gray) 
                    ax[kimg,1].set_title('centered patch', fontsize=10)
                    xmajor_ticks = [int(xi) for xi in [xinit,xinit+self.ptcsize[0]]]
                    ymajor_ticks = [int(yi) for yi in [yinit,yinit+self.ptcsize[1]]]
                    ax[kimg,1].set_xlabel(str(xmajor_ticks)) 
                    ax[kimg,1].set_ylabel(str(ymajor_ticks))
                                                   
                    # finally append all elemnet arrays
                    xp1 = int(x)
                    yp1 = int(y)
                    xp2 = int(x + self.ptcsize[0])
                    yp2 = int(y + self.ptcsize[1])
                
                    # subpatch
                    subpatchk = imgslicestime[kimg][xp1:xp2,yp1:yp2]
                    allsubpatches = np.insert(allsubpatches, len(allsubpatches), subpatchk.flatten())
                    
                    #plot subpatch
                    ax[kimg,kp].imshow(subpatchk, cmap=plt.cm.gray) 
                    ax[kimg,kp].set_title(str(k)+'_da_'+str(int(math.degrees(angle))), fontsize=10)
                    
                    xpmajor_ticks = [int(xi) for xi in [xp1,xp2]]
                    ypmajor_ticks = [int(yi) for yi in [yp1,yp2]]
                    
                    ax[kimg,kp].set_xlabel(str(xpmajor_ticks)) 
                    ax[kimg,kp].set_ylabel(str(ypmajor_ticks))
            
                # append extracted patch of size npratio * 8directions (30x30x4) and append 
                allLpatches.append(allsubpatches)
                kp += 1
                print len(allLpatches)
                
        # display
        # Fine-tune figure; make subplots close to each other and hide x ticks for
        # all but bottom plot.
        plt.setp([a.get_xticklabels() for a in fig.axes], visible=False)
        plt.setp([a.get_yticklabels() for a in fig.axes], visible=False)
        
        mng = plt.get_current_fig_manager()
        mng.window.showMaximized()
        fig.tight_layout()
        
        plt.savefig(snapshot_loc+os.sep+name+'.pdf')
        plt.show(block=False) 
        plt.close()
                
        return allLpatches

       
        