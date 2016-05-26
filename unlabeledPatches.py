# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 12:00:07 2016

Thresholding is the simplest way to segment objects from a background. 
If that background is relatively uniform, then you can use a global threshold value to binarize the image 
by pixel-intensity. If there’s large variation in the background intensity, however, 
adaptive thresholding (a.k.a. local or dynamic thresholding) may produce better results.

Here, we binarize an image using the threshold_adaptive function, which calculates thresholds in regions 
of size block_size surrounding each pixel (i.e. local neighborhoods). Each threshold value is the weighted mean 
of the local neighborhood minus an offset value.
source: http://scikit-image.org/docs/dev/auto_examples/plot_threshold_adaptive.html


@ Copyright (C) Cristina Gallego, University of Toronto, 2016
----------------------------------------------------------------------
"""
import os, os.path
import fnmatch
import sys
import numpy as np
import pandas as pd
import unicodedata

import matplotlib.pyplot as plt
from matplotlib.image import AxesImage
import SimpleITK as sitk
from dictionaries import mha_data_loc, data_loc
from query_localdatabase import *

from Imgpick import *
from patchify import *
from add_newpatches import *

import seaborn as sns
sns.set(color_codes=True)
from dictionaries import snapshot_loc
import six.moves.cPickle as pickle
import gzip

def ignore_exception(IgnoreException=Exception,DefaultVal=None):
    """ Decorator for ignoring exception from a function
    e.g.   @ignore_exception(DivideByZero)
    e.g.2. ignore_exception(DivideByZero)(Divide)(2/0)
    """
    def dec(function):
        def _dec(*args, **kwargs): 
            try:
                return function(*args, **kwargs)
            except IgnoreException:
                return DefaultVal
        return _dec
    return dec
    
sint = ignore_exception(ValueError)(int)
    
def filen_patt_match(mha_data_loc, patternstr):
    
    selfilenm = []
    for filenm in os.listdir(mha_data_loc):
        if fnmatch.fnmatch(filenm, patternstr):
            selfilenm = filenm      
        
    return selfilenm
    

if __name__ == '__main__':
    # Get Root folder ( the directory of the script being run)
    path_rootFolder = os.path.dirname(os.path.abspath(__file__))
    
    ## if fail at some lesion id, reload current pickled pathces and continue adding
    #allUpatches = []
    # or
    finitUpatches = gzip.open('allUpatches.pklz', 'rb')
    allUpatches = pickle.load(finitUpatches)
    finitUpatches.close()
            
    lesion_id = 571
    
    while ( lesion_id ) :            
        print("====================")
        print lesion_id
        print("====================")
        #############################
        ###### Query local databse
        #############################
        print " Querying local databse..."
        querylocal = Querylocal()
        dflesion = querylocal.queryby_lesionidwpatch(lesion_id)
        
        # process query, make sure there's a record of a lesion
        if(dflesion):
            lesion_record = [dflesion.Lesion_record.__dict__]
            dfinfo = pd.DataFrame.from_records(lesion_record)
            nmlesion_record = [dflesion.Nonmass_record.__dict__]
            dfnmlesion = pd.DataFrame.from_records(nmlesion_record)
            
            #############################
            ###### load mhas to memory
            #############################
            StudyID = str(int(dfinfo.iloc[0]['cad_pt_no_txt']))
            AccessionN = dfinfo.iloc[0]['exam_a_number_txt']
            preDynSeries = dfnmlesion.iloc[0]['DynSeries_id']
            print(dfinfo.iloc[0]['comment_txt'])
            print(dfinfo.iloc[0]['original_report'])
            
            try:
                s_img_no = int(dfnmlesion.iloc[0]['start_image_no_int'] )  
            except:
                s_img_no=None
            try:
                f_img_no = sint(dfnmlesion.iloc[0]['finish_image_no_int'] )
            except:
                f_img_no=None
                
            print("====================")
            print('s_img_no: ', s_img_no)
            print('f_img_no: ', f_img_no)
            nmenh_dist = dfnmlesion.iloc[0]['mri_nonmass_dist_int']
            nmenh_intenh = dfnmlesion.iloc[0]['mri_nonmass_int_enh_int']
            print('nmenh_dist: ', nmenh_dist)
            print('nmenh_intenh: ', nmenh_intenh)
            print("====================")
    
    
            ## or read from prior patch
            patch_record = pd.DataFrame.from_records([dflesion.lesion_patch.__dict__])
            c = patch_record.iloc[0]['centroid']
            centroid = c[c.find("(")+1:c.find(")")].split(',')
            sno=int(centroid[2])
            print('centroid patch: ', centroid)
            
            # define patch sample height and width: this is base on analysis of patch size distributions
            ha = 30
            wa = 30
            print("====================")
                    
            
            if(centroid): #sno
                # Save entire array as patch_set, patch ind to db (labelled patches of entire 1200 datasets)
                # if patch size is a x b x c, patch_set will be flatten array concatenated as axbxc, 
                # to reshape flatten array as image patch a.reshape(c,b,a)
                
                ## ====================                
                ## 1) read Image (right now needs to use filename pattern matching)
                ## ====================
                strtomatch = StudyID+'_'+AccessionN+'_'+str(int(preDynSeries)+1)
                filen_Vol1 = filen_patt_match(mha_data_loc, strtomatch+'*.mha')
                Vol1 = sitk.ReadImage(mha_data_loc+os.sep+filen_Vol1)
                strtomatch = StudyID+'_'+AccessionN+'_'+str(int(preDynSeries)+2)
                filen_Vol2 = filen_patt_match(mha_data_loc, strtomatch+'*.mha')
                Vol2 = sitk.ReadImage(mha_data_loc+os.sep+filen_Vol2)
                strtomatch = StudyID+'_'+AccessionN+'_'+str(int(preDynSeries)+3)
                filen_Vol3 = filen_patt_match(mha_data_loc, strtomatch+'*.mha')
                Vol3 = sitk.ReadImage(mha_data_loc+os.sep+filen_Vol3)
                strtomatch = StudyID+'_'+AccessionN+'_'+str(int(preDynSeries)+4)
                filen_Vol4 = filen_patt_match(mha_data_loc, strtomatch+'*.mha')
                Vol4 = sitk.ReadImage(mha_data_loc+os.sep+filen_Vol4)
        
                # reformat Vol slices as Float32
                Vol1 = sitk.Cast(Vol1,sitk.sitkFloat32)
                Vol2 = sitk.Cast(Vol2,sitk.sitkFloat32)
                Vol3 = sitk.Cast(Vol3,sitk.sitkFloat32)
                Vol4 = sitk.Cast(Vol4,sitk.sitkFloat32)
                
                # Get number of slices
                nslices = Vol1.GetDepth()
                nx = Vol1.GetWidth()
                ny = Vol1.GetHeight()
                
                if(f_img_no):
                    rndslices = [x for x in xrange(nslices) if x not in range(sno-1,f_img_no)]
                else:
                    rndslices = [x for x in xrange(nslices) if x not in [sno-2,sno-1,sno]]
                    
                selectedsl = random.choice(rndslices, size=5)
                # make sure to account for patch size so pixel init does not produce a patch outside image
                selectednx = random.choice(range(0,nx-ha), size=10)
                selectedny = random.choice(range(0,ny-wa), size=10)
                
                ## ====================                
                ## 2) Sample from random locations in the image and selected slices
                ## For 10 random slices, sample 50 random locations for a total of 500 unlabel patches per case
                ## ====================              
                ## Sample patch from Image
    #                3.1) convert sitk to np from int16 to uint8
    #                3.2) contrast stretch each image to include all intensities within 2nd and 98th percentils --> mapped to [0-255]
    #                # refs; http://scikit-image.org/docs/dev/auto_examples/plot_equalize.html; http://homepages.inf.ed.ac.uk/rbf/HIPR2/stretch.htm
    #                3.3) extract patches in x1, y1, x2, y2
                ## ====================
                extrP = Patchify()
                
                ## 3.2) contrast stretch
                pimg1 = extrP.contrast_stretch( sitk.GetArrayFromImage( Vol1) )
                pimg2 = extrP.contrast_stretch( sitk.GetArrayFromImage( Vol2) )
                pimg3 = extrP.contrast_stretch( sitk.GetArrayFromImage( Vol3) )
                pimg4 = extrP.contrast_stretch( sitk.GetArrayFromImage( Vol4) )
                
                patchname = str(lesion_id)+'_'+StudyID+'_'+AccessionN+'_'+preDynSeries+'_unlabeled'
                
                ## 3.3) extract patches
                Upatches = extrP.extractUnlabeledPatches( pimg1, pimg2, pimg3, pimg4, 
                                                      selectedsl, selectednx, selectedny, ha, wa, patchname)   
    
                
                ## ====================                
                ## 3) Append allUpatches Summarize and compare statistics
                ## ====================  
                # Append for all cases
                allUpatches.extend(Upatches)
                
                if( lesion_id % 2 == 0 ):
                    dfUpatches = [pd.DataFrame(pd.Series(up).describe()).transpose() for up in allUpatches]
                    desUpatches = pd.concat(dfUpatches)
                
                    sns.pairplot(desUpatches.ix[:,1:8]).savefig(snapshot_loc+os.sep+patchname+"_stats.pdf")
    
                    ## ====================                
                    ## 4) and save to file 
                    ## ====================            
                    # save to file 
                    fUpatches = gzip.open('allUpatches.pklz', 'wb')
                    pickle.dump(allUpatches, fUpatches, protocol=pickle.HIGHEST_PROTOCOL)
                    fUpatches.close()
    
    
                # TODO: sample N patches from other slices in the img (unl_patch), unl_patch cannot share pixels with patch ind
                # will have to implement Adaptive thesholding http://scikit-image.org/docs/dev/auto_examples/plot_threshold_adaptive.html    
                # Apply threshold. 
    #            block_size = 10
    #            bw = threshold_adaptive(imgArray[k], block_size)
    #            a[2*k+1].imshow(bw, cmap=plt.cm.gray) 
    #            
                # Histogram. 
                #values, bins = np.histogram(imgArray[k], bins=np.arange(256))
                #ax.plot(bins[:-1], values, lw=2, c='k')
  
            
        #############################
        ## continue to next case
        lesion_id += 1

        
        if lesion_id > 576:
            lesion_id = []
        








