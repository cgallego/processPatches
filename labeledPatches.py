# -*- coding: utf-8 -*-
"""
Created on Sat May 21 09:41:41 2016

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
#    allLpatches = []
#    allLabels = pd.DataFrame()
    # or
    finitLpatches = gzip.open('allLpatches.pklz', 'rb')
    allLpatches = pickle.load(finitLpatches)
    finitLpatches.close() 
    # for labesl
    finitLables = gzip.open('allLabels.pklz', 'rb')
    allLabels = pickle.load(finitLables)
    finitLables.close()
                    
    lesion_id = 31
    
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
            p = patch_record.iloc[0]['patch_size']
            patch_size = p[p.find("(")+1:p.find(")")].split(',')
            print('centroid patch: ', centroid)
            print('patch_size: ', patch_size) 
            p1 = patch_record.iloc[0]['patch_diag1'] 
            patch_diag1 = p1[p1.find("(")+1:p1.find(")")].split(',')
            p2 = patch_record.iloc[0]['patch_diag2']
            patch_diag2 = p2[p2.find("(")+1:p2.find(")")].split(',')
            print("====================")
                                
            if(centroid): #sno
                # Save entire array as patch_set, patch ind to db (labelled patches of entire 1200 datasets)
                # if patch size is a x b x c, patch_set will be flatten array concatenated as axbxc, 
                # to reshape flatten array as image patch a.reshape(c,b,a)
                # define patch sample height and width: this is base on analysis of patch size distributions
                ha = 30
                wa = 30
                
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
                
                ## ====================                
                ## 2) Sample from the centroid of the patch, patches of size [ha,wa] 
                ## Sampling scheme that scales with the size of the original path Hp,Wp
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
                
                patchname = str(lesion_id)+'_'+StudyID+'_'+AccessionN+'_'+preDynSeries+'_lesionpatches'
                
                ## 3.3) extract patches
                Lpatches = extrP.extractROIPatches( pimg1, pimg2, pimg3, pimg4, centroid, patch_size,
                                                   patch_diag1, patch_diag2, ha, wa, patchname)   
             
                ## ====================                
                ## 3) Query labels, append allpatches and Summarize and compare statistics
                ## ====================  
                # create labels for corresponding patches                                
                dfl = pd.DataFrame(pd.Series({'lesion_id': lesion_id, 'nmenh_dist': nmenh_dist, 'nmenh_intenh': nmenh_intenh})).transpose()
                # append
                dflabels = pd.DataFrame()
                dflabels = dflabels.append([dfl]*len(Lpatches), ignore_index=True)
                
                ## ====================  
                # APPEND: Note: only do once for both patches and lables
                # Append for all cases
                allLpatches.extend(Lpatches)
                allLabels = allLabels.append(dflabels)
                ## ====================                  
                
                if( lesion_id % 10 == 0 ):
                    dfLpatches = [pd.DataFrame(pd.Series(lesion).describe()).transpose() for lesion in allLpatches]
                    desLpatches = pd.concat(dfLpatches)
                
                    sns.pairplot(desLpatches.ix[:,1:8]).savefig(snapshot_loc+os.sep+patchname+"_stats.pdf")
    
                    ## ====================                
                    ## 4) and save to files, patches and labels will have same dim 1
                    ## ====================            
                    # save to file allLpatches.pklz for patches
                    fLpatches = gzip.open('allLpatches.pklz', 'wb')
                    pickle.dump(allLpatches, fLpatches, protocol=pickle.HIGHEST_PROTOCOL)
                    fLpatches.close()
                    # save to file allLabels.pklz for labesl
                    fLables = gzip.open('allLabels.pklz', 'wb')
                    pickle.dump(allLabels, fLables, protocol=pickle.HIGHEST_PROTOCOL)
                    fLables.close()
    
        #############################
        ## continue to next case
        lesion_id += 1

        
        if lesion_id > 576:
            lesion_id = []
        








