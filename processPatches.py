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
import annot
import processDicoms


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
    lesion_id = 570
    
    while ( lesion_id ) :            
        #############################
        ###### Query local databse
        #############################
        print " Querying local databse..."
        querylocal = Querylocal()
        dflesion = querylocal.queryby_lesionid(lesion_id)
        
        # process query
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

        # Deal with no img info  
        # case: 1 1
        if s_img_no and f_img_no:  
            sno = s_img_no
        # case: 1 0
        if s_img_no and not f_img_no:
            sno = s_img_no
        # case: 0 1
        if not s_img_no and f_img_no:
            sno = f_img_no
        # case: 0 0
        if not s_img_no and not f_img_no:
            sno = None
            Sliceno = raw_input("Enter Slice no or 'x' to skip case:")
            if(Sliceno != 'x'):
                sno = int(Sliceno)
            ## or read from prior patch
#            patch_record = pd.DataFrame.from_records([dflesion.lesion_patch.__dict__])
#            c = patch_record.iloc[0]['centroid']
#            centroid = c[c.find("(")+1:c.find(")")].split(',')
#            sno=int(centroid[2])
#            print('centroid patch: ', centroid)
#            print("====================")
                
        #############################
        # Reveal annotations                                      
        #############################
        annotflag = False 
        annotationsfound = []
        study_loc = data_loc+os.sep+str(int(StudyID))+os.sep+AccessionN
        SeriesIDall =  [preDynSeries, str(int(preDynSeries)+1), str(int(preDynSeries)+2), str(int(preDynSeries)+3) ,str(int(preDynSeries)+4)]### for all processDicoms.get_immediate_subdirectories(study_loc)
        for iSer in SeriesIDall:
            exam_loc = data_loc+os.sep+str(int(StudyID))+os.sep+AccessionN+os.sep+str(iSer)
            print "Path Series annotation inspection: %s" % exam_loc
            annotationsfound, annotflag = annot.list_ann(exam_loc, annotflag, annotationsfound)             
        
        print annotationsfound
        
        if(sno): #sno
            # Save entire array as patch_set, patch ind to db (labelled patches of entire 1200 datasets)
            # if patch size is a x b x c, patch_set will be flatten array concatenated as axbxc, 
            # to reshape flatten array as image patch a.reshape(c,b,a)
            
            ## ====================                
            ## 1) read Image (right now needs to use filename pattern matching)
            ## ====================
            strtomatch = StudyID+'_'+AccessionN+'_'+preDynSeries
            filen_preVol = filen_patt_match(mha_data_loc, strtomatch+'*.mha')
            preVol = sitk.ReadImage(mha_data_loc+os.sep+filen_preVol)
            
            # deal with the rest:
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
            preVol = sitk.Cast(preVol,sitk.sitkFloat32)
            Vol1 = sitk.Cast(Vol1,sitk.sitkFloat32)
            Vol2 = sitk.Cast(Vol2,sitk.sitkFloat32)
            Vol3 = sitk.Cast(Vol3,sitk.sitkFloat32)
            Vol4 = sitk.Cast(Vol4,sitk.sitkFloat32)
            
            # plot slice location over time
            slicestime = [ Vol1[:,:,sno-1], 
                           Vol2[:,:,sno-1],
                           Vol3[:,:,sno-1], 
                           Vol4[:,:,sno-1] ]  
            
            ## ====================                
            ## 2) select patch, load annotation if found
            ## ====================
            # Initialize display picker
            selcPatch = Imgpick()            
            annots_dict_list = selcPatch.extract_annot(annotationsfound)
              
            # display to pick pixels
            selcPatch.selectDiagonal(slicestime, sno, annots_dict_list)        
            print('diagonal: ', selcPatch.diag)
            
            # x1,y1 pixel coords from selcPatch.diag[0]
            # x2,y2 pixel coords from selcPatch.diag[1]
            diagx1, diagy1 = [round(x) for x in selcPatch.diag[0]]
            diagx2, diagy2 = [round(x) for x in selcPatch.diag[1]]
            
            # select diag xs
            if(diagx2 < diagx1):
                x1 = diagx2-10
                x2 = diagx1+10
            else:
                x1 = diagx1-10
                x2 = diagx2+10
            
            # similarly do for diag ys
            if(diagy2 < diagy1):
                y1 = diagy2-10
                y2 = diagy1+10
            else:
                y1 = diagy1-10
                y2 = diagy2+10
        
            
                
            ## ====================                
            ## 3) Sample patch from Image
#                3.1) convert sitk to np from int16 to uint8
#                3.2) contrast stretch each image to include all intensities within 2nd and 98th percentils --> mapped to [0-255]
#                # refs; http://scikit-image.org/docs/dev/auto_examples/plot_equalize.html; http://homepages.inf.ed.ac.uk/rbf/HIPR2/stretch.htm
#                3.2) extract patches in x1, y1, x2, y2
            ## ====================
            extrP = Patchify()
            
            pimg0 = extrP.contrast_stretch( sitk.GetArrayFromImage(  preVol[:,:,sno-1] ))
            pimg1 = extrP.contrast_stretch( sitk.GetArrayFromImage(  Vol1[:,:,sno-1] ))
            pimg2 = extrP.contrast_stretch( sitk.GetArrayFromImage(  Vol2[:,:,sno-1] ))
            pimg3 = extrP.contrast_stretch( sitk.GetArrayFromImage(  Vol3[:,:,sno-1] ))
            pimg4 = extrP.contrast_stretch( sitk.GetArrayFromImage(  Vol4[:,:,sno-1] ))
            
            # extract patches   
            imgslicestime = [ pimg0, pimg1, pimg2, pimg3, pimg4]
            
            patchname = str(lesion_id)+'_'+StudyID+'_'+AccessionN+'_'+preDynSeries
            allpatches = extrP.extractPatches(imgslicestime, x1, y1, x2, y2, patchname)             
            print("====================")
            print('ROI x1 to x2: ', [x1,x2])
            print('ROI y1 to y2: ', [y1,y2])
            print('img_size: ',  extrP.img_size)
            print('ptcsize: ', extrP.ptcsize)
            print('resulting ptcsize = 5 x ptcsize[1] x ptcsize[0]: ', allpatches.shape)
            print("====================")
            
            # calculate patch centroid
            centroid = [x2-x1/2, y2-y1/2, sno]
            print(allpatche)
            ## ====================                
            ## 3) Send patch to db
            record = AddNewPatches()
            record.patch_2DB(lesion_id, str(tuple(extrP.img_size)), str(tuple(extrP.ptcsize)), str(tuple(centroid)), str(tuple(selcPatch.diag[0])), str(tuple(selcPatch.diag[1])), str(tuple(allpatches)) )
            
            # add annotation if any
            if(selcPatch.annot_ax):
                for annot in selcPatch.sel_annots:
                    record.gtannot_2DB(lesion_id, annot['AccessionNumber'], annot['SeriesDate'], annot['SeriesNumber'], annot['SliceLocation'], annot['SeriesDescription'],
                                            annot['PatientID'], annot['StudyID'], annot['SeriesInstanceUID'],
                                            annot['note'], annot['xi'], annot['yi'], annot['xf'], annot['yf'])
                 
                 
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
        print lesion_id
        
        if lesion_id > 1215:
            lesion_id = []
        








