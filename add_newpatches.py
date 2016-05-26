# -*- coding: utf-8 -*-
"""
USAGE: 
=============
from add_newPatches import *
record = AddNewPatches()
record.patch_2DB(image, ...)

Class Methods:
=============
dicomTransform(image, image_pos_pat, image_ori_pat)
addSegment(lesion3D)
subImage(Images2Sub, timep)                  
visualize(images, image_pos_pat, image_ori_pat, sub, postS, interact)

Class Instance Attributes:
===============

Created on Mon May 09 16:10:17 2016

@ author (C) Cristina Gallego, University of Toronto
--------------------------------------------------------------------
 """
import os, os.path
import sys
import string
from sys import argv, stderr, exit
from numpy import *
import pandas as pd

from sqlalchemy.orm import sessionmaker
from base import localengine
import mylocaldatabase

class AddNewPatches(object):
    """
    USAGE:
    =============
    record = AddNewPatches()
    """
    def __init__(self): 
        """ initialize database session """           
        #  create a top level Session configuration which can then be used throughout
        # Create the Session
        self.Session = sessionmaker()
        self.Session.configure(bind=localengine)  # once engine is available
        
    def __call__(self):       
        """ Turn Class into a callable object """
        AddNewPatches() 

    def patch_2DB(self, lesion_id, img_size, ptcsize, centroid, diag1, diag2, allpatches):      
        
        self.session = self.Session() #instantiate a Session
        # Send to database lesion info
        patch_info = mylocaldatabase.lesion_patch(lesion_id, img_size, ptcsize, centroid, diag1, diag2, allpatches)
        self.session.add(patch_info)
        
        # Finally send records to database
        try:
            self.session.commit()
        except:
            self.session.rollback()
            raise
        finally:
            self.session.close()
            
        return
        
    def gtannot_2DB(self, lesion_id, AccessionNumber, SeriesDate, SeriesNumber, SliceLocation, SeriesDescription,
                 PatientID, StudyID, SeriesInstanceUID, note, xi_coord, yi_coord, xf_coord, yf_coord):      
        
        self.session = self.Session() #instantiate a Session
        # Send to database lesion info
        annot_info = mylocaldatabase.Annot_record(lesion_id, AccessionNumber, SeriesDate, SeriesNumber, SliceLocation, SeriesDescription,
                 PatientID, StudyID, SeriesInstanceUID, note, xi_coord, yi_coord, xf_coord, yf_coord)
                 
        self.session.add(annot_info)
        
        # Finally send records to database
        try:
            self.session.commit()
        except:
            self.session.rollback()
            raise
        finally:
            self.session.close()
            
        return
