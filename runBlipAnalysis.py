##############################################################################
#*######################                              ######################*#
#######################     Ryan's Blip Analysis Here  #######################
#*######################                              ######################*#
#*##########################################################################*#
#Blip Analysis Code Used to examine low energy signals in our LArTPC Detector#

import math
import itertools
from math import exp, pi
from ROOT import TFile, TTree

import numpy as np
import matplotlib.pylab as plt
import matplotlib.pyplot as plot
import matplotlib.axes as axes
from matplotlib.colors import LogNorm

import pandas as pd
import scipy as sp
from scipy.optimize import curve_fit
from scipy.stats import norm
import matplotlib.mlab as mlab
from pylab import rcParams
import seaborn as sns
sns.set()

from root_numpy import root2array, tree2array,testdata
#from root_pandas import read_root
from glob import glob
from matplotlib.colors import LogNorm
from scipy.interpolate import interp1d
from scipy.misc import factorial
plt.rcParams.update({'font.size': 20})

#*##########################################################################*#
#################  Begin Importing Data, Add a few functions  #################
#*##########################################################################*#
from blipFunctions import *

#Import the file below - split into event branches. Import from the root file to a pandas DF
blip_file1='/data/users/rdorrill/BlipAna_20220712_RadonData_Phase1_BlipTree.root'
fullFileTree='/blipana/anatree'
newTree='/blipanaTrkMask/anatree'
EVT_BRANCH = '/blipana/anatree/event'
splitTree_ANA = 'anatree'
blipTree= 'bliptree'
newFileTree='/blipana/bliptree'
df_blip_file1 = pd.DataFrame( root2array(blip_file1,blipTree) )
print("Length of data frame and a short sample of events: ")
print((len(df_blip_file2)))
df_blip_file2.head()

#QOL Functions to plot the initial location of blips and energy

#*##########################################################################*#
#################  Cuts on the Data in the df  #################
#*##########################################################################*#

#G10_region is in our hotspot regions. Rest_det is everything outside those areas. Complimentary is an equal volume to the G10 points, outside them

G10_region='((blip_y < -90.0 or blip_y > 90.0) and ((blip_z>10 and blip_z<50) or (blip_z>90 and blip_z<130) or (blip_z>210 and blip_z<250) or (blip_z>320 and blip_z<360) or (blip_z>440 and blip_z<470) or (blip_z>560 and blip_z<600) or (blip_z>660 and blip_z<700) or (blip_z>780 and blip_z<810) or (blip_z>890 and blip_z<930) or (blip_z>980 and blip_z<1100) )) or ((blip_y < 60.0 and blip_y > 40.0 ) and ((blip_z>0 and blip_z<40) or (blip_z>1000 and blip_z<1200))) or ((blip_y < -40.0 and blip_y > -60.0 ) and ((blip_z>0 and blip_z<40)or (blip_z>1000 and blip_z<1200))) '
rest_det_region='(blip_y > -90.0 and blip_y < 90.0) and (blip_z>20 and blip_z<1000)'
complementary_region=' ( (blip_y < -90.0 or blip_y > 90.0) and ( (blip_z>50 and blip_z<90) or (blip_z>150 and blip_z<195) or (blip_z>265 and blip_z<305) or (blip_z>375 and blip_z<420) or (blip_z>488 and blip_z<540) or (blip_z>610 and blip_z<645) or (blip_z>720 and blip_z<765) or (blip_z>830 and blip_z<875) or (blip_z>945 and blip_z<980) or (blip_z>1100 ) )) or ( (blip_y < 85.0 and blip_y > 60.0) and ((blip_z>10 and blip_z<50) or (blip_z>1000 and blip_z<1036.5 ) ))  or  ( (blip_y > -85.0 and blip_y < -60.0) and ((blip_z>10 and blip_z<50) or (blip_z>1000 and blip_z<1036.5 ) )  )'
#Volumes equivalent to the G10 spots, but in the center of the det, and near the outside:
G10_equivalent_Vol_Center='(blip_y > -27.0 and blip_y < 27.0) and (blip_z>150 and blip_z<592.8) '
G10_equivalent_Vol_NearG10='(blip_y > -90.0 and blip_y < -36.0) and (blip_z>150 and blip_z<593.8) '

df_blip_file1_3planes=df_blip_file1.query('(blip_nplanes==3) and (blip_sigmayz < 2) and (blip_incylinder == False)')
df_blip_file1_3planes_g10=df_blip_file1_3planes.query(G10_region)
df_blip_file1_3planes_complementary=df_blip_file1_3planes.query(complementary_region)
df_blip_file1_3planes_equivalentVol=df_blip_file1_3planes.query(G10_equivalent_Vol_NearG10)

#*##########################################################################*#
#################  Begin Main  #################
#*##########################################################################*#

#Example of Data before the cut, plotting in 2D by location, and the Energies
plotAllRegionBlips(df_blip_file1_3planes,df_blip_file1_3planes_g10,df_blip_file1_3planes_complementary,df_blip_file1_3planes_equivalentVol)
compareEnergiesByRegion(df_blip_file1_3planes_g10,df_blip_file1_3planes_complementary,df_blip_file1_3planes_equivalentVol)
plot3MevBlipLocations(df_blip_file1_3planes)
