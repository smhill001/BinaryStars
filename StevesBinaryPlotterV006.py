# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 07:19:54 2014

@author: s.hill
"""
### FUNCTION:
# 1) Read, parse and search WDS for desired star
# 2) Take user inputs/preferences, e.g., axis scale and intervals
# 3) Computes ephemerides based on WDS elements
# 4) Reads and parses single star observations file (Needs to be object)
#Needs to read and parse variable length files of single star observations
#Needs to handle headers and multiple star observation input files
# X) Creates two orbital plots (WIDE and ZOOM) with annotations and legend

#Needs to annotate with meta-data, e.g., orbital period, other identifiers?

# SHOULD CREATE A PARENT CLASS FOR STELLAR ORBIT AND LET OBS AND PRD INHERIT
# BOTH SUBCLASSES SHOULD BE ABLE TO PROVIDE POSITIONAL INFORMATION
# AS X-Y PAIRS AND THETA-RHO PAIRS


#NEED TO FIGURE OUT HOW TO AVOID SCREEN SIZE LIMITS ON FIGURE BOUNDS
#Needs to provide nice, parseable output text

import BinStarLibV003 as BSL
import matplotlib.pyplot as pl
import numpy as np


# Create Observations and Predicted Orbit Objects for Star
#StarIdentifierDD='STF1110AB' # Castor
#StarIdentifierDD='STF2272AB' # 70 Oph
#StarIdentifierDD='STF2055AB' # Lam Oph
#StarIdentifierDD='STF2382AB' # Eps 1 Lyr
#StarIdentifierDD='STF2383Cc-D' # Eps 2 Lyr
#StarIdentifierDD='BU  648AB'
#StarIdentifierDD='STF2084' #Zet Her
#StarIdentifierDD='STF  60AB' #Eta Cas
#StarIdentifierDD='STF2130AB' #Mu Dra
#StarIdentifierDD='STF 668A BC' #Rigel
#StarIdentifierDD='STF 774Aa B' #Alnitak
StarIdentifierDD='DA    5AB' #Eta Ori

CastorAncParams=BSL.ancillary_params(StarIdentifierDD)
print ''
print 'CommonName = ',CastorAncParams.CommonName
print ''
if CastorAncParams.Catalog=='ORB6':
    CastorPrd=BSL.StarOrbitPrediction(StarIdentifierDD) ###This is where I need different catlogs
elif CastorAncParams.Catalog=='WDS':
    CastorPrd=BSL.VisualBinaryPrediction(StarIdentifierDD) ###This is where I need different catlogs

CastorObs=BSL.StarOrbitObservations(StarIdentifierDD)

CastorWide=BSL.plot_params(StarIdentifierDD,'Wide')
CastorZoom=BSL.plot_params(StarIdentifierDD,'Zoom')
CastorPASep=BSL.plot_params(StarIdentifierDD,'PASep')

# COMPUTE FITS AND BOUNDARIES
#if CastorAncParams.Catalog=='ORB6': 
Sep_PrdObs,PA_PrdObs=CastorPrd.CalcEphemV2(np.array(CastorObs.TObs),'Polar')

SepFit,BoundSep,SepCoefs,PAFit,BoundPA,PACoefs = CastorObs.ObsStats()

X_Fit=SepFit*np.cos((PAFit+90.)*np.pi/180.)
Y_Fit=SepFit*np.sin((PAFit+90.)*np.pi/180.)
X_Fit_Bound=BoundSep*np.cos((BoundPA+90.)*np.pi/180.)
Y_Fit_Bound=BoundSep*np.sin((BoundPA+90.)*np.pi/180.)
print X_Fit,Y_Fit
#####################

pl.figure(figsize=(6.5, 6.5), dpi=150,facecolor="white")

# WIDE PLOT

print CastorWide.T0,CastorWide.T1,CastorWide.dTF

X_fine,Y_fine=CastorPrd.CalcEphemV2(CastorWide.TimeGridFine(),'Cartesian')
X_coarse,Y_coarse=CastorPrd.CalcEphemV2(CastorWide.TimeGridCoarse(),'Cartesian')
BSL.DrawOrbitPlotV0(StarIdentifierDD,CastorWide.TimeGridFine(),X_fine,Y_fine, \
                    CastorWide.TimeGridCoarse(),X_coarse,Y_coarse,\
                    CastorObs,X_Fit,Y_Fit,X_Fit_Bound,Y_Fit_Bound,SepCoefs,PACoefs,CastorWide,CastorAncParams,CastorObs.NObs,1)

# ZOOM PLOT
X_fine,Y_fine=CastorPrd.CalcEphemV2(CastorZoom.TimeGridFine(),'Cartesian')
X_coarse,Y_coarse=CastorPrd.CalcEphemV2(CastorZoom.TimeGridCoarse(),'Cartesian')
BSL.DrawOrbitPlotV0(StarIdentifierDD,CastorZoom.TimeGridFine(),X_fine,Y_fine, \
                    CastorZoom.TimeGridCoarse(),X_coarse,Y_coarse, \
                    CastorObs,X_Fit,Y_Fit,X_Fit_Bound,Y_Fit_Bound,SepCoefs,PACoefs,CastorZoom,CastorAncParams,CastorObs.NObs,2)

BSL.DrawOrbitPlotV0(StarIdentifierDD,CastorObs.TObs,Sep_PrdObs,PA_PrdObs, \
                    CastorObs.TObs,Sep_PrdObs,PA_PrdObs, \
                    CastorObs,SepFit,PAFit,BoundSep,BoundPA,SepCoefs,PACoefs,CastorPASep,CastorAncParams,CastorObs.NObs,3)
                    
BSL.DrawOrbitPlotV0(StarIdentifierDD,CastorObs.TObs,Sep_PrdObs,PA_PrdObs, \
                    CastorObs.TObs,Sep_PrdObs,PA_PrdObs, \
                    CastorObs,SepFit,PAFit,BoundSep,BoundPA,SepCoefs,PACoefs,CastorPASep,CastorAncParams,CastorObs.NObs,4)

pl.savefig(StarIdentifierDD+'-'+CastorAncParams.CommonName+'.png', format='png',dpi=300)