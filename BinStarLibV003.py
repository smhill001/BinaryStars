# -*- coding: utf-8 -*-
"""
Created on Sat Dec 13 14:18:23 2014

@author: Steven

This library contains the five functions and classes needed by the main program,
StevesBinaryPlotterV03.py. The program generates plots of predicted and observed
orbital positions of binary stars. The functions and classes are listed below:

    Class plot_params
    Class ancillary_params
    Class StarOrbitObservations
    Class StarOrbitPrediction
    Class VisualBinaryPrediction
    Function DrawPlotV0
    Class StarRainbow
    
A GOOD MOVE TO PRACTICE INHERITENCE WOULD BE TO CREATE A STAR ORBIT CLASS
FROM WHICH BOTH THE PREDICTIONS AND OBSERVATIONS INHERIT.

Last updated 2014-12-16 SMH

"""

class plot_params:
    def __init__(self,StarIdentifierDD,View):
        
        self.T0=0.
        self.T1=0.
        self.dTF=0.
        self.dTC=0.
        self.X0=0.
        self.X1=0.
        self.NX=0.
        self.Y0=0.
        self.Y1=0.
        self.NY=0.
        self.ANSCL=0.

        CfgFile=open('StarPlotConfig.txt','r')
        CfgLines=CfgFile.readlines()
        CfgFile.close()
        nrecords=len(CfgLines)
        #print CfgLines

        for recordindex in range(1,nrecords):
            fields=CfgLines[recordindex].split(',')
            #print fields[0], fields[1]
            if fields[0] == StarIdentifierDD:
                print "In first if, fields[1]",fields[1],fields[0]
                if fields[1] == View:
                    print "In 2nd if, fields[1]",fields[1]
                    self.T0=float(fields[2])
                    self.T1=float(fields[3])
                    self.dTF=float(fields[4])
                    self.dTC=float(fields[5])
                    self.X0=float(fields[6])
                    self.X1=float(fields[7])
                    self.NX=float(fields[8])
                    self.Y0=float(fields[9])
                    self.Y1=float(fields[10])
                    self.NY=float(fields[11])
                    self.ANSCL=float(fields[12])
                    print self.dTF
                    
    def TimeGridFine(self):
        import numpy as np        
        print "In TimeGridFine", self.T0,self.T1,self.dTF
        return np.arange(self.T0, self.T1, self.dTF)
        
    def TimeGridCoarse(self):
        import numpy as np        
        return np.arange(self.T0, self.T1, self.dTC)

class ancillary_params:
    def __init__(self,StarIdentifierDD):
        
        self.CommonName=''
        self.Catalog=''
        self.DIST=0.

        CfgFile=open('StarAncillaryData.txt','r')
        CfgLines=CfgFile.readlines()
        CfgFile.close()
        nrecords=len(CfgLines)
        #print CfgLines

        for recordindex in range(1,nrecords):
            fields=CfgLines[recordindex].split(',')
            print fields[0], fields[1]
            if fields[0] == StarIdentifierDD:
                print "In first if, fields[1]",fields[1],StarIdentifierDD
                self.CommonName=str(fields[1])
                self.Catalog=str(fields[2])
                self.DIST=float(fields[3])
                    
class StarOrbitObservations:
    def __init__(self,StarIdentifierDD):
        import numpy as np
        ObsFile=open('Observations_SMH.txt','r')
        ObsLines=ObsFile.readlines()
        ObsFile.close()
        nrecords=len(ObsLines)
        self.TObs=[0.]
        self.XObs=[0.]
        self.YObs=[0.]
        self.PAObs=[0.]
        self.SepObs=[0.]
        self.PAErr=[0.]
        self.SepErr=[0.]
        self.NObs=0
        self.FirstTime=True

        for recordindex in range(1,nrecords):
            fields=ObsLines[recordindex].split(',')
            if fields[0] == StarIdentifierDD:
                #print fields
                if self.FirstTime:
                    self.TObs[0]=float(fields[1])
                    self.PAObs[0]=float(fields[2])
                    self.SepObs[0]=float(fields[3])
                    self.PAErr[0]=float(fields[4])
                    self.SepErr[0]=float(fields[5])
                    self.XObs[0]=self.SepObs[0]*np.cos((self.PAObs[0]+90.)*np.pi/180.)
                    self.YObs[0]=self.SepObs[0]*np.sin((self.PAObs[0]+90.)*np.pi/180.)
                    self.FirstTime=False
                    self.NObs=1
                else:
                    self.TObs.extend([float(fields[1])])
                    self.PAObs.extend([float(fields[2])])
                    self.SepObs.extend([float(fields[3])])
                    self.PAErr.extend([float(fields[4])])
                    self.SepErr.extend([float(fields[5])])
                    self.NObs=self.NObs+1                    
                    self.XObs.extend([self.SepObs[self.NObs-1]*np.cos((self.PAObs[self.NObs-1]+90.)*np.pi/180.)])
                    self.YObs.extend([self.SepObs[self.NObs-1]*np.sin((self.PAObs[self.NObs-1]+90.)*np.pi/180.)])
                #print self.XObs

    def ObsStats(self):
        import numpy as np
        # http://docs.scipy.org/doc/numpy/reference/generated/numpy.polyfit.html
        if self.NObs > 2:        
            SepCoefs=np.polyfit(self.TObs, self.SepObs, 1, \
                        w=1./(np.array(self.SepErr)),cov=True)
            SepFit=SepCoefs[0][0]*np.array(self.TObs)+SepCoefs[0][1]
            
            PACoefs=np.polyfit(self.TObs, self.PAObs, 1, \
                w=1./(np.array(self.PAErr)),cov=True)
            PAFit=PACoefs[0][1]+PACoefs[0][0]*np.array(self.TObs)
        else:
            SepCoefs=np.polyfit(self.TObs, self.SepObs, 1, \
                        w=1./(np.array(self.SepErr)),cov=False)
            SepFit=SepCoefs[0]*np.array(self.TObs)+SepCoefs[1]
            
            PACoefs=np.polyfit(self.TObs, self.PAObs, 1, \
                w=1./(np.array(self.PAErr)),cov=False)
            PAFit=PACoefs[1]+PACoefs[0]*np.array(self.TObs)
       
        
        sgn=[-1.,1.]
        index=0
        BoundSep=np.zeros((2,2))
        BoundPA=np.zeros((2,2))
        if self.NObs > 2:        
            for index in range(0,2):  
                dSint=SepCoefs[0][1]+sgn[index]*np.sqrt(SepCoefs[1][1,1])
                dSslp=SepCoefs[0][0]-sgn[index]*np.sqrt(SepCoefs[1][0,0])
                BoundSep[0,index]= dSslp*np.array(self.TObs).min()+dSint
                BoundSep[1,1-index]= dSslp*np.array(self.TObs).max()+dSint
                
                dPint=PACoefs[0][1]+sgn[index]*np.sqrt(PACoefs[1][1,1])
                dPslp=PACoefs[0][0]-sgn[index]*np.sqrt(PACoefs[1][0,0])
                BoundPA[0,index]= dPslp*np.array(self.TObs).min()+dPint
                BoundPA[1,1-index]= dPslp*np.array(self.TObs).max()+dPint
            print BoundSep
        
        return SepFit,BoundSep,SepCoefs,PAFit,BoundPA,PACoefs
              
                    
class StarOrbitPrediction:
    """A PREDICTED ORBIT object for a single star
        Attributes:
            Orbital elements
        Methods:
            __init__: load orbital elements from the WD6 catalog file
            __CalcEphemV2__: computing ephemerides for a given time grid
            """
    def __init__(self,StarIdentifierDD):

        WD6File=open('BinaryOrbitCat6CSVDel.txt','r')
        WD6Lines=WD6File.readlines()
        WD6File.close()
        nrecords=len(WD6Lines)
        
        for recordindex in range(1,nrecords-1):
            fields=WD6Lines[recordindex].split(',')
            #print fields[3]
            if fields[3] == StarIdentifierDD:
                StarIndex=recordindex

        StarField=WD6Lines[StarIndex].split(',')     
        
        self.Period=float(StarField[9])
        self.T=float(StarField[19])
        self.eccentricity=float(StarField[22])
        self.sma=float(StarField[12])
        self.inclination=float(StarField[15])
        self.AscNode=float(StarField[17])
        self.LongofPerhelion=float(StarField[24])
        self.MeanMotion=360./self.Period
        
    def displayElements(self):
        print "Elements"
        print self.Period,self.T, \
              self.eccentricity, self.sma, \
              self.inclination,self.AscNode, \
              self.LongofPerhelion,self.MeanMotion

    def CalcEphemV2(self,time,type):
        import numpy as np
    #    import pyasl as pyasl
    
        # Instantiate the solver
    #    ks = pyasl.MarkleyKESolver()
    
        # Solves Kepler's Equation for a set
        # of mean anomaly and eccentricity.
        # Uses the algorithm presented by
        # Markley 1995.
    #    M = 0.75
    #    e = 0.3
    #    print "Eccentric anomaly: ", ks.getE(M, e)
    # COMPUTE EPHEMERIDES - EQUATIONS BASED ON EXCEL SPREADSHEET
    
        LongofPerhelionRad=self.LongofPerhelion*np.pi/180.
        inclinationRad=self.inclination*np.pi/180.
        meananomaly=(self.MeanMotion*(time-self.T)%360.)/180.*np.pi #modulo and in radians
        E0=meananomaly
        #print E0
        #Iterate on mean anomaly
        for i in range (6):
            E0=E0+(meananomaly+self.eccentricity*np.sin(E0)-E0)/(1.-self.eccentricity*np.cos(E0))
        
        E6=E0
        
        RadiusVector=self.sma*(1.-self.eccentricity*np.cos(E6)) #arcsec
        TanVover2=np.sqrt((1.+self.eccentricity)/(1.-self.eccentricity))*np.tan(E6/2.)
        V=np.arctan(TanVover2)*2.0
        #Note: arctan2 has parameters reversed from MS Excel, e.g., arctan2(y,x) rather than atan2(x,y)
        ThetaminusAscNode=np.arctan2(np.sin(V+LongofPerhelionRad)*np.cos(inclinationRad),
                                     np.cos(V+LongofPerhelionRad))
        ThetaminusAscNodeDeg=ThetaminusAscNode*180./np.pi
        ApparentPADeg=(ThetaminusAscNodeDeg+self.AscNode)%360.
        ApparentRhoArcsec=RadiusVector*np.cos(V+LongofPerhelionRad)/np.cos(ThetaminusAscNode)
        X=ApparentRhoArcsec*np.cos((ApparentPADeg+90.)*np.pi/180.)
        Y=ApparentRhoArcsec*np.sin((ApparentPADeg+90.)*np.pi/180.)
        
        if type == 'Polar':      
          Eph1=ApparentRhoArcsec
          Eph2=ApparentPADeg
        if type == 'Cartesian':
            Eph1=X
            Eph2=Y
        return Eph1,Eph2


class VisualBinaryPrediction:
    """A PREDICTED LOCATION object for a single star
        Attributes:
            Orbital elements
        Methods:
            __init__: load orbital elements from the WD6 catalog file
            __CalcEphemV2__: computing ephemerides for a given time grid
            """
    def __init__(self,StarIdentifierDD):

        #WD6File=open('WDSCatalog20150812CSVDel.txt','r')
        WD6File=open('wdsweb_summ20150824CSV.txt','r')
        WD6Lines=WD6File.readlines()
        WD6File.close()
        nrecords=len(WD6Lines)
        
        for recordindex in range(1,nrecords-1):
            fields=WD6Lines[recordindex].split(',')
            #print fields[3]
            if fields[1] == StarIdentifierDD:
                StarIndex=recordindex

        StarField=WD6Lines[StarIndex].split(',')     
        
        self.LastEpoch=float(StarField[3])
        self.LastPA=float(StarField[6])
        self.LastSep=float(StarField[8])
        
    def displayElements(self):
        print "Elements"
        print self.LastEpoch,self.LastPA,self.LastSep

    def CalcEphemV2(self,time,type):
        import numpy as np
    
        t=np.zeros(time.size)+1.
        X=self.LastSep*np.cos((self.LastPA+90.)*np.pi/180.)*t
        Y=self.LastSep*np.sin((self.LastPA+90.)*np.pi/180.)*t
        
        if type == 'Polar':      
          Eph1=self.LastSep*t
          Eph2=self.LastPA*t
        if type == 'Cartesian':
            Eph1=X
            Eph2=Y
        return Eph1,Eph2
    

    # -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 14:37:28 2014

@author: steven.hill
"""
def DrawOrbitPlotV0(StarIdentifierDD, time_fine,X_fine,Y_fine,time_coarse,X_coarse,Y_coarse,Observations, \
    X_Fit,Y_Fit,X_Fit_Bound,Y_Fit_Bound,SepCoefs,PACoefs,params,anc_params,NumObs,CallNumber):
    import numpy as np
    import matplotlib.pyplot as pl
    pl.subplot(2, 2, CallNumber,axisbg="white")
###############       
    if CallNumber <= 2:    
        # Plot coarse prediction, fine prediction, and primary star
        pl.scatter(X_coarse, Y_coarse, color='black', s=20.,linewidth=0.0, linestyle="-")
        pl.plot(X_fine,Y_fine,color='black',linewidth=0.5)
        pl.scatter(0,0,s=50.,linewidth=0.0,color='black')
        if CallNumber == 2:
            print "***Call 2***"
            pl.plot(X_Fit,Y_Fit,label='Fit to Obs.',c='r')
            pl.plot(X_Fit_Bound,Y_Fit_Bound,label='Bounds',c='0.7')
            print X_Fit,Y_Fit
        # Create color array for observations - THIS SHOULD BE A CLASS!!!!!
        # Loop through observations and plot each in a different color
        # using modulo stepping through the color array. Create a legend
        # if this is the Zoom (second) plot
    
        clrs=StarRainbow()
        for j in range(0,len(Observations.XObs)):
            clr=clrs.c1[j % 6,:]
            pl.scatter(Observations.XObs[j],Observations.YObs[j],color=clr, \
                marker='+',label=str(Observations.TObs[j]))
        if CallNumber == 2:
            pl.legend(loc=1,fontsize=5,scatterpoints=1,ncol=1)
        
        pl.xlim(params.X0,params.X1)
        pl.xticks(np.linspace(params.X0,params.X1,params.NX, endpoint=True)) # Could be a param method...
        pl.ylim(params.Y0,params.Y1)
        pl.yticks(np.linspace(params.Y0,params.Y1,params.NY, endpoint=True)) # Could be a param method...
    
        pl.tick_params(axis='both', which='major', labelsize=6)
        pl.grid()
        pl.xlabel("Right Ascension (arcsec)",fontsize=6)
        pl.ylabel("Declination (arcsec)",fontsize=6)
    
        # Annotate the plots by labeling with the year of each coarse, predicted
        # point.
        for i in range(0,time_coarse.size): 
            pl.annotate('%.0f' % time_coarse[i],xy=(X_coarse[i],Y_coarse[i]), \
                        xytext=(X_coarse[i]*params.ANSCL,Y_coarse[i]*params.ANSCL), \
                        xycoords='data', \
                        textcoords='data',color='black',size=6,
                        horizontalalignment='center',verticalalignment='center')  
    
        # Create secondary axes by cloning the original axes
        ax2=pl.twinx()
        ay2=pl.twiny()
    
        pl.xlim(params.X0*anc_params.DIST,params.X1*anc_params.DIST)
        pl.xticks(np.linspace(params.X0*anc_params.DIST,params.X1*anc_params.DIST, params.NX, endpoint=True), \
            np.round(np.linspace(params.X0*anc_params.DIST,params.X1*anc_params.DIST, params.NX, endpoint=True),0))
        pl.xticks()
        pl.xlabel("Right Ascension (AU)",fontsize=6)
        
        pl.tick_params(axis='both', which='major', labelsize=6)
        pl.tick_params()
        
        pl.ylim(params.Y0*anc_params.DIST,params.Y1*anc_params.DIST)
        pl.yticks(np.linspace(params.Y0*anc_params.DIST,params.Y1*anc_params.DIST, params.NY, endpoint=True), \
            np.round(np.linspace(params.Y0*anc_params.DIST,params.Y1*anc_params.DIST, params.NY, endpoint=True),0))
        pl.ylabel("Declination (AU)",fontsize=7)
        
        ax2.tick_params(axis='both', which='major', labelsize=6)
        ay2.tick_params(axis='both', which='major', labelsize=6)
    
        print "ax2=",ax2
        pl.suptitle(StarIdentifierDD+' - '+anc_params.CommonName,fontsize=8, fontweight='bold')
        #Plot origin lines
        pl.axhline(y=0,color='grey')
        pl.axvline(x=0,color='grey')
    
        # Save figure using 72 dots per inch
        # savefig("exercice_2.png", dpi=72)
        # Show result on screen
        #pl.show()
##################    
    if CallNumber == 3:
        clrs=StarRainbow()
        for j in range(0,len(Observations.XObs)):
            clr=clrs.c1[j % 6,:]
            pl.scatter(Observations.TObs[j],Observations.SepObs[j],color=clr, \
                marker='+')
        
        pl.plot(Observations.TObs,X_Fit,label='Fit to Obs.',c='r')
        pl.scatter(Observations.TObs,X_fine, color='black', s=20.,linewidth=0.0, linestyle="-",label='Predicted')
        pl.plot(Observations.TObs,X_fine,color='black',linewidth=0.5)
        pl.plot([np.array(Observations.TObs).min(),np.array(Observations.TObs).max()], \
            X_Fit_Bound[:,0],label='Lower Bound',c='0.7')
        pl.plot([np.array(Observations.TObs).min(),np.array(Observations.TObs).max()], \
            X_Fit_Bound[:,1],label='Upper Bound',c='0.7')
        
        
        pl.xlim(params.T0,params.T1)  #COULD HARDWIRE X-AXIS HERE, WON'T REALLY CHANGE
        print params.T1,params.T0,params.dTC
        xtks=(params.T1-params.T0)/params.dTC+1
        print xtks
        pl.xticks(np.linspace(params.T0,params.T1,xtks,endpoint=True)) # Could be a param method...
        pl.tick_params(axis='both', which='major', labelsize=6)
        pl.grid()
        
        pl.ylim(params.X0,params.X1) 
        pl.yticks(np.linspace(params.X0,params.X1,params.NX, endpoint=True)) # Could be a param method...
        
        pl.ylabel("Separation - Rho (arcsec)",fontsize=6)
        pl.legend(loc=1,fontsize=5,scatterpoints=1,ncol=1)
        
        
        print '***', (X_fine[Observations.NObs-1]-X_fine[0])
        print '***', (Observations.TObs[Observations.NObs-1]-Observations.TObs[0])
        
        PrdSepRate=(X_fine[Observations.NObs-1]-X_fine[0]) / (Observations.TObs[Observations.NObs-1]-Observations.TObs[0])
        TimeCornerForText=0.97*(params.T1-params.T0)+params.T0
        VertCornerForText=0.03*(params.X1-params.X0)+params.X0
        if NumObs >2:
            pl.text(TimeCornerForText,VertCornerForText,'O  ='+str(np.mean(Observations.SepObs))[0:5]+' asec\n' \
                            +'P  ='+str(np.mean(X_fine))[0:5]+' asec\n' \
                            +'O-P='+str(np.mean(Observations.SepObs)-np.mean(X_fine))[0:6]+' asec\n\n'\
                            +'dO  ='+str(SepCoefs[0][0])[0:5]+' asec/yr\n' \
                            +'dP  ='+str(PrdSepRate)[0:5]+' asec/yr\n' \
                            +'dO-dP='+str(SepCoefs[0][0]-PrdSepRate)[0:6]+' asec/yr'\
                            ,fontsize=5,verticalalignment='bottom',horizontalalignment='right', \
                            bbox=dict(facecolor='white', alpha=1))

    if CallNumber == 4:
        clrs=StarRainbow()
        for j in range(0,len(Observations.XObs)):
            clr=clrs.c1[j % 6,:]
            pl.scatter(Observations.TObs[j],Observations.PAObs[j],color=clr, \
                marker='+')
        
        pl.plot(Observations.TObs,Y_Fit,label='Fit to Obs.',c='r')
        pl.scatter(Observations.TObs,Y_fine, color='black', s=20.,linewidth=0.0, linestyle="-",label='Predicted')
        pl.plot(Observations.TObs,Y_fine,color='black',linewidth=0.5)
        pl.plot([np.array(Observations.TObs).min(),np.array(Observations.TObs).max()], \
            Y_Fit_Bound[:,0],label='Lower Bound',c='0.7')
        pl.plot([np.array(Observations.TObs).min(),np.array(Observations.TObs).max()], \
            Y_Fit_Bound[:,1],label='Upper Bound',c='0.7')
        
        
        pl.xlim(params.T0,params.T1)  #COULD HARDWIRE X-AXIS HERE, WON'T REALLY CHANGE
        xtks=(params.T1-params.T0)/params.dTC+1
        print xtks
        pl.xticks(np.linspace(params.T0,params.T1,xtks,endpoint=True)) # Could be a param method...
        pl.tick_params(axis='both', which='major', labelsize=6)
        pl.grid()
        
        pl.ylim(params.Y0,params.Y1) 
        pl.yticks(np.linspace(params.Y0,params.Y1,params.NY, endpoint=True)) # Could be a param method...
        
        pl.ylabel("Position Angle - Theta (degrees)",fontsize=6)
        pl.legend(loc=1,fontsize=5,scatterpoints=1,ncol=1)
        
        
        print '***', (Y_fine[Observations.NObs-1]-X_fine[0])
        print '***', (Observations.TObs[Observations.NObs-1]-Observations.TObs[0])
        
        PrdPARate=(Y_fine[Observations.NObs-1]-Y_fine[0]) / (Observations.TObs[Observations.NObs-1]-Observations.TObs[0])
        TimeCornerForText=0.97*(params.T1-params.T0)+params.T0
        VertCornerForText=0.03*(params.Y1-params.Y0)+params.Y0
        if NumObs >2:
            pl.text(TimeCornerForText,VertCornerForText,'O  ='+str(np.mean(Observations.PAObs))[0:5]+' deg\n' \
                            +'P  ='+str(np.mean(Y_fine))[0:5]+' deg\n' \
                            +'O-P='+str(np.mean(Observations.PAObs)-np.mean(Y_fine))[0:6]+' deg\n\n'\
                            +'dO  ='+str(PACoefs[0][0])[0:5]+' deg/yr\n' \
                            +'dP  ='+str(PrdPARate)[0:5]+' deg\yr\n' \
                            +'dO-dP='+str(PACoefs[0][0]-PrdPARate)[0:6]+' deg/yr'\
                            ,fontsize=5,verticalalignment='bottom',horizontalalignment='right', \
                            bbox=dict(facecolor='white', alpha=1))
       
    return 1
    
    # -*- coding: utf-8 -*-

                    
class StarRainbow:
    def __init__(self):
        import numpy as np
        self.c1=np.array([[0.7,0.,0.],[0.7,0.7,0.0],[0.,0.7,0.], \
                             [0.7,0.,0.7],[0.,0.7,0.7],[0.,0.,0.7]])
