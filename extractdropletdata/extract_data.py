#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
from scipy.optimize import curve_fit

class droplet:
    def __init__(self):
        self.__dropletcount = 0
        
        self.__dropletmax  = None
        self.__dropletbase = None
        self.__droplettime = None
        self.__dropletdataindex = None
        
        self.__fitdata = {}
        self.__fitkeys = []

    
    def add_droplet(self,maxval,baseval,time,idx = 0):
        if self.__dropletcount > 0:
            self.__dropletmax  = np.concatenate((self.__dropletmax ,np.array([maxval])  ))
            self.__dropletbase = np.concatenate((self.__dropletbase,np.array([baseval]) ))
            self.__droplettime = np.concatenate((self.__droplettime,np.array([time])    ))
            self.__dropletdataindex = np.concatenate((self.__dropletdataindex,np.array([idx]) ))
        else:
            self.__dropletmax  = np.array([maxval])
            self.__dropletbase = np.array([baseval])
            self.__droplettime = np.array([time])
            self.__dropletdataindex = np.array([idx])
        self.__dropletcount += 1
    
    
    def add_droplet_fitinfo(self,key,dropletID,values):
        c = len(values)
        if not key in self.__fitkeys:
            self.__fitdata[key] = np.zeros((self.__dropletcount,c))
            self.__fitkeys.append(key)
        if dropletID < self.__dropletcount:
            self.__fitdata[key][dropletID] = values
    
    def __getitem__(self,key):
        if key in self.__fitkeys:
            return self.__fitdata[key]
        elif key == "time":
            return self.__droplettime
        elif key == "max":
            return self.__dropletmax
        elif key == "base":
            return self.__dropletbase
        elif key == "index":
            return self.__dropletdataindex
        elif key == "count":
            return self.__dropletcount
    
    def __int__(self):
        return self.__dropletcount
    
    def get_max(self):
        return self.__dropletmax
    def get_times(self):
        return self.__droplettime
    def get_base(self):
        return self.__dropletbase
    def get_count(self):
        return self.__dropletcount
    def get_data(self):
        return self.__droplettime,self.__dropletmax,self.__dropletbase
                        
    def __iter__(self):
        for key in self.__fitkeys:
            yield self[key]


def gaussoffset(x,mean,stddev,sqrtmaxval,sqrtbaseval):
    return np.exp(-0.5*(x-mean)**2/stddev**2)*sqrtmaxval**2+sqrtbaseval**2

def stepoffset(x,t1,t2,sqrtmaxval,sqrtminval):
    r = np.ones(len(x))*sqrtminval**2
    i1 = ((x-t1)**2).argmin()
    i2 = ((x-t2)**2).argmin()
    r[i1:i2] = sqrtmaxval**2
    return r

def sigmoidoffset(x,mean,sqrtcrossover,sqrtsteepness,sqrtheight,sqrtbase):
    r = np.zeros(len(x))
    i1 = ((x-mean)**2).argmin()
    r[:i1] = sqrtheight**2 /(1.+np.exp(-(x[:i1]-mean+sqrtcrossover**2)/sqrtsteepness**2)) + sqrtbase**2
    r[i1:] = sqrtheight**2 /(1.+np.exp( (x[i1:]-mean-sqrtcrossover**2)/sqrtsteepness**2)) + sqrtbase**2
    return r


class timeseries:
    def __init__(self,filename,minmax_lowerthreshold = .15, minmax_upperthreshold = .3,maxfev=5000):
        self.__filename = filename
        try:
            self.__timeseriesdata = np.genfromtxt(filename)
        except:
            raise IOError
        
        self.__droplets = droplet()
        self.__datalength = len(self.__timeseriesdata[:,0])
        self.__minmax_lowerthreshold = minmax_lowerthreshold
        self.__minmax_upperhreshold = minmax_upperthreshold
        
        self.__maxfev = maxfev
    
    def get_droplet_data_minmax(self):
        inlower = True
        inupper = False
        enteredlower = 0
        for i in range(self.__datalength):
            if inlower and self.__timeseriesdata[i,1] > self.__minmax_upperhreshold:
                inupper = True
                inlower = False
                enteredupper = i
            if inupper and self.__timeseriesdata[i,1] < self.__minmax_lowerthreshold:
                inupper = False
                inlower = True
                
                baseidx = self.__timeseriesdata[enteredlower:i,1].argmin()+enteredlower
                baseval = self.__timeseriesdata[baseidx,1]
                maxidx  = self.__timeseriesdata[enteredlower:i,1].argmax()+enteredlower
                maxval  = self.__timeseriesdata[maxidx,1]
                time    = self.__timeseriesdata[maxidx,0]
                
                self.__droplets.add_droplet(maxval,baseval,time,maxidx)
                enteredlower = i
    
    def fit_step(self,startindex,finalindex):
        t = self.__timeseriesdata[startindex:finalindex,0]
        p = self.__timeseriesdata[startindex:finalindex,1]
        
        maxidx = p.argmax()
        minidx = p.argmin()
        
        p0 = np.array([t[maxidx]-0.03,t[maxidx]+0.03,np.sqrt(p[maxidx]),np.sqrt(p[minidx] if p[minidx]>0 else 0)])
        
        parameters,cov = curve_fit(stepoffset,t,p,p0=p0,maxfev=self.__maxfev)
        fitparam = np.array([0.5*(parameters[0]+parameters[1]),parameters[2]**2,parameters[3]**2,abs(parameters[1]-parameters[0])])
        return fitparam
        
    
    def fit_gaussian(self,startindex,finalindex):
        t = self.__timeseriesdata[startindex:finalindex,0]
        p = self.__timeseriesdata[startindex:finalindex,1]
        
        maxidx = p.argmax()
        minidx = p.argmin()
        
        p0 = np.array([t[maxidx],0.03,np.sqrt(p[maxidx]-p[minidx]),np.sqrt(p[minidx] if p[minidx]>0 else 0)])
        
        parameters,cov = curve_fit(gaussoffset,t,p,p0=p0,maxfev=self.__maxfev)
        fitparam = np.array([parameters[0],parameters[2]**2+parameters[3]**2,parameters[3]**2,abs(parameters[1])])
        return fitparam
    
    
    def fit_sigmoid(self,startindex,finalindex):
        t = self.__timeseriesdata[startindex:finalindex,0]
        p = self.__timeseriesdata[startindex:finalindex,1]
        
        maxidx = p.argmax()
        minidx = p.argmin()

        p0 = np.array([t[maxidx],np.sqrt(0.03),np.sqrt(0.003),np.sqrt(p[maxidx]-p[minidx]),np.sqrt(p[minidx] if p[minidx]>0 else 0)])
        
        parameters,cov = curve_fit(sigmoidoffset,t,p,p0=p0,maxfev=self.__maxfev)
        fitparam = np.array([parameters[0],parameters[3]**2+parameters[4]**2,parameters[4]**2,parameters[1]**2,parameters[2]**2])
        return fitparam
        
    
    def get_droplet_data_fit(self):
        for i in range(self.__droplets["count"]):
            if i > 0:                       
                start = int((self.__droplets["index"][i-1] + self.__droplets["index"][i])/2)
            else:
                start = 0
            if i < self.__droplets["count"] - 1:
                final = int((self.__droplets["index"][i+1] + self.__droplets["index"][i])/2)
            else:
                final = len(self.__timeseriesdata)
            
            valuesgauss = self.fit_gaussian(start,final)            
            self.__droplets.add_droplet_fitinfo("gauss",i,valuesgauss)
            
            valuesstep = self.fit_step(start,final)
            self.__droplets.add_droplet_fitinfo("step",i,valuesstep)
            
            valuessigmoid = self.fit_sigmoid(start,final)
            self.__droplets.add_droplet_fitinfo("sigmoid",i,valuessigmoid)

    def __str__(self):
        n = self.__droplets.get_count()
        if n == 0:
            return "no data"
        else:
            dtime,dmax,dbase = self.__droplets.get_data()
            
            s = ["%5d %e %e %e"%(i,dtime[i],dmax[i],dbase[i]) for i in range(n)]
            for fitdata in iter(self.__droplets):
                for i in range(n):
                    for j in range(len(fitdata[i])):
                        s[i] += " %e"%fitdata[i,j]
            
            return "\n".join(s)
    



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infile")
    parser.add_argument("-m","--minthreshold",type=float,default=.15)
    parser.add_argument("-M","--maxthreshold",type=float,default=.3)
    args = parser.parse_args()

    data = timeseries(args.infile,args.minthreshold,args.maxthreshold)
    data.get_droplet_data_minmax()
    data.get_droplet_data_fit()
    print data

if __name__ == "__main__":
    main()
