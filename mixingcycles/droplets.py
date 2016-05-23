#!/usr/bin/env python


import numpy as np
import argparse
import sys,math
from scipy import stats

class population:
    def __init__(self,droplets = 100,initialpool = (1e5,1e5), growthrates = (1.,2.), yields = (2.,1.), substrate = 1e5, dilution = 1e-4,mixingtime = 1e3,debug = True, dilutiontype = 0,precision = 1e-10):
        self.__droplets = droplets
        self.__numstrains = len(initialpool)
        self.__pool = np.array(initialpool)
        self.__growthrates = np.array(growthrates,dtype=float)
        self.__yields = np.array(yields,dtype=float)
        assert len(initialpool) == len(growthrates) == len(yields)
        self.__substrate = substrate
        self.__remaining_substrate_in_pool = 0
        self.__dilution = dilution
        self.__dilutiontype = dilutiontype
        self.__mixingtime = mixingtime
        self.__populations = np.outer(np.ones(self.__droplets),self.__pool*self.__dilution)
        self.__precision = precision
        self.__debug = debug
        

    def get_time_to_substrate_depletion(self,coefficients,rates,substrate,alpha=1.,maxsteps = 10000):
        z0 = 0
        if np.sum(coefficients) > 0:
            p = coefficients.argmax()
            z1 = np.log(substrate/coefficients[p]+1.)/rates[p]
            i = 0
            while ((z1-z0)/z1)**2 > self.__precision**2:
                z0 = z1
                z1 += alpha*(substrate-np.sum(coefficients*(np.exp(rates*z1)-1.)))/(np.sum(coefficients*rates*np.exp(rates*z1)))
                i+=1
                if i > maxsteps:
                    raise ValueError
            return min(z1,self.__mixingtime)
        else:
            return 0
        
        
    def dilute_fraction(self):
        self.__populations = np.random.poisson(lam = self.__pool*self.__dilution/(1.*self.__droplets),size = (self.__droplets,self.__numstrains))

    def dilute_single(self):
        self.__populations = np.zeros((self.__droplets,self.__numstrains))
        rv = stats.rv_discrete(values=(np.arange(self.__numstrains), self.__pool/np.sum(self.__pool)))
        a = rv.rvs(size=self.__droplets)
        for i in range(self.__droplets):
            self.__populations[i,a[i]] = 1
        
    def grow(self):
        if self.__dilutiontype == 0:
            self.dilute_fraction()
        elif self.__dilutiontype == 1:
            self.dilute_single()
        
        #self.__seedsize = np.sum(self.__populations,axis=1)
        #self.__seedhisto,self.__seedbins = np.histogram(self.__seedsize)
        
        self.__pool = np.zeros(self.__numstrains)
        
        self.__averagetime = 0
        self.__remaining_substrate_in_pool = 0
        n = 0
        
        for i in range(self.__droplets):
            t = self.get_time_to_substrate_depletion(self.__populations[i]/self.__yields,self.__growthrates,self.__substrate)
            self.__pool += self.__populations[i]*np.exp(self.__growthrates*t)
            self.__remaining_substrate_in_pool += self.__substrate - np.sum(self.__populations[i]/self.__yields * (np.exp(self.__growthrates*t) - 1))
            
            if t>0:
                self.__averagetime += t
                n += 1
            if self.__debug: print "    %3d %.4e %6d %6d %.4e %.4e"%(i,t,self.__populations[i][0],self.__populations[i][1],self.__populations[i][0]*np.exp(self.__growthrates[0]*t),self.__populations[i][1]*np.exp(self.__growthrates[1]*t))

        if n > 0:
            self.__averagetime /= n

        
        
    def get_populations(self):
        return np.sum(self.__populations,axis=0)/(1.*self.__droplets)
    
    
    def get_ratios(self):
        totalpool = np.sum(self.__pool)
        if totalpool > 0:
            return self.__pool/totalpool
        else:
            return np.zeros(self.__numstrains)
        
    def get_average_time(self):
        return self.__averagetime

    def get_seedhisto(self):
        b = 0.5*(self.__seedbins[1:] + self.__seedbins[:-1])
        return b,self.__seedhisto
    
    def get_remaining_substrate(self):
        return self.__remaining_substrate_in_pool





def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-n","--droplets",type=int,default=100)
    parser.add_argument("-d","--dilution",type=float,default=1e-4)
    parser.add_argument("-S","--substrate",type=float,default=1e5)
    parser.add_argument("-g","--growthrates",type=float,nargs="*",default=[1,2])
    parser.add_argument("-y","--yields",type=float,nargs="*",default=[2,1])
    parser.add_argument("-N","--initialpopulation",type=float,nargs="*",default=[1e5,1e5])
    parser.add_argument("-T","--mixingtime",type=float,default=1e5)
    parser.add_argument("-M","--maxsteps",type=int,default=10000)
    parser.add_argument("-O","--outputsteps",type=int,default=100)
    parser.add_argument("-v","--verbose",action="store_true",default=False)
    parser.add_argument("-q","--debug",action="store_true",default=False)

    args = parser.parse_args()
    
    s = population(args.droplets,args.initialpopulation,args.growthrates,args.yields,args.substrate,args.dilution/args.droplets,args.mixingtime,debug = args.debug,dilutiontype = 0)
    for i in range(args.maxsteps):
        s.grow()
        if i%args.outputsteps == 0:
            print "%6d %14.6e %14.6e %14.6e %6.2f %6.2f"%(i,s.get_ratios()[0],s.get_average_time(),s.get_remaining_substrate(),s.get_populations()[0],s.get_populations()[1])
            if args.verbose:
                b,h = s.get_seedhisto()
                for j in range(len(h)):
                    print >> sys.stderr,i,b[j],h[j]
                print >> sys.stderr,""
        
    
if __name__ == "__main__":
    main()
    
