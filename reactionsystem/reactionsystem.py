#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
from scipy import stats

class reactionsystem:
    def __init__(self,indexset = ""):
        self.__indexset = indexset.replace("0","")
        self.__numpops = len(self.__indexset)
        self.__n = {}
        for r in indexset:
            self.__n[r] = 0

        self.__recreaterv = True
        self.__reactionrates = np.array((0))
        self.__reactants = np.array(("0"))
        self.__products = np.array(("0"))
        self.__numreactions = 1

        self.__time = 0.
        self.__steps = 0
    
    def load_populations_from_file(self,filename = None,permissive = False):
        try:
            data = np.genfromtxt(filename,dtype=(str,int))
        except:
            raise IOError("Could not load populations from file '%s'"%filename)
        for i in range(len(data[:,0])):
            p = data[i,0]
            n = data[i,1]
            self.set_population(p,n,permissive)
    
    def load_reactions_from_file(self,filename = None,permissive = False):
        try:
            data = np.genfromtxt(filename)
        except:
            raise IOError("Could not load reactions from file '%s'"%filename)
        for i in range(len(data[:,0])):
            a = data[i,0] # reactants
            b = data[i,1] # products
            r = data[i,2] # rate
            self.add_reaction(a,b,r,permissive)
    
    
    def set_population(self,population = "0",value = 0,permissive = False):
        if isinstance(population,str):
            if len(population.replace("0","")) >= 1:
                for p in population.replace("0",""):
                    if p in self.__indexset:
                        self.__n[p] = value
                    else:
                        if permissive and p != "0":
                            self.__indexset += p
                            self.__n[p] = value
        
    
    def add_reaction(self,reactants,products,rate,permissive = False):
        correctreaction = True
        if rate <= 0:
            correctreaction = False
        if correctreaction:
            for r in reactants:
                if not r in self.__indexset.replace("0",""):
                    if permissive:
                        self.__indexset += r
                        self.__n[r] = 0
                    else:
                        correctreaction = False
            for r in products:
                if not r in self.__indexset.replace("0",""):
                    if permissive:
                        self.__indexset += r
                        self.__n[r] = 0
                    else:
                        correctreaction = False
        if correctreaction:
            self.__reactants     = np.append(self.__reactants,reactants)
            self.__products      = np.append(self.__products,products)
            self.__reactionrates = np.append(self.__reactionrates,rate)
            self.__numreactions += 1
            self.__recreaterv = True
    
    
    def isavailable(self,populations = "0"):
        a = True
        if populations != "0":
            for p in populations.replace("0",""):
                if self.__n[p] == 0:
                    a = False
        else:
            a = False
        return a
    
    def recreate_RV(self):
        self.__totalrate = np.sum(self.__reactionrates)
        if self.__totalrate > 0.:
            rk = self.__reactionrates/self.__totalrate
            ak = np.arange(self.__numreactions)
            self.__nextreaction = stats.rv_discrete(values = (ak,rk))
            self.__recreaterv = False
            return True
        else:
            return False
    
    def step(self):
        if self.__recreaterv:
            if not self.recreate_RV():
                raise ValueError("Total rate of all reactions is 0")
            
        nextreaction = 0
        i = 0
        while not self.isavailable(self.__reactants[nextreaction]):
            nextreaction = self.__nextreaction.rvs(size=1)[0]
            
        for r in self.__reactants[nextreaction].replace("0",""):
            self.__n[r] -= 1
        for r in self.__products[nextreaction].replace("0",""):
            self.__n[r] += 1
        
        self.__steps += 1
        self.__time  += np.random.exponential(1/self.__totalrate)
        return self.__steps


    def get_populations(self, populations = None):
        a = np.array((),dtype = int)
        if populations == None:
            listpops = self.__indexset
        else:
            listpops = populations.replace("0","")
        for r in listpops:
            a = np.append(a,self.__n[r])
        return a
    
    def get_time(self):
        return self.__time
    
    def set_time(self,time):
        self.__time = time
    
    def get_step(self):
        return self.__steps
    
    def print_reactions(self):
        if self.__numreactions > 1:
            print "# Reactants\tProducts\trate"
            print "# ============================================="
            for i in range(1,self.__numreactions):
                print "# %s\t->\t%s\t%e"%(self.__reactants[i],self.__products[i],self.__reactionrates[i])
        else:
            print "# No reactions defined"
    
    
    def is_present(self,reactant):
        a = True
        if isinstance(reactant,str):
            for r in reactant:
                if r in self.__indexset.replace("0",""):
                    if self.__n[r] == 0:
                        a = False
                else:
                    a = False
        else:
            a = False
        return a



def main():
    
    r = reactionsystem(indexset = "PNRAG")
    r.set_population("PNRA",100)
    r.add_reaction("PR","PP",1.)
    r.add_reaction("NR","NN",1.)
    r.add_reaction("P","PG",1.)
    r.add_reaction("GA","0",1.)
    r.add_reaction("PRA","R",1.)
    r.add_reaction("NRA","R",1.)
    r.print_reactions()
    
    while r.is_present("R"):
        s =  r.step()
        if s%10 == 0:
            print r.get_time(),r.get_populations()
    
if __name__ == "__main__":
    main()
    
    
    
    
    
        
