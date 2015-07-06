'''
Created on 4 Jul 2015

@author: NicholasKlokkou
'''
import math
import numpy as np
import pylab


   
class PairCorrelation(object):
    
    def __init__(self, scalar):
        
        self.gr_bins = []
        self.scalar = scalar
        self.gr_count = 0.0
        self.wallLength = 500.0
        self.radius = 30.0
        self.N = 10
        self.rho = float(self.N) / self.wallLength**2
        
        self.max_r = 5* self.radius
        self.d_r = 1.0
        
    def get_shell_size(self, pos, distance):
        theta = 0.0
        phi = 0.0
        alpha = 0.0
        shell = 2*math.pi*distance
        
        overlap = False
        if distance > (self.wallLength/2.0 - abs(pos[0])):
            theta = math.acos((self.wallLength/2.0-abs(pos[0]))/distance)
            overlap = True
                    
        if distance > (self.wallLength/2.0 - abs(pos[1])):
            phi = math.acos((self.wallLength/2.0-abs(pos[1]))/distance)
            overlap = True
                    
        if overlap:
            alpha = 2*theta*phi/math.pi
            shell = distance*(2*theta +2*phi - alpha)
       
        return shell

    
        
    def auto_correlation_function(self, hardDisks):
        
        for i in range(self.N):
            
            r = self.radius
            distances = []
            
            for j in range(self.N):
                distances.append(np.linalg.norm(hardDisks[i].position-hardDisks[j].position))
                
            while r <= self.max_r:
                n = 0.0
                shell = self.get_shell_size(hardDisks[i].position, r)
                for dist in distances:
                    if dist > r and dist <= (r+self.d_r):
                        n+=1
                g_r = n/(self.rho * shell)
                
                self.sample(g_r, r)
                
                self.gr_count +=1
                
                r += self.d_r
                
    def sample(self, gr, r):
        r = float(self.scalar)*r 
        try:
            self.gr_bins[int(round(r, 0))] += gr
        except IndexError:
            while len(self.gr_bins) <= int(round(gr, 0)):
                self.gr_bins.append(0.0)
            self.gr_bins[-1] += gr
        
    def normalise(self):
        for a in range(len(self.gr_bins)):
            self.gr_bins[a]=self.gr_bins[a]/self.gr_count
        
    def plot_results(self):
        #plot results
        pylab.plot(self.gr_bins)
        pylab.show()
