'''
Created on 4 Jul 2015

@author: NicholasKlokkou
'''
import math
import numpy as np
import pylab


   
class PairCorrelation(object):
    
    def __init__(self):
        
        self.gr_bins = []
        self.scalar = 4.0
        self.gr_count = 0.0
        self.wallLength = 500.0
        self.radius = 25.0
        self.N = 20
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
        shell = shell*self.d_r
        return shell

    def in_middle(self, pos):
        middle = True
        
        for xy in pos:
            if (self.wallLength/2) -abs(xy) < self.max_r:
                middle = False
        return middle

    def auto_correlation_function(self, hardDisks, t):
        for j in range(self.N):
            if j > 0:
                hardDisk0Position = hardDisks[0].position+hardDisks[0].velocity*t
                
                distance = np.linalg.norm(hardDisk0Position - (hardDisks[j].position+hardDisks[j].velocity*t))
                if distance <= self.max_r and self.in_middle(hardDisk0Position):
                    shell = self.get_shell_size(hardDisk0Position, distance)
                    g_r = 1.0/(self.rho * shell )
                    
                    self.sample(g_r, distance)
        

    def sample(self, gr, r):
        r = float(self.scalar)*r 
        try:
            self.gr_bins[int(round(r, 0))] += gr
        except IndexError:
            while len(self.gr_bins) <= int(round(r, 0)):
                self.gr_bins.append(0.0)
            self.gr_bins[-1] += gr
        
    def normalise(self):
        for a in range(len(self.gr_bins)):
            self.gr_bins[a]=self.gr_bins[a]/self.gr_count
        
    def plot_results(self):
        #plot results
        pylab.plot(self.gr_bins)
        pylab.show()

    

"""
        
    def auto_correlation_function(self, hardDisks):
        
        for i in range(self.N):
            i=randint(self.N)
            r = self.radius
            distances = []
            
            for j in range(self.N):
                if i!=j:
                    
                    distances.append(np.linalg.norm(hardDisks[i].position-hardDisks[j].position))
                
            while r <= self.max_r:
                n = 0.0
                shell = self.get_shell_size(hardDisks[i].position, r)
                for dist in distances:
                    if dist > r and dist <= (r+self.d_r):
                        n+=1.0
                g_r = n/(self.rho * shell)
                
                self.sample(g_r, r)
                
                self.gr_count +=1
                
                r += self.d_r
            break
"""