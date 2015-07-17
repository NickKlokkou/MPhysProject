'''
Created on 23 Apr 2015

@author: NicholasKlokkou
'''

import pylab

class Measure(object):
    '''
    Class responsible for storing and processing measurements from simulation
    '''
    def __init__(self, scaler, radius, N):
        self.bins = []
        self.scaler = scaler
        
    def sample(self, x):
        x = float(self.scaler)*x
        
        try:
            self.bins[int(round(x, 0))] += 1.0
        except IndexError:
            while len(self.bins) <= int(round(x, 0)):
                self.bins.append(0.0)
            self.bins[-1] += 1.0
            
    def normalise(self):
        #sum all bins and normaliseCharList
        totalSamples = sum(self.bins)
        for i in range(len(self.bins)):
            self.bins[i] = self.bins[i]/totalSamples
        
    def plot_results(self):
        #plot results
        pylab.plot(self.bins)
        pylab.show()
        
        
    
