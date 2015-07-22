'''
Created on 22 Jul 2015

@author: NicholasKlokkou
'''

'''
Created on 17 Jul 2015

@author: NicholasKlokkou
'''
import time


import pygame
import numpy as np
from math import *
import Measurements

animation = True

#VelocityDistribution = Measurements.Measure(100.0, 1.0, 1.0)

black =(0,0,0)
white = (255, 255, 255)
grey = (100, 100, 100)

def animate():
    if animation:
        pygame.display.flip()
        screen.fill(white)

def forces(ux,uy,N):
    
   # rx = np.array([[0.0 for i in range(N)] for i in range(N)])
    #ry = np.array([[0.0 for i in range(N)] for i in range(N)])
    #for i in range(N):
    #    for j in range(N):
    #        rx[i][j] = ux[i] - ux[j]
    #        ry[i][j] = uy[i] - uy[j]
    rx=ux-ux[:,None] #obscure way to get matrix rx[i,j] = ux[i]-ux[j]
    ry=uy-uy[:,None]
    d = rx*rx + ry*ry + np.eye(N)
    dinv = 1.0/d - np.eye(N)
    
    twdinv8 = 12.0*dinv**8
    twdinv14 = 12.0*dinv**14
    
    return (rx*twdinv14-rx*twdinv8).sum(axis=0),(ry*twdinv14-ry*twdinv8).sum(axis=0)+0.00

def constrain_update(s_ux, s_uy, ux,uy,pairs,N,d):
    for ij in pairs:
        i=ij[0]
        j=ij[1]
       
        si = np.array([s_ux[i],s_uy[i]])
        sj = np.array([s_ux[j],s_uy[j]])
        ri = np.array([ux[i],uy[i]])
        rj = np.array([ux[j],uy[j]])
        
        dr = rj-ri
        ds = sj-si
        constraint = (np.dot(ds,ds)-d**2) / (4.0*np.dot(ds,dr))
        
        rinext = si + constraint*(rj-ri)
        rjnext = sj - constraint*(rj-ri)
       
        ux[i] = rinext[0]
        uy[i] = rinext[1]
        ux[j] = rjnext[0]
        uy[j] = rjnext[1]
        
    return ux, uy

d = 1.0#distance dimers are apart
dt = 0.0009
N=4
L = 12.0
xmax = 0.5*L - 0.5
t=0.0
ux = np.zeros(N,float)
uy = np.zeros(N,float)
vx = np.zeros(N,float)
vy = np.zeros(N,float)

pairs = []

for k in range(int(float(N/2))):
    pairs.append([2*k,2*k+1])

print pairs

halfNsqrt=ceil(sqrt(float(N)/2))

for ij in range(N/2):
    i,j = divmod(ij, halfNsqrt)
    ux[2*ij]=1.3*i - L/2.0 + 1 -d/2
    uy[2*ij]=1.3*j - L/2.0 + 1
    ux[2*ij+1]=1.3*i - L/2.0 + 1 +d/2
    uy[2*ij+1]=1.3*j - L/2.0 + 1 
pygame.init()
Lscreen = 500.0
screen_size = 500.0
sigma = 0.02
circle_size = int(screen_size*sigma+0.5)
screen = pygame.display.set_mode([int(screen_size),int(screen_size)])
pygame.display.set_caption('Molecules with Lennard-Jones potential. N = '+str(N))
clock = pygame.time.Clock()

vx[0]=7.0 #initial conditions: one particle moves
vy[0]=5.0
ax,ay=forces(ux,uy,N)



done = False
while done == False:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            done = True
    for nsteps in range(40):
        t= t+dt
        
        
        #ux=ux+vx*dt+ax*dt*dt/2.
        #uy=uy+vy*dt+ay*dt*dt/2.
        
        s_ux=ux+vx*dt+ax*dt*dt/2.
        s_uy=uy+vy*dt+ay*dt*dt/2.
        ux, uy =constrain_update(s_ux, s_uy, ux,uy,pairs,N,d)
        vx=vx+ax*dt/2.
        vy=vy+ay*dt/2.
        
        ax,ay=forces(ux,uy,N)
        vx=vx+ax*dt/2.
        vy=vy+ay*dt/2.
        
        vx[np.abs(ux)>xmax]*=-1.0 #bounce off walls
        vy[np.abs(uy)>xmax]*=-1.0
        
        
     
    for i in range(N):
        #if t > 1.0: VelocityDistribution.sample(sqrt(vx[i]**2 + vy[i]**2))
        if animation: 
            pygame.draw.circle(screen, black, (int((Lscreen*(ux[i]+L/2)/L)),int((Lscreen*(uy[i]+L/2)/L))), int(circle_size), 0)
    animate()
    

#VelocityDistribution.plot_results()