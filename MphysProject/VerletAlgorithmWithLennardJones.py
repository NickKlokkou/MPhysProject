'''
Created on 17 Jul 2015

@author: NicholasKlokkou
'''



import pygame
import numpy as np
from math import *

animation = True


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
    
    return (rx*twdinv14-rx*twdinv8).sum(axis=0),(ry*twdinv14-ry*twdinv8).sum(axis=0)+0.50
 
dt = 0.0009
N=60
L = 12
xmax = 0.5*L - 0.5
t=0.0
ux = np.zeros(N,float)
uy = np.zeros(N,float)
vx = np.zeros(N,float)
vy = np.zeros(N,float)
Nsqrt=ceil(sqrt(N))

for ij in range(N):
    i,j = divmod(ij, Nsqrt)
    ux[ij]=1.3*i - L/2.0 + 1
    uy[ij]=1.3*j - L/2.0 + 1
  
pygame.init()
Lscreen = 500
screen_size = 500
sigma = 0.02
circle_size = int(screen_size*sigma+0.5)
screen = pygame.display.set_mode([screen_size,screen_size])
pygame.display.set_caption('Molecules with Lennard-Jones potential. N = '+str(N))
clock = pygame.time.Clock()

vx[0]=0.0 #initial conditions: one particle moves
vy[0]=0.0
ax,ay=forces(ux,uy,N)



done = False
while done == False:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            done = True
    for nsteps in range(40):
        t= t+dt
        
        ux=ux+vx*dt+ax*dt*dt/2.
        uy=uy+vy*dt+ay*dt*dt/2.
        vx=vx+ax*dt/2.
        vy=vy+ay*dt/2.
        ax,ay=forces(ux,uy,N)
        vx=vx+ax*dt/2.
        vy=vy+ay*dt/2.
        
        vx[np.abs(ux)>xmax]*=-1.0 #bounce off walls
        vy[np.abs(uy)>xmax]*=-1.0
        
        
    if animation: 
        for i in range(N):
            pygame.draw.circle(screen, black, (int((Lscreen*(ux[i]+L/2)/L)),int((Lscreen*(uy[i]+L/2)/L))), int(circle_size), 0)
    animate()
