'''
Created on 17 Jul 2015

@author: NicholasKlokkou
'''


scaling = False


import pygame
import numpy as np
from math import *
import Measurements
import pylab



dimer = True
d = 0.8
animation = True
auto_correlate = True
scalar = 10.0
gr_bins = []
max_r = 6
VelocityDistribution = Measurements.Measure(100.0, 1.0, 1.0)

black =(0,0,0)
white = (255, 255, 255)
grey = (100, 100, 100)

def animate():
    if animation:
        pygame.display.flip()
        screen.fill(white)

def forces(ux,uy,N, pairs_x, dimers):
    
    rx=ux-ux[:,None] #obscure way to get matrix rx[i,j] = ux[i]-ux[j]
    ry=uy-uy[:,None]

    rx = rx - L*np.round(rx/L,0)
    ry = ry - L*np.round(ry/L,0)
    
    d = np.sqrt(rx**2+ry**2)

    d = d + np.eye(N)
    dinv = 1.0/d - np.eye(N)
    
    if dimers: dinv=dinv*pairs_x
    
    twdinv8 = 12.0*dinv**8
    twdinv14 = 12.0*dinv**14
    global PE
    #PE = 0.5*np.sum((1.0/d)**12-2*(1.0/d)**6)
    PE = 0.5*np.sum((dinv)**12-2*(dinv)**6)
    return (rx*twdinv14-rx*twdinv8).sum(axis=0),(ry*twdinv14-ry*twdinv8).sum(axis=0)+0.00


def scaling_factor(vx,vy,N):
    
    TD = 65.0/119.8
    sum_v2 = np.sum(vx**2 + vy**2)
    factor = sqrt((N)*2.0*TD/sum_v2)#sqrt((N-1)*3*kB*TD/sum_v2)
    vx*=factor
    vy*=factor
    return vx, vy

def constrain_update(s_ux, s_uy, ux,uy,pairs,N,d,L):
    for ij in pairs:
        i=ij[0]
        j=ij[1]
        
        si = np.array([s_ux[i],s_uy[i]])
        sj = np.array([s_ux[j],s_uy[j]])
        ri = np.array([ux[i],uy[i]])
        rj = np.array([ux[j],uy[j]])
        
        dr = rj-ri
        
        dr = dr - L*np.round(dr/L,0)

        ds = sj-si
        
        ds = ds - L*np.round(ds/L,0)

        constraint = (np.dot(ds,ds)-d**2) / (4.0*np.dot(ds,dr))
        
        rinext = si + constraint*(dr)
        rjnext = sj - constraint*(dr)
        
        if abs(rinext[0]) > L:
            rinext[0] = rinext[0] - L*rinext[0]/abs(rinext[0])
        if abs(rinext[1]) > L:
            rinext[1] = rinext[1] - L*rinext[1]/abs(rinext[1])
        if abs(rjnext[0]) > L:
            rjnext[0] = rjnext[0] - L*rjnext[0]/abs(rjnext[0])
        if abs(rjnext[1]) > L:
            rjnext[1] = rjnext[1] - L*rjnext[1]/abs(rjnext[1])
        
        ux[i] = rinext[0]
        uy[i] = rinext[1]
        ux[j] = rjnext[0]
        uy[j] = rjnext[1]
        
    return ux, uy

PE_bins = []
KE_bins = []
Total_bins = []

dt = 0.00009#0.0009
N = 72
L = 12
xmax = 0.5*L - 0.5
t=0.0
ux = np.zeros(N,float)
uy = np.zeros(N,float)
vx = np.zeros(N,float)
vy = np.zeros(N,float)
Nsqrt=ceil(sqrt(N))
halfNsqrt = ceil(sqrt(N/2))
'''
for ij in range(N):
    i,j = divmod(ij, Nsqrt)
    ux[ij]=1.3*i - L/2.0 + 1
    uy[ij]=1.3*j - L/2.0 + 1
'''
pairs = []
if dimer:
    for ij in range(N/2):
        i,j = divmod(ij, halfNsqrt)
        ux[2*ij]=(i+0.5)*L/halfNsqrt-L/2 -d/2
        uy[2*ij]=(j+0.5)*L/halfNsqrt-L/2
        ux[2*ij+1]=(i+0.5)*L/halfNsqrt-L/2 +d/2
        uy[2*ij+1]=(j+0.5)*L/halfNsqrt-L/2
    

    for k in range(int(float(N/2))):
        pairs.append([2*k,2*k+1])
        pairs.append([2*k+1,2*k])
    
    pairs_x =  np.array([2.0]*N)-np.array([1.0]*N)[:,None]
    
    
    for pair in pairs:
        pairs_x[pair[0],pair[1]] = 0.0
        pairs_x[pair[1],pair[0]] = 0.0
        
else:
    for ij in range(N):
        i,j = divmod(ij, Nsqrt)
    
        ux[ij]=(i+0.5)*(L)/Nsqrt-L/2
        uy[ij]=(j+0.5)*(L)/Nsqrt-L/2
     
pygame.init()
Lscreen = 500
screen_size = 500
sigma = 0.02
circle_size = int(screen_size*sigma+0.5)
screen = pygame.display.set_mode([screen_size,screen_size])
pygame.display.set_caption('Molecules with Lennard-Jones potential. N = '+str(N))
clock = pygame.time.Clock()

vx[0]=14.0 #initial conditions: one particle moves
vy[0]=14.0
vx[1]=14.0
vy[1]=14.0
ax,ay=forces(ux,uy,N, pairs_x, dimer)


Temp_count= 0.0
Temperature = 0.0
done = False
while done == False:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            done = True
    for nsteps in range(40):
        t= t+dt
        
        
        if dimer:
            s_ux=ux+vx*dt+ax*dt*dt/2.
            s_uy=uy+vy*dt+ay*dt*dt/2.
            ux, uy = constrain_update(s_ux, s_uy, ux, uy, pairs, N, d, L)
        else:
            ux=ux+vx*dt+ax*dt*dt/2.
            uy=uy+vy*dt+ay*dt*dt/2.
        vx=vx+ax*dt/2.
        vy=vy+ay*dt/2.
        ax,ay=forces(ux,uy,N, pairs_x, dimer)
        vx=vx+ax*dt/2.
        vy=vy+ay*dt/2.
        if dimer and t > 20: 
            dimer = False
            scaling = False 
            print 'Constraint broken'
        if scaling and int(t*1000)%7==0 and t > 2: 
            vx, vy = scaling_factor(vx, vy, N)
            if t > 200:
                scaling = False
                print 'scaling stopped'
    
        if not scaling and auto_correlate and int(t*1000)%7==0:
            rx=ux-ux[:,None]
            ry=uy-uy[:,None]
            rx = rx - L*np.round(rx/L,0)
            ry = ry - L*np.round(ry/L,0)
            distances = np.sqrt(rx**2+ry**2)
            for r_group in distances:
                for r in r_group:
                    if r <= max_r and r > 0.0:
                        r = float(scalar)*r 
                        try:
                            gr_bins[int(round(r, 0))] += 1.0
                        except IndexError:
                            while len(gr_bins) <= int(round(r, 0)):
                                gr_bins.append(0.0)
                                gr_bins[-1] += 1.0
                
    
        #if t>1 and t<20:
            #vx=vx*0.9999
            #vy=vy*0.9999
            #print t
        #vx[np.abs(ux)>xmax]*=-1.0 #bounce off walls
        #vy[np.abs(uy)>xmax]*=-1.0
        for i in range(N):
            if np.abs(ux[i])>L/2:               
                ux[i]-=(L)*(np.abs(ux[i])/ux[i])
                
            if np.abs(uy[i])>L/2:
                uy[i]-=(L)*(np.abs(uy[i])/uy[i])
          
          
        KE = np.sum(vx**2+vy**2)*0.5
        KE_bins.append(KE)
        PE_bins.append(PE)
        Total_bins.append(KE+PE)
        if not scaling:
            Temperature += KE/N 
            Temp_count+=1
        '''
        if t > 1:
            KE = np.sum(vx**2+vy**2)*0.5
            print 'Kinetic energy: ' + str(KE)
            print 'Temperature: ' + str(KE/N)
            print 'Total: ' + str(PE + KE)
           ''' 
     
    if animation:
        for i in range(N):
            #if t > 1.0: VelocityDistribution.sample(sqrt(vx[i]**2 + vy[i]**2))
            pygame.draw.circle(screen, black, (int((Lscreen*(ux[i]+L/2)/L)),int((Lscreen*(uy[i]+L/2)/L))), int(circle_size), 0)
            
    animate()

for i in range(len(gr_bins)):
    gr_bins[i] = gr_bins[i] /((L**2/N)*2*pi*(float(i+1)/scalar))
    
pylab.plot(gr_bins)
pylab.show()
#VelocityDistribution.plot_results()
pylab.plot(KE_bins)
pylab.plot(PE_bins)
pylab.plot(Total_bins)
pylab.show()

print 119.8*Temperature/Temp_count