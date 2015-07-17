
import HardDisk
import pygame
import Measurements
import AutoCorrelationFunction
import numpy as np
import HardDiskConfig
from math import sqrt



animation = False    
allMeasurements = True

positionSample = True
velocitySample = False
radialDistributionSample = False

if allMeasurements:
    positionSample = True
    velocitySample = True
    radialDistributionSample = True

N =HardDiskConfig.N  #25
radius = HardDiskConfig.radius        #0.08
wallLength = HardDiskConfig.wallLength



E = N
T = 0.0
wallCollisionTimes = []
diskCollisionTimes = {}
disk = {}


pygame.init()
screenScalar = 500
size = (int(screenScalar*wallLength), int(screenScalar*wallLength))
screen = pygame.display.set_mode(size)
clock = pygame.time.Clock()
pygame.display.set_caption("Hard Disk Model")

def animate():
    if animation:
        pygame.display.flip()
        clock.tick(60)
        screen.fill(white)

black =(0,0,0)
white = (255, 255, 255)
grey = (100, 100, 100)

positionMeasurements = Measurements.Measure(float(screenScalar),radius, N)
velocityMeasurements = Measurements.Measure(float(screenScalar),radius, N)
diskDiskDistanceMeasurements = AutoCorrelationFunction.PairCorrelation(float(screenScalar), radius, N, wallLength)


for i in range(N):
    
    disk[i] = HardDisk.Disk()
    disk[i].set_radius(radius)
    disk[i].set_velocity([0.0, 0.0])
    if i == 0:
        disk[i].set_velocity([sqrt(2.001*E)/sqrt(2.0),sqrt(2.0*E)/sqrt(2.0)])
    #disk[i].set_lattice_position(wallLength, i, N)
    disk[i].set_random_position(wallLength)
    
    
 
    disk[i].get_wall_collision_time(wallLength, T)
    wallCollisionTimes.append(round(float(disk[i].wallCollision),8))

print 'placed'
for i in range(N):
    for j in range(N):
        if i != j:
            disk[i].get_disk_disk_collision_time(disk[j], T)
            diskCollisionTimes[(i,j)] = round(float(disk[i].diskCollision), 8)
    

done = False
while done == False:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            done = True
    nextT = min(wallCollisionTimes)
    minDiskTime = float('inf')
    for pair in diskCollisionTimes:
        if diskCollisionTimes[pair] < nextT:
            nextT = diskCollisionTimes[pair]
    
    for i in range(N):
        disk[i].get_previous_velocity()
        disk[i].get_previous_position()
        
        disk[i].propagate(T, nextT)
    
    deltaT = 100*screenScalar*(nextT- T)
    T = nextT
    
    for i in range(N):
        if wallCollisionTimes[i] == T:
            disk[i].wall_collision(wallLength)
            disk[i].get_wall_collision_time(wallLength, T)
            wallCollisionTimes[i] = round(float(T + disk[i].wallCollision), 8)
            for pair in diskCollisionTimes:
                if pair[0] == i or pair[1] == i:
                    disk[pair[0]].get_disk_disk_collision_time(disk[pair[1]], T)
                    diskCollisionTimes[pair] = round(float(T + disk[pair[0]].diskCollision), 8)
                    
    for pair in diskCollisionTimes:
        if diskCollisionTimes[pair] == T:
            disk[pair[0]].disk_disk_collision(disk[pair[1]])
            disk[pair[0]].get_wall_collision_time(wallLength, T)
            wallCollisionTimes[pair[0]] = round(float(T + disk[pair[0]].wallCollision), 8)
            disk[pair[1]].get_wall_collision_time(wallLength, T)
            wallCollisionTimes[pair[1]] = round(float(T + disk[pair[1]].wallCollision), 8)
            
            for nextPair in diskCollisionTimes:
                if nextPair[0] == pair[0] or nextPair[1] == pair[0] or nextPair[0] == pair[1] or nextPair[1] == pair[1]:
                    if nextPair != pair and nextPair != (pair[1],pair[0]):
                        disk[nextPair[0]].get_disk_disk_collision_time(disk[nextPair[1]], T) 
                        diskCollisionTimes[nextPair] = round(float(T + disk[nextPair[0]].diskCollision), 8)
                    else:
                        disk[nextPair[0]].diskCollision = float('inf')
                        disk[nextPair[1]].diskCollision = float('inf')
                        diskCollisionTimes[nextPair] = float('inf')
    
   
    if velocitySample:
        for i in range(N):
            velocityMeasurements.sample(np.linalg.norm(disk[i].velocity))
    for t in range(int(deltaT)):
        if radialDistributionSample and int(screenScalar*(T+t))%37 == 0 and T > 2:
            diskDiskDistanceMeasurements.auto_correlation_function(disk, t/(100*screenScalar))
        for i in range(N):
            
        
            pos = disk[i].pygame_propagate(t/(100*screenScalar), wallLength)     
            if positionSample and int(screenScalar*(T+t))%37 == 0  and T > 1: 
                for i in range(N):
                    positionMeasurements.sample(pos[0])
                
        
                
            if animation: pygame.draw.circle(screen, black, pos.astype(int)*screenScalar, int(radius*screenScalar), 0)
      
        animate()
pygame.quit

if positionSample: 
    positionMeasurements.normalise()
    positionMeasurements.plot_results()
if velocityMeasurements:
    velocityMeasurements.normalise()
    velocityMeasurements.plot_results()
if radialDistributionSample:
    #diskDiskDistanceMeasurements.normalise()
    diskDiskDistanceMeasurements.plot_results()

 


