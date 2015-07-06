import HardDisk
import pygame
import Measurements
import AutoCorrelationFunction
import numpy as np
from math import sqrt

animation = False
positionSample = True
velocitySample = True
radialDistributionSample = True


N =  10
radius = 30
wallLength = 500.0



E = N
T = 0.0
wallCollisionTimes = []
diskCollisionTimes = {}
disk = {}


pygame.init()
size = (500, 500)
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

positionMeasurements = Measurements.Measure(1.0)
velocityMeasurements = Measurements.Measure(50.0)
diskDiskDistanceMeasurements = AutoCorrelationFunction.PairCorrelation(1.0)


for i in range(N):
    
    disk[i] = HardDisk.Disk()
    disk[i].set_radius(radius)
    disk[i].set_velocity([0.0, 0.0])
    if i == 0:
        disk[i].set_velocity([sqrt(2.0*E), 0.0])
    disk[i].set_random_position(wallLength)
    
 
    disk[i].get_wall_collision_time(wallLength, T)
    wallCollisionTimes.append(round(disk[i].wallCollision,6))

        
for i in range(N):
    for j in range(N):
        if i != j:
            disk[i].get_disk_disk_collision_time(disk[j], T)
            diskCollisionTimes[(i,j)] = round(disk[i].diskCollision, 6)
    

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
    
    deltaT = nextT- T
    T = nextT
    
    for i in range(N):
        if wallCollisionTimes[i] == T:
            disk[i].wall_collision(wallLength)
            disk[i].get_wall_collision_time(wallLength, T)
            wallCollisionTimes[i] = round(T + disk[i].wallCollision, 6)
            for pair in diskCollisionTimes:
                if pair[0] == i or pair[1] == i:
                    disk[pair[0]].get_disk_disk_collision_time(disk[pair[1]], T)
                    diskCollisionTimes[pair] = round(T + disk[pair[0]].diskCollision, 6)
                    
    for pair in diskCollisionTimes:
        if diskCollisionTimes[pair] == T:
            disk[pair[0]].disk_disk_collision(disk[pair[1]])
            disk[pair[0]].get_wall_collision_time(wallLength, T)
            wallCollisionTimes[pair[0]] = round(T + disk[pair[0]].wallCollision, 6)
            disk[pair[1]].get_wall_collision_time(wallLength, T)
            wallCollisionTimes[pair[1]] = round(T + disk[pair[1]].wallCollision, 6)
            
            for nextPair in diskCollisionTimes:
                if nextPair[0] == pair[0] or nextPair[1] == pair[0] or nextPair[0] == pair[1] or nextPair[1] == pair[1]:
                    if nextPair != pair and nextPair != (pair[1],pair[0]):
                        disk[nextPair[0]].get_disk_disk_collision_time(disk[nextPair[1]], T) 
                        diskCollisionTimes[nextPair] = round(T + disk[nextPair[0]].diskCollision)
                    else:
                        disk[nextPair[0]].diskCollision = float('inf')
                        disk[nextPair[1]].diskCollision = float('inf')
                        diskCollisionTimes[nextPair] = float('inf')
    
   
    if velocitySample:
        for i in range(N):
            velocityMeasurements.sample(np.linalg.norm(disk[i].velocity))
    for t in range(int(deltaT)):     
        for i in range(N):
            
            pos = disk[i].pygame_propagate(t, wallLength)
            if (T+t)%30 == 0:
                if positionSample:
                    positionMeasurements.sample(pos[0])
                
                if radialDistributionSample:
                    diskDiskDistanceMeasurements.auto_correlation_function(disk)
            
                
            if animation: pygame.draw.circle(screen, black, pos.astype(int), int(radius), 0)
      
        animate()
pygame.quit

if positionSample: 
    positionMeasurements.normalise()
    positionMeasurements.plot_results()
if velocityMeasurements:
    velocityMeasurements.normalise()
    velocityMeasurements.plot_results()
if radialDistributionSample:
    diskDiskDistanceMeasurements.normalise()
    diskDiskDistanceMeasurements.plot_results()

 


