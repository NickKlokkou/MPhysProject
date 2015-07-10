'''
A class for Hard Disk parameters
'''
import numpy as np
from math import sqrt, ceil
 
class Disk(object):
  
    def __init__(self):
        self.position = np.array([0.0, 0.0])
        self.velocity = np.array([0.0, 0.0])
        self.mass = 1.0
        self.radius = 1.0
        self.wallCollision = float('inf')
        self.diskCollision = float('inf')  #tempory
        self.previousPosition = self.position
        self.previousVelocity = self.velocity
        
        
    def set_position(self,positionVector):
        self.position = np.array(positionVector)
        
        
    def set_velocity(self, velocityVector):
        self.velocity = np.array(velocityVector)
        
    def get_previous_velocity(self): 
        self.previousVelocity = np.array([0.0,0.0]) + self.velocity
        
    def get_previous_position(self): 
        self.previousPosition = np.array([0.0,0.0]) + self.position
        
    def set_radius(self, radiusFloat):
        self.radius = radiusFloat
        
    def set_mass(self, massFloat):
        self.mass = massFloat
        
    def set_random_position(self, boxLength):
        maxL = boxLength - self.radius
        self.position = np.array([(np.random.random()-0.5)*maxL,
                                  (np.random.random()-0.5)*maxL])
        
    def set_lattice_position(self, boxLength, i, N):
        noOfDiskInAxis = ceil(sqrt(float(N)))
        space = (boxLength - 2*self.radius*1.001)
        y = float(int(float(i+1)/noOfDiskInAxis))
        x = (float(i+1)/noOfDiskInAxis) - y
        x = x*space - space/2
        y = y/noOfDiskInAxis*space - space/2
        self.position = np.array([x, y])
        
        
    def get_wall_collision_time(self, boxLength, t):
        v = self.velocity
        x = self.position
        r = self.radius
        velDir = np.array([
                            v[0]/abs(v[0]),
                            v[1]/abs(v[1])
                            ])
        tVector = np.round ((((boxLength/2.0-r)*velDir - x)/v), 8)
        
        tMin = float('inf')
    
        for t in tVector:
            if abs(t) != 0.0 and t< tMin:
                tMin = t
            
        if tMin <=0.0: print 'negetive time', tMin, x, v
        elif tMin == float('inf'): print 'infinite time', tMin, x, v
        tMin = round(tMin, 6)
        self.wallCollision = tMin
    
    def wall_collision(self, boxLength):
        HalfWall = boxLength/2.0
        x = np.around(self.position, 6)
        v = self.velocity
        r = self.radius
        for xy in range(2):
            if abs(x[xy]) > (1.00*HalfWall - r):
                x[xy] = (HalfWall - r)*(x[xy]/abs(x[xy]))
            xy_dist = abs(abs(x[xy]) - HalfWall)
            if xy_dist < (r+0.05):
                v[xy] = -v[xy]
        
        self.velocity = v
        self.position = x
        
    def get_disk_disk_collision_time(self, disk2, t):
        del_x = self.position-disk2.position
        del_v = self.velocity-disk2.velocity
        gamma = np.dot(del_x, del_v)**2 - (np.linalg.norm(del_v)**2) * ((np.linalg.norm(del_x)**2) - 4*(self.radius**2))
        if gamma > 0 and np.dot(del_x, del_v) < 0:
            t_min = round(-((np.dot(del_x, del_v)+sqrt(gamma))/(np.linalg.norm(del_v)**2)),8)
            if t_min < 0.0: 
                t_min = float('inf')
        else:
            t_min = float('inf')
        t_min = round(t_min, 6)
        self.diskCollision = t_min
        
    def disk_disk_collision(self, disk2):
        del_x = self.position - disk2.position
        del_v = self.velocity - disk2.velocity
        e_hat = del_x/np.linalg.norm(del_x)
        self.velocity = self.velocity - e_hat*np.dot(del_v, e_hat)
    
        disk2.velocity = disk2.velocity + e_hat*np.dot(del_v, e_hat)
            
    def pygame_propagate(self, t, wallLength):
        return self.previousPosition + self.previousVelocity*float(t) + wallLength/2.0
        
    def propagate(self, tNow, tNext):
        x = self.position
        v = self.velocity
        self.position = x + (v * (tNext - tNow))
        

