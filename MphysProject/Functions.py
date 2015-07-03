import numpy as np
from math import sqrt

def get_disk_disk_collision_time(disk1, disk2):
    
    
    del_x = disk1.position-disk2.position
    del_v = disk1.velocity-disk2.velocity
    
    gamma = np.dot(del_x, del_v)**2 - (np.linalg.norm(del_v)**2) * ((np.linalg.norm(del_x)**2) - 4*(disk1.radius**2))
    
    if gamma > 0 and np.dot(del_x, del_v) < 0:
        t_min = round(-((np.dot(del_x, del_v)+sqrt(gamma))/(np.linalg.norm(del_v)**2)),8)
        if t_min < 0.0: 
            t_min = float('inf')
    else:
        t_min = float('inf')
   
    return t_min