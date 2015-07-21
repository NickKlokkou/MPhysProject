'''
Created on 11 Jul 2015

@author: Nick
'''
import sys
try:
    N = int(sys.argv[1])
    radius = float(sys.argv[2])
except IndexError:
    print 'Not enough arguments. 1st arg is N, 2nd is radius'
    print 'taking defult values'
    N = 16
    radius = 0.1


wallLength = 1.0
