'''
Created on 11 Apr 2015

@author: NicholasKlokkou
'''
import pygame

class PygameOutput(object):
    '''
    classdocs
    '''


    def __init__(self, screenSize):
    
        pygame.init()
        #size = (500, 500)
        self.screen = pygame.display.set_mode(screenSize)
        self.clock = pygame.time.Clock()
        pygame.display.set_caption("Hard Disk Model")
        