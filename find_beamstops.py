# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 12:19:27 2014

@author: eastwoodknight
"""

import read_files

def find(filenames):
    frame = read_files.read_files(filenames,strlength=14)
    beam_events = frame[ (frame['beam_marker']>0)]
    
    beam_events.to_csv('beam_events.csv')
    import numpy as np
    f = lambda x: np.hstack([[0], np.array( x.index[1:]) - np.array( x.index[:-1])])
    
    beam_events = beam_events[ f(beam_events)>1 ]
    print 'amout_of_stops = ',len(beam_events)
    beam_events.to_csv('beam_stops.csv')
    return beam_events
    
if __name__ == '__main__':
    print find('tsn.397').tail()