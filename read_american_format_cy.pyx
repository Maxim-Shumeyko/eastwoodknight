# -*- coding: utf-8 -*-
"""
Created on Thu Aug 28 13:05:29 2014

@author: eastwood
"""
import struct
import numpy as np
from pandas import DataFrame as DF, Series as Series
from time import time as show_time
import datetime
import os



def open_am_file(size,fileobj,chunk_size=10000):
    """
    Read events from binary file.
    Each event includes
        	id: 1-front, 3-side, 4-back
        	strip: (1,48)-front, (1,6)-side 
        	channel - amplitudes
        	time 
        	b_strip 1 - strip number for back detectors (1,128)-back
        	beam_marker
        	tof
        	veto
    """
    
    start = show_time()
    #fileobj = open(filename,'rb')
    cdef int block_size = 32
    cdef long int num_blocks = 0
    cdef long int offset = 0
    cdef int i, size_msv  
    
    def read_block(offset,file_object):
        data = file_object.read(32)
        values = struct.unpack('qqqhhi',data[0:32])
        time = values[0]
        time_sync = values[1]
        time_inner = values[2]
        beam_on = values[3]
        count_tof_camers = values[4]
        count_events = values[5]
        return time,time_sync,time_inner,beam_on,count_tof_camers,count_events
        
    def read_events(offset,count_events,file_object):
        size_events =count_events*8
        data = file_object.read(size_events)
        frm = str(2*count_events)+'i'
        if offset == size:
            raise ValueError('\n Value error: The end of the file')
        if size_events > 0:
            values = struct.unpack(frm,data[:size_events])
            amplitudes = values[1::2]
            positions = values[::2]
        else:
            amplitudes,positions = [],[]
        return amplitudes,positions
    
    
    list_id = []
    list_strips = []
    list_channels = []
    list_time = []
    list_tof = []
    list_back_strips = []
    list_beam = []
    list_veto = []
    
    #adding data to DataFrame 
    cdef int pos, chan     
    
    while offset < size:
        time,time_sync,time_inner,beam_on,count_tof_camers,count_events = read_block(offset,fileobj)
        offset += 32 
        try:
            amplitudes,positions = read_events(offset,count_events,fileobj)
        except ValueError:
            amplitudes,positions = [],[]
            print 'There is no data about events at the end of the file'
        offset +=count_events*8
        num_blocks += 1
        #print offset,size
        #glue lists of events
        tof = False
        veto = False
        #print offset, size
        size_msv = len(amplitudes)
        for i in range(size_msv):
            pos = positions[i]
            chan = amplitudes[i]
            if (pos == 121) or (pos==122):
                tof = True
            elif pos == 112:
                veto = True
            elif pos != 0:
                list_channels.append(chan)
                list_time.append(datetime.datetime.fromtimestamp(time))
                list_tof.append(tof)
                list_beam.append(bool(beam_on))
                list_veto.append(veto)
                if pos >= 103 and pos <= 108: #side_detectors
                    list_id.append(3)
                    list_strips.append(pos-102)
                    list_back_strips.append(0)
                elif pos> 0 and pos < 49: #front detectors
                    list_id.append(1)
                    list_strips.append(pos)
                    list_back_strips.append(0)
                elif pos < 0: #back_detectors
                    list_id.append(4)
                    list_strips.append(0)
                    list_back_strips.append(-pos)
            
        if ( len(list_id)>0 ) & (num_blocks % chunk_size == 0):
            #print num_blocks
		   
            frame = DF( {'id': list_id,'strip':list_strips,'channel':list_channels,
                      'time': list_time,
                      'b_strip 1':list_back_strips,
                      'beam_marker':list_beam,
                      'tof':list_tof,
                      'veto':list_veto
                      })
            list_id = []
            list_strips = []
            list_channels = []
            list_time = []
            list_tof = []
            list_back_strips = []
            list_beam = []  
            list_veto = []
            yield frame
                     
    #add the last part of data                
    if len(list_id)>0 :
        frame = DF( {'id': list_id,'strip':list_strips,'channel':list_channels,
                  'time': list_time,
                  'b_strip 1':list_back_strips,
                  'beam_marker':list_beam,
                  'tof':list_tof,
                  'veto':list_veto
                  })
        
        print 'convert to DF completed, time:',show_time() - start
        yield frame
    
    
    fileobj.close()   

""" 
    for name in list_files[:31]:
       Frames = read_american_format_cy.OpenAmFile(name)
       hist1,sum1 = Frames.apply_func(read_files.get_front_spectrs) 
       hist = hist.add(hist1,fill_value=0)
       sum += sum1        
            if num_blocks % chunk_size == 0:
                to_datetime = lambda x: datetime.datetime.fromtimestamp(x/1000000)
                list_time = map(to_datetime,list_time) #Series(list_time).apply(to_datetime)
               
                if len(list_id)>0:
                    frame = DF( {'id': list_id,'strip':list_strips,'channel':list_channels,
                              'time': list_time,
                              'b_strip 1':list_back_strips,
                              'beam_marker':list_beam,
                              'tof':list_tof,
                              'veto':list_veto
                              })
                    list_id = []
                    list_strips = []
                    list_channels = []
                    list_time = []
                    list_tof = []
                    list_back_strips = []
                    list_beam = []  
                    list_veto = []
                    yield frame
                else:
                    raise ValueError('No data')
"""
   
      
class OpenAmFile:
    """ 
    Class is for getting dataframes from big files.   
    Return a generator which contains chunks of 10^7 (self.chunk_size) blocks  (each contains a group of events approx. 5-50).
    Also may apply function and collect out put from each chunk (apply_func)
    """
    def __init__(self,filename,chunk_size=100000):
        f = open(filename,'rb')
        if f:
            self.size = os.stat(filename).st_size
            self.fileobj = f
            print 'File: %s; size: %dB is ready to load' %(filename,self.size)
            print 'It contains roughly about %d chunks' %(self.size/(chunk_size*56)+1)  
            self.chunk_size = chunk_size
        else:
            print 'Error while openning. Please, check up the data file.'
               
    def get_frames(self):
        return open_am_file(self.size,self.fileobj,chunk_size=self.chunk_size)
        
    def apply_func(self,func):
        frames = self.get_frames()
        output = func( frames.next() )
        for i in frames:
            output += func(i)
        return output

