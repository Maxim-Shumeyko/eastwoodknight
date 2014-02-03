# -*- coding: cp1251 -*-
import numpy as np

class Opener:
    """
    Opener(filename) -> Открывает файл и записывает данные в виде списка 1х16 в
    соответствии с известным форматом. Данны доступны через поле .mtr или через
    вызов самого оператора ()
    Методы:
    save(filename) -> сохранить в известном формате в файл
    """
    def __init__(self,filename):
        
        if str(filename).isdigit():
            filename = 'tsn.' + filename
        elif 'tsn.' in filename:
            pass
        else:
            raise ValueError, 'wrong filename'
        self.mtr = np.fromfile(filename, dtype = 'uint16', count = -1, sep = '')
        self.mtr.shape = (-1,14)
        

    def save(self,filename):
        """
        .save(filename) - save to file
        """
        np.savetxt(filename,self.mtr,'%5d',delimiter = ' ', newline = '\n')

    def count_2events(self):
        events = np.where( (self.mtr[:,1] >0) & (self.mtr[:,2]>0) &
                               (self.mtr[:,7] > 0) & (self.mtr[:,8] > 0) )
        #np.savetxt('double_events',self.mtr[events,9],'%5d',delimiter = ' ', newline = '\n')
        print float(len(events))/self.mtr.shape[0]*100,' %'


        
    def timeof(self):
        """
        .timeof() - calculate the whole time of the file
        """
        time = self.mtr[-1,4] - self.mtr[0,4]
        if time < 0:
            time += 65536
        time += ( self.mtr[-1,3] - self.mtr[0,3] )*65536
        time_hours = (time / 1000000)/3600
        time_min = ((time / 1000000)%3600)/60
        time_sec = (time / 1000000)%60
        print 'hours = %d; minutes = %d; seconds = %d' % (time_hours,
                                                          time_min, time_sec)
    def __call__(self):
        return self.mtr

	def __getitem__(self,i):
		return self.mtr[i]

#############################################################################################

class Spectr:
    """
    Создает спектр по матрице со списком событий, получаемой с помощью Opener.
    Spectr( mtr, # матрица с исходным списком
            e_type = 'alpha', # энергетический тип спектра,
                              # например альфа-частицы {'A','a','alpha'}
                              # или спектр осколков {'F','f','fission'}
            dislocation = 'front', # определяет набор детекторов, где регистрируются частицы
                                   # фронтальные детекторы фокальной сборки {'front'}
                                   # задние детекторы фокальной сборки {'back'}
                                   # боковые детекторы {'side'}
            tof_mark = True, # наличие сигнала о срабатывании времяпролетной камеры {True,False}
            **argkw
          )
    Данные доступны через поле .mtr или через стандартный вызов экземпляра ().
    """
            
    def __init__(self, mtr , e_type = 'alpha', dislocation = 'front', tof_mark = True, **argkw):
        
        print 'step1'
        if 'energy' in argkw:
            pass

        if 'visual' in argkw:
            pass

        if 'savetxt' in argkw:
            pass
        
        print 'step2'
        # initialization a matrix shape from parameters
        alpha_types = 'a','alpha'
        fission_types = 'F','f','fission'
        disl = 'side','front','back'
        if e_type in alpha_types:
            if dislocation == 'front':
                self.scale_size = 8192
                self.strips = 48
                self.method = 'af'
            if dislocation == 'back':
                self.scale_size = 8192
                self.strips = 128
                self.method = 'ab'
            if dislocation == 'side':
                self.scale_size = 8192
                self.strips = 16
                self.method = 'as'
        elif e_type in fission_types:
            if dislocation == 'front':
                self.scale_size = 4096
                self.strips = 48
                self.method = 'ff'
            if dislocation == 'back':
                self.scale_size = 8192
                self.strips = 128
                self.method = 'fb'
            if dislocation == 'side':
                self.scale_size = 4096
                self.strips = 16
                self.method = 'fs'
        else:
            raise ValueError,'wrong parameters'
        #initialize tof_mark condition
        if tof_mark:
            tof_cond = '(mtr[:,6] > 0)'
        else:
            tof_cond = '(mtr[:,6] == 0)'    
        
        self.mtr = np.zeros((self.strips,self.scale_size))
        print 'step3:',self.strips,self.scale_size
        # filing matrix (self.mtr) for alpha-front spectr
        if self.method == 'af':
            print 'step4:af'
            #events = np.where( (mtr[:,0]%16 <= 2) & (mtr[:,1] > 0) ) & eval(tof_cond))
            print 'count events in spectr:',len(events[0])
            #print tof_cond
            for i in events[0]:
                strip = mtr[i,0]%16*16 + mtr[i,2]/4096
                chan = mtr[i,1] % 8192
                if (strip > 48) or (chan > 8192):
                    print strip, chan
                self.mtr[strip,chan] += 1 
                

        # filing matrix (self.mtr) for alpha-side spectr
        if self.method == 'as':
            events = np.where( (mtr[:,0]%16 == 3) & (mtr[:,1] > 0)& eval(tof_cond) )
            for i in events[0]:
                strip = mtr[i,2]/4096
                chan = mtr[i,1] % 8192
                self.mtr[strip,chan] += 1

        # filing matrix (self.mtr) for alpha-back spectr
        # функция преобразования стрипа-канала(номер кодировщика) в физический стрип (номер детектора)
        def strip_convert(s_ch, Id): 
            if Id%2 == 0:
                return 2*s_ch - 1 + 16*Id
            else:
                return 2*( (Id - 1)*8 + s_ch )

        """def id_det(x):
            Id = []
            count = 1
            for bits in xrange(8):
                if x & count > 0:
                    Id.append(bits)
                count *= 2
            return Id
        """

        if self.method == 'ab':
            events = np.where( (mtr[:,7] > 0 ) & (mtr[:,8] > 0) & eval(tof_cond)) & ((mtr[:,0]>>8) % 16)
            for i in events[0]:
                # определение, какие разряды заполнены
                Id = ((mtr[:,0]>>8) % 16)#mtr[i,0] % 16 
                # Выбираем старший разряд (может быть неверно!)
                strip_ch = Id[-1] * 16 + 1 + mtr[i,8]/4096
                strip = strip_convert( strip_ch, Id[-1] ) - 1 
                chan = mtr[i,8] % 4096
                self.mtr[strip,chan] += 1

        # filing matrix (self.mtr) for fission-front spectr
        if self.method == 'ff':
            events = np.where( (mtr[:,0]%16 <= 2) & (mtr[:,2] > 0) & eval(tof_cond))
            for i in events[0]:
                strip = mtr[i,0]*16 + mtr[i,2]/4096
                chan = mtr[i,2] % 4096
                self.mtr[strip,chan] += 1

        # filing matrix (self.mtr) for fission-side spectr
        if self.method == 'fs':
            events = np.where( (mtr[:,0]%16 == 3) & (mtr[:,2] > 0) & eval(tof_cond))
            for i in events[0]:
                strip = mtr[i,2]/4096
                chan = mtr[i,2] % 4096
                self.mtr[strip,chan] += 1

        # filing matrix (self.mtr) for fission-back spectr             
        if self.method == 'fb':
            events = np.where( (mtr[:,7] > 0 ) & (mtr[:,8] > 0) & eval(tof_cond))
            for i in events[0]:
                # определение, какие разряды заполнены
                Id = mtr[i,7] >> 4
                # Выбираем старший разряд (может быть неверно!)
                strip_ch = Id[-1] * 16 + 1 + mtr[i,8]/4096
                strip = strip_convert( strip_ch, Id[-1] ) - 1
                chan = mtr[i,8] % 4096
                self.mtr[strip,chan] += 1

    def __call__(self, strip = None):
        if strip == None:
            return self.mtr
        else:
            return self.mtr[strip+1]

    def __getitem__(self,i):
        return self.mtr[i]


        


        
    '''
    def transform(self,mtr):
        mtr = self.mtr
        self.events = np.array([],dtype = 'uint16)
        #line = np.zeros( (1,13) , dtype = 'uint16')
        for i in self.mtr:
            ### Indentificate a type of an event
            # front and side detectors                  
            ID = mtr[i,0]
            is_event = True # the signal is in front or side detectors
            if ID in [0,2] : # a of F type
                if mtr[i,2] / 4096  == 0 and mtr[i,1] > 0:
                    ID = 0 # a
                elif mtr[i,2] > 0:
                    ID = 1 # f
                else: is_event = False
            elif ID == 3:  # side
                if mtr[i,2] / 4096  == 0 and mtr[i,1] > 0:
                    ID = 2 # a_side
                elif mtr[i,2] > 0:
                    ID = 3 # f_side
                else: is_event = False
            # back detectors
            ID_b = mtr[i,8]
            if ID_b > 0:
                if mtr[i,9] > 0:
                    ID_b = 4 # a_back
                elif mtr[i,10] > 0:
                    ID_b = 5 # F_back
                elif mtr[i,11] > 0:
                    ID_b = 6
                elif mtr[i,12] > 0:
                    ID_b = 70
        '''
                

if __name__ == '__main__':
    #open('probe1','w').write( Spectr(Opener('tsn.118').mtr,'alpha','front')() )
    #Spectr(Opener('tsn.114')(),'alpha','front')()
    #np.savetxt('probe1', spectr ,'%5d',delimiter = ' ', newline = '\n')
    #Opener('tsn.99').timelist('probe1')
    Opener('tsn.398').save('yura5.txt')
    #Opener('tsn.251').count_2events()
    #Opener('tsn.124').timeof()
    
    
        
