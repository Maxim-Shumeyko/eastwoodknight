c_list = []
#front coef
f = open('cal1','r')
flist = f.readlines()
for i in range(2,419,9):
    try:
        s = flist[i].split()
        c_list.append( (s[2],s[4]) )
    except: print 'front ',i
f.close()
#back coef
#f = open('zcl2.o','r')
#flist = f.readlines()
for i in range(446,1594,9):
    try:
        s = flist[i].split()
        c_list.append( (s[2],s[4]) )
    except: print 'back ', i,' ',flist[i]    
#f.close()
#side coef
#f = open('zcl1.o','r')
#flist = f.readlines()
for i in range(1641,1661,7):
    try:
        s = flist[i].split()
        c_list.append( (s[2],s[4]) )
    except: print 'side ', i,' ',flist[i]  
#
for i in range(1662,1679,8):
    try:
        s = flist[i].split()
        c_list.append( (s[2],s[4]) )
    except: print 'side ', i,' ',flist[i]    
#f.close()

f = open('coef_alpha','w')
for i in c_list:
    f.write( i[0] + ' '+i[1]+'\n')
f.close()
    
    
