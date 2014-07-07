from datetime import datetime
import numpy as np
import ehablib as ehab

eco_list0 = np.genfromtxt('pas/ecoregs.csv',dtype='int') # crear este archivo en subpas!
eco_list = np.unique(eco_list0)
print eco_list

m = len(eco_list)
#print pa_list[1]
for pm in range(0,m): # 3,m # 0 without the negative ecoregs!
 ecor = eco_list[pm]
 print ecor
 ehab.ehabitat(ecor,'')
print str(datetime.now())
print "BATCH END"
