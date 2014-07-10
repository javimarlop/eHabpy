import os
import sys
from ehab import *
import string
from network_mapping import *

alpha = list(string.ascii_uppercase[::-1])

resultsdir = 'E:/bastilu/data/eHabitat/eHabpy'
ecor = '60165'
    
# Map a drive letter to my directory:
# TODO - make this safer (letter P may not be available) and see whether it's possible to do it wthout my password
try:
    mapped = False
    letter_index = 0
    while not mapped:
        drive_letter = alpha[letter_index]
        #mapping_command = "rNET USE " + drive_letter + ": \\netapp2\H05_Homes\bastilu metaJRC7 /USER:bastilu"

        if (mapNetworkDrive(drive_letter, "\\\\netapp2\\H05_Homes\\bastilu", "bastilu", "metaJRC7") != -1):
            mapped = True
            print("Mapped to " + drive_letter)
            indir = drive_letter + ":"

        else:
            letter_index += 1
            if (letter_index == 20):
                break
    if (mapped):      
        ehabitat(ecor, indir, resultsdir)
        unmapNetworkDrive(drive_letter)
    else:
        print "Couldn't map a drive"
except:
    unmapNetworkDrive(drive_letter)
        
