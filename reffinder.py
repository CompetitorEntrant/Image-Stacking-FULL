#Last updated: 8/28/2017
#Description: This script identifies reference quasars for Mg II absorbers and downloads the imaging data files. 

from shutil import copyfile
from astropy.io import fits
from astropy.table import Table
import os
import os.path
import numpy as np
import scipy as sp
import math
import string
import linecache
import multiprocessing
from multiprocessing import Pool

"""
for i in range(8373, 50001):
    filename = 'MG II Dataset/' + str(i) + '.fit'
    try:
        hdulist = fits.open(filename)
        counter = 0
        while counter < 4:
            
    except:
        print(str(i))
        try:
            copyfile('Final Data Extract/' + str(i) + '.fit', 'Reference Dataset/' + str(i) + '_REF.fit')
        except:
            print(str(i))
"""




def begin(index):
    print(index)
    color = index.split('-')[1]
    index = int(index.split('-')[0])
    line = linecache.getline('Full Data.txt', index)
    mgdata = line.split()
    redshift = float(mgdata[3])
    mggmag = float(mgdata[6])
    count = 0

    for i in range(1, 105784):
        if i == index:
            continue

        
        line = linecache.getline('Full Data.txt', i)
        refdata = line.split()
        z = float(refdata[3])
        refgmag = float(refdata[6])

        if abs(redshift - z) < 0.1 and abs(mggmag - refgmag) < 0.5 and (i not in usedref):
            breaker = False
            for j in range(len(mgtable)):
                if abs(mgtable['RA'][j] - float(refdata[1])) < 0.0001:
                    breaker = True
                    break

            if breaker == True:
                continue
            
            try:
                if 'FLUX20' not in fits.open('/data/marvels/jerryxu/dr7/catalog/' + str(i) + '-g.fit')[0].header.keys():
                    continue
                copyfile('/data/marvels/jerryxu/dr7/raw_catalog/' + str(i) + '.fit', '/data/marvels/billzhu/Reference Obj/0.37 - 0.55/' + str(i) + '.fit')

                copyfile('/data/marvels/jerryxu/dr7/catalog/' + str(i) + '-g.fit', '/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/g/' + str(i) + '-' + str(index) + '-g.fit')
                

                """
                if (str(i) + '_SUB.fit') in refcutdirs:
                    print(i)
                    os.rename('/data/marvels/billzhu/Reference PSF Subtract 2/0.37 - 0.55/' + str(i) + '_SUB.fit', '/data/marvels/billzhu/Reference PSF Subtract 2/0.37 - 0.55/' + str(i) + '_' + str(index) + '_SUB.fit')

                """
                # Rename the reference PSF final, rescaled images to include their MG II QSO index for easy access later


                count += 1
                usedref.append(i)

            except:
                continue

        if count == 5:
            return
        

        
if __name__ == '__main__':
    #multiprocessing.set_start_method('spawn')

    #try:

    dirs = os.listdir('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/g/')
    mgtable = Table.read('QSObased_Trimmed_SDSS_DR7_107.fits', hdu=1)
    rangelist = []
    usedref = []
    
    #mgdirs = os.listdir('/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/g/')
    #refdirs = os.listdir('/data/marvels/billzhu/Reference Dataset/')
    #refcutdirs = os.listdir('/data/marvels/billzhu/Reference PSF Subtract 2/0.37 - 0.55')
    
    #print(refcutdirs)
    for d in dirs:
        index = d.split('_')[0]
        rangelist.append(index)
        begin(index)
    
    


    # Not using multiprocessing due to shared data, will resolve later

    """
    pool = Pool(os.cpu_count())
    pool.map(begin, rangelist)
    #except:
    #print("Error: Unable to process file")     
    """
    
    
