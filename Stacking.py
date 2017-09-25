#Last updated: 9/23/2017
#Description: This script mean combines the PSF subtracted frames, applies galactic extinction corrections based on the individual QSO parameters given in the QSO catalogs, and intensity calibrates all files to a uniform scale. Blocked out sections can also further binn the images by redshift / Rest Equivalent width of the absorption lines.

from astropy.io import fits
from astropy.table import Column
from astropy.table import Table
import scipy as sp
import scipy.optimize as opt
import numpy as np
from numpy import *
import string
import decimal
import matplotlib.pyplot as plt
import fileinput
import math
import linecache

import photutils
from photutils import DAOStarFinder
from photutils import subtract_psf
from astropy.stats import mad_std
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
import os
import os.path


 #dirs = os.listdir('Reference Subtract/')


# Binn the data by Signal of 2796 A line

"""
table_list = Table.read('QSObased_Trimmed_SDSS_DR7_107.fits')
dirs = os.listdir('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/')
refcutdirs = os.listdir('/data/marvels/billzhu/Reference PSF Subtract/0.37 - 0.55/')


for j in range(1, 7):
    bin1 = []
    bin2 = []
    bin3 = []
    refbin1 = []
    refbin2 = []
    refbin3 = []
    lower_bound = round(0.37 + 0.03 * (float(j) - 1.), 2)
    upper_bound = round(0.37 + 0.03 * float(j), 2)

    print(len(table_list))
    for i in range(len(table_list)):
        if table_list['REW_MGII_2796'][i][0] >= 0.8 and table_list['REW_MGII_2796'][i][0] < 1.12 and table_list['ZABS'][i][0] >= lower_bound and table_list['ZABS'][i][0] < upper_bound and (str(table_list['INDEX_QSO'][i] + 1) + '_SUB.fit') in dirs:

            if 'FLUX20' not in fits.open('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + (str(table_list['INDEX_QSO'][i] + 1) + '_SUB.fit'))[0].header.keys():
                continue
                #print(str(table_list['INDEX_QSO'][i] + 1))
            
            bin1.append(i)
            
            for k in range(len(refcutdirs)):
                if str(table_list['INDEX_QSO'][i] + 1) in refcutdirs[k] and 'FLUX20' in fits.open('/data/marvels/billzhu/Reference PSF Subtract/0.37 - 0.55/' + refcutdirs[k])[0].header.keys():
                    refbin1.append(refcutdirs[k].split('_')[0])


        if table_list['REW_MGII_2796'][i][0] >= 1.12 and table_list['REW_MGII_2796'][i][0] < 1.58 and table_list['ZABS'][i][0] >= lower_bound and table_list['ZABS'][i][0] < upper_bound and (str(table_list['INDEX_QSO'][i] + 1) + '_SUB.fit') in dirs:

            if 'FLUX20' not in fits.open('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + (str(table_list['INDEX_QSO'][i] + 1) + '_SUB.fit'))[0].header.keys():
                continue         

            bin2.append(i)

            for k in range(len(refcutdirs)):
                if str(table_list['INDEX_QSO'][i] + 1) in refcutdirs[k] and 'FLUX20' in fits.open('/data/marvels/billzhu/Reference PSF Subtract/0.37 - 0.55/' + refcutdirs[k])[0].header.keys():
                    refbin2.append(refcutdirs[k].split('_')[0])


        if table_list['REW_MGII_2796'][i][0] >= 1.58 and table_list['ZABS'][i][0] >= lower_bound and table_list['ZABS'][i][0] < upper_bound and (str(table_list['INDEX_QSO'][i] + 1) + '_SUB.fit') in dirs:

            if 'FLUX20' not in fits.open('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + (str(table_list['INDEX_QSO'][i] + 1) + '_SUB.fit'))[0].header.keys():
                continue        

            bin3.append(i)

            for k in range(len(refcutdirs)):
                if str(table_list['INDEX_QSO'][i] + 1) in refcutdirs[k] and 'FLUX20' in fits.open('/data/marvels/billzhu/Reference PSF Subtract/0.37 - 0.55/' + refcutdirs[k])[0].header.keys():
                    refbin3.append(refcutdirs[k].split('_')[0])

                    

    scitot080_112 = np.zeros((42, 42, len(bin1)))
    scitot112_158 = np.zeros((42, 42, len(bin2)))
    scitot158 = np.zeros((42, 42, len(bin3)))


    counter = 0
    for i in bin1:
        print(str(table_list['INDEX_QSO'][i] + 1))
        filename = '/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + str(table_list['INDEX_QSO'][i] + 1) + '_SUB.fit'
        hdulist = fits.open(filename)
        scidata = hdulist[0].data.astype(float) / hdulist[0].header['FLUX20'] * 2000.
        
        scitot080_112[:, :, counter] = scidata
        counter += 1

    counter = 0
    for i in bin2:
        filename = '/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + str(table_list['INDEX_QSO'][i] + 1) + '_SUB.fit'
        hdulist = fits.open(filename)
        scidata = hdulist[0].data.astype(float) / hdulist[0].header['FLUX20'] * 2000.
        scitot112_158[:, :, counter] = scidata
        counter += 1
            
    counter = 0
    for i in bin3:
        filename = '/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + str(table_list['INDEX_QSO'][i] + 1) + '_SUB.fit'
        hdulist = fits.open(filename)
        scidata = hdulist[0].data.astype(float) / hdulist[0].header['FLUX20'] * 2000.
        scitot158[:, :, counter] = scidata
        counter += 1

    
    reftot080_112 = np.zeros((42, 42, len(refbin1)))
    reftot112_158 = np.zeros((42, 42, len(refbin2)))
    reftot158 = np.zeros((42, 42, len(refbin3)))

    for i in range(len(refbin1)):
        filename = '/data/marvels/billzhu/Reference PSF Subtract/0.37 - 0.55/' + refbin1[i] + '_SUB.fit'
        hdulist = fits.open(filename)
        scidata = hdulist[0].data.astype(float) / hdulist[0].header['FLUX20'] * 2000.
        reftot080_112[:, :, i] = scidata[0 : 42, 0 : 42]

    for i in range(len(refbin2)):
        filename = '/data/marvels/billzhu/Reference PSF Subtract/0.37 - 0.55/' + refbin2[i] + '_SUB.fit'
        hdulist = fits.open(filename)
        scidata = hdulist[0].data.astype(float) / hdulist[0].header['FLUX20'] * 2000.
        reftot112_158[:, :, i] = scidata[0 : 42, 0 : 42]

    for i in range(len(refbin3)):
        filename = '/data/marvels/billzhu/Reference PSF Subtract/0.37 - 0.55/' + refbin3[i] + '_SUB.fit'
        hdulist = fits.open(filename)
        scidata = hdulist[0].data.astype(float) / hdulist[0].header['FLUX20'] * 2000.
        reftot158[:, :, i] = scidata[0 : 42, 0 : 42]
     

    scitot080_112_final = np.zeros((42, 42))
    scitot112_158_final = np.zeros((42, 42))
    scitot158_final = np.zeros((42, 42))
    reftot080_112_final = np.zeros((42, 42))
    reftot112_158_final = np.zeros((42, 42))
    reftot158_final = np.zeros((42, 42))
    
    for i in range(42):
        for k in range(42):
            scitot080_112_final[i][k] = np.mean(scitot080_112[i][k][:])
            scitot112_158_final[i][k] = np.mean(scitot112_158[i][k][:])
            scitot158_final[i][k] = np.mean(scitot158[i][k][:])
            reftot080_112_final[i][k] = np.mean(reftot080_112[i][k][:])
            reftot112_158_final[i][k] = np.mean(reftot112_158[i][k][:])
            reftot158_final[i][k] = np.mean(reftot158[i][k][:])

            
    fits.writeto('080_112_' + str(lower_bound) + '_' + str(upper_bound) + '_MGcomb.fit', scitot080_112_final, clobber = True)
    fits.writeto('112_158_' + str(lower_bound) + '_' + str(upper_bound) + '_MGcomb.fit', scitot112_158_final, clobber = True)
    fits.writeto('158_' + str(lower_bound) + '_' + str(upper_bound) + '_MGcomb.fit', scitot158_final, clobber = True)

    
    fits.writeto('080_112_' + str(lower_bound) + '_' + str(upper_bound) + '_REFcomb.fit', reftot080_112_final, clobber = True)
    fits.writeto('112_158_' + str(lower_bound) + '_' + str(upper_bound) + '_REFcomb.fit', reftot112_158_final, clobber = True)
    fits.writeto('158_' + str(lower_bound) + '_' + str(upper_bound) + '_REFcomb.fit', reftot158_final, clobber = True)

    list1 = list()
    header1 = fits.Header(list1)
    header1.append(('MGIIQSO', len(bin1), 'Number of MG II QSO stacked, not number of absorbers'))
    header1.append(('REFQSO', len(refbin1), 'Number of REF QSO stacked'))
    fits.writeto('/data/marvels/billzhu/stacked/080_112_' + str(lower_bound) + '_' + str(upper_bound) + '_SUBcomb.fit', scitot080_112_final - reftot080_112_final, header1, clobber = True)

    list1 = list()
    header1 = fits.Header(list1)
    header1.append(('MGIIQSO', len(bin2), 'Number of MG II QSO stacked, not number of absorbers'))
    header1.append(('REFQSO', len(refbin2), 'Number of REF QSO stacked'))
    fits.writeto('/data/marvels/billzhu/stacked/112_158_' + str(lower_bound) + '_' + str(upper_bound) + '_SUBcomb.fit', scitot112_158_final - reftot112_158_final, header1, clobber = True)

    list1 = list()
    header1 = fits.Header(list1)
    header1.append(('MGIIQSO', len(bin3), 'Number of MG II QSO stacked, not number of absorbers'))
    header1.append(('REFQS)', len(refbin3), 'Number of REF QSO stacked'))
    fits.writeto('/data/marvels/billzhu/stacked/158_' + str(lower_bound) + '_' + str(upper_bound) + '_SUBcomb.fit', scitot158_final - reftot158_final, header1, clobber = True)
    


        



"""    

dirs = os.listdir('/data/marvels/billzhu/2175 PSF Subtract/g/')
count = len(dirs)
print(str(count))
#scitot1 = np.zeros((42, 42, len(dirs)))

scitot = np.zeros((101, 101))
#scitot2 = np.zeros((42, 42, len(dirs)//2))
j = 0

#tableDR12 = Table.read('DR12Q.fits')

counter = 0
mean = []
num = 0
for i in range(len(dirs)):
    #print(str(i))

    try:
        filename = '/data/marvels/billzhu/2175 PSF Subtract/g/' + dirs[i]
        hdulist = fits.open(filename)
        scidata = hdulist[0].data.astype(float)

        """
        if np.mean(scidata) < -0.75:
            continue
        
        if 'FLUX20' not in hdulist[0].header.keys():
            scidata /= 1900
            scidata *= 2000
        else:
            scidata /= hdulist[0].header['FLUX20']
            #print(hdulist[0].header['FLUX20'])
            mean.append(hdulist[0].header['FLUX20'])
            num += 1
            scidata *= 2000
        
        #if scidata[50, 41] > 80:
        #    print(dirs[i])

        #scidata *= hdulist[0].header['NMGY']
        
        
        linedata = linecache.getline('Full Data.txt', int(dirs[i].split('-')[0])).split()

        
        color = 'r'
        multiplier = 0
        if color == 'g':
            multiplier = 10 ** (0.4 * 0.736 * float(linedata[14]))
        if color == 'r':
            multiplier = 10 ** (0.4 * 0.534 * float(linedata[14]))
        if color == 'i':
            multiplier = 10 ** (0.4 * 0.405 * float(linedata[14]))
        if color == 'z':
            multiplier = 10 ** (0.4 * 0.287 * float(linedata[14]))

        scidata *= multiplier
        """
        
        if np.max(scidata) < 1:
            scitot += scidata
            counter += 1
        #else:
        #    print("%s, %f" % (dirs[i], np.max(scidata)))

        #if scidata[57, 48] > 60:
        #    print(dirs[i])
        
        """
        if scidata[50, 50] > -20:
            scitot[:, :, j] = scidata
            j += 1
        else:
            scitot[:, :, j] = np.zeros((101, 101))
        """
        
    except:
        print("Error: File not found")
        j += 1

print(np.median(mean))
print(counter)

"""
j = 0   
for i in range(len(dirs)//2, len(dirs)):
    print(str(i))
    try:
        filename = '/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + dirs[i]
        hdulist = fits.open(filename)
        scidata = hdulist[0].data.astype(float)
        scitot2[:, :, j] = scidata
        j += 1
    
    except:
        scitot2[:, :, j] = np.zeros((42, 42))
        print("Error: File not found")
        j += 1
"""

#print(scitot)
i = 0
j = 0

medcomb = np.zeros((101, 101))
#medcomb2 = np.zeros((42, 42))


for i in range(len(medcomb)):
    for j in range(len(medcomb)):
        #mean, median, stddev = sigma_clipped_stats(scitot[i][j][:], sigma_lower=-3.0, iters=5)
        #print("%f, %f, %f" % (mean, median, stddev))
        medcomb[i][j] += scitot[i][j] / counter

"""
for i in range(42):
    for j in range(42):
        #mean, median, stddev = sigma_clipped_stats(scitot[i][j][:], sigma_lower=-3.0, iters=5)
        #print("%f, %f, %f" % (mean, median, stddev))
        medcomb2[i][j] += np.mean(scitot2[i][j][:])
"""

#fits.writeto('RefHalf1.fit', medcomb, clobber = True)
#fits.writeto('RefHalf2.fit', medcomb2, clobber = True)
fits.writeto('2175SUBComb22-g.fit', medcomb, clobber = True)

"""
mgres = fits.open('MGCombine.fit')
mgdata = mgres[0].data.astype(float)

res = mgdata - medcomb

fits.writeto('Result.fit', res, clobber = True)
"""
