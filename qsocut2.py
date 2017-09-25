#Last updated: 9/23/2017
#Description: Second version of qso cutting that uses the fpObj file and quasar catalogs for optimal centroid measurements and locating in the field images.

import numpy as np
import scipy as sp
from scipy import interpolate
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.stats import mad_std
from astropy.table import Table
import photutils
from photutils import DAOStarFinder
from photutils.centroids import centroid_2dg
import os
import linecache
import string
import math
import multiprocessing
from multiprocessing import Pool
import sys
sys.setrecursionlimit(15000)
#import sdsspy
#import sdsspy.atlas
#import sdsspy.yanny



# Calculates the distance btween two given points

def distance(x, y, x1, y1):
    return math.sqrt((x - x1)**2 + (y - y1)**2)



def inbounds(x, y):
    return x < 2048 and x >= 0 and y < 1489 and y >= 0



# Method that checks if a source > 1 sigma is outside of a certain radius, and if so, masks it by putting mean value

def checkOutter(data, mean, std):
    count = 0
    for i in range(len(data)):
        for j in range(len(data[0])):
            if data[i][j] > mean + 3 * std and distance(i, j, 50, 50) > 50:
                data[i][j] = mean
    return data



# Method that checks if a sources > 3 sigma is inside of a certain radius, and if so, rejects the image

def checkInner(data, sources, qX, qY, mean, stddev, pointer):
    count = 0
    for i in range(len(sources)):
        if distance(sources['colc'][i][pointer] - qX + 50, sources['rowc'][i][pointer] - qY + 50, 50, 50) > 8:
            visited = np.zeros((len(data), len(data)), dtype=bool)
            data = floodfill(data, int(sources['colc'][i][pointer] - qX + 50), int(sources['rowc'][i][pointer] - qY + 50), mean, mean + 2 * stddev, visited)
            #print(True)
            #print("%f, %f" % (sources['colc'][i][pointer] - qX + 50, sources['rowc'][i][pointer] - qY + 50))

            """
            for j in range(-6, 6):
                for k in range(-6, 6):
                    if distance(sources['colc'][i][pointer], sources['rowc'][i][pointer], sources['colc'][i][pointer] + j, sources['rowc'][i][pointer] + k) < 6 and int(sources['rowc'][i][pointer] - qY + 50 + k) >= 0 and int(sources['rowc'][i][pointer] - qY + 50 + k) < 101 and int(sources['colc'][i][pointer] - qX + 50 + j) >= 0 and int(sources['colc'][i][pointer] - qX + 50 + j) < 101 and data[int(sources['rowc'][i][pointer] - qY + 50 + k)][int(sources['colc'][i][pointer] - qX + 50 + j)] > mean + 3 * stddev:
                        data[int(sources['rowc'][i][pointer] - qY + 50 + k)][int(sources['colc'][i][pointer] - qX + 50 + j)] = mean
                        #print(data[int(sources['colc'][i][pointer] - qX + 50 + j)][int(sources['rowc'][i][pointer] - qY + 50 + k)])
            """

    return data





def perimeter(data, mean, threshold):
    visited = np.zeros((len(data), len(data)), dtype=bool)
    for i in range(len(data)):
        if data[len(data) - 1][i] > threshold:
            data = floodfill(data, i, len(data) - 1, mean, threshold, visited)
        if data[i][len(data) - 1] > threshold:
            data = floodfill(data, len(data) - 1, i, mean, threshold, visited)
        if data[0][i] > threshold:
            data = floodfill(data, i, 0, mean, threshold, visited)
        if data[i][0] > threshold:
            data = floodfill(data, 0, i, mean, threshold, visited)
            
    return data



            


def floodfill(data, x, y, mean, threshold, visited):
    if x >= 0 and x < len(data) and y >= 0 and y < len(data) and data[y][x] >= threshold and distance(x, y, 150, 150) > 8 and visited[y][x] == False:
        data[y][x] = mean
        visited[y][x] = True
    else:
        return data

    floodfill(data, x - 1, y, mean, threshold, visited)
    floodfill(data, x + 1, y, mean, threshold, visited)
    floodfill(data, x, y - 1, mean, threshold, visited)
    floodfill(data, x, y + 1, mean, threshold, visited)

    return data
    

# Disect the objid to gain useful information like id

def objid_extract(obj_id, full=False):
    masks={'sky_version':0x7800000000000000,
           'rerun':0x07FF000000000000,
           'run':0x0000FFFF00000000,
           'camcol':0x00000000E0000000,
           'first_field':0x0000000010000000,
           'field':0x000000000FFF0000,
           'id':0x000000000000FFFF}

    run=(obj_id & masks['run']) >> 32
    rerun=(obj_id & masks['rerun']) >> 48
    camcol=(obj_id & masks['camcol']) >> 29
    field=(obj_id & masks['field']) >> 16
    id=(obj_id & masks['id']) >> 0
    sky_version=(obj_id & masks['sky_version']) >> 59
    first_field=(obj_id & masks['first_field']) >> 28

    return {'run':run,
            'rerun':rerun,
            'camcol':camcol,
            'field':field,
            'id':id,
            'first_field':first_field,
            'sky_version':sky_version}




def begin(index):
    print(index)
    i = int(index.split('-')[0])
    mgi = int(index.split('-')[1])
    color = index.split('-')[2]

    try:
        #hdulist = fits.open('/data/marvels/billzhu/2175 Reference Dataset/' + color + '/' + str(index) + '.fit')
        hdulist = fits.open('/data/marvels/billzhu/2175 Dataset/' + color + '/' + str(index) + '.fit')
        #hdulist = fits.open('/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '.fit')
        #hdulist = fits.open('/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/' + color + '/' + str(i) + '-' + str(mgi) + '-' + color + '.fit')
        scidata = hdulist[0].data.astype(float)

        #obj_table = Table.read('/data/marvels/billzhu/2175 Reference Obj/' + str(i) + '-' + str(mgi) + '.fit', hdu=1) 
        obj_table = Table.read('/data/marvels/billzhu/2175 Obj/' + str(i) + '-' + str(mgi) + '.fit', hdu=1) 
        #obj_table = Table.read('/data/marvels/billzhu/MG II Obj/0.37 - 0.55/' + str(i) + '.fit', hdu=1)
        #obj_table = Table.read('/data/marvels/billzhu/Reference Obj/0.37 - 0.55/' + str(i) + '.fit', hdu=1)
    except:
        print(str(i) + ' Can not get table')
        return
    
    #line_data = linecache.getline('Full Data.txt', i).split()
    #obj_id = int(line_data[52])
    
    print(dr12_table[mgi]['OBJ_ID'])
    obj_id = objid_extract(int(dr12_table[mgi]['OBJ_ID']), False)['id']
    print("%f, %f" % (obj_id, len(obj_table)))
    if obj_id >= len(obj_table):
        print('obj_id over table length')
        return
    quasar = obj_table[obj_id - 1]


    #print(quasar)

    pointer = 0
    if color == 'g':
        pointer = 1
    if color == 'r':
        pointer = 2
    if color == 'i':
        pointer = 3
    if color == 'z':
        pointer = 4
    if color == 'u':
        pointer = 0

    chunk_size = 50
    
    try:
        hdulist[0].header.append(('XCOORD', int(quasar['colc'][pointer]), 'x coordinate of quasar in image'), end = True)
        hdulist[0].header.append(('YCOORD', int(quasar['rowc'][pointer]), 'y coordinate of quasar in image'), end = True)
        hdulist[0].header.append(('ID', obj_id, 'id of the quasar in the fpObj file'), end = True) 
    except:
        print(str(i) + ' Unable to get coords')
        return

    #print("%f, %f" % (quasar['colc'][pointer], quasar['rowc'][pointer]))
    if inbounds(quasar['colc'][pointer] + chunk_size + 6, quasar['rowc'][pointer] + chunk_size + 6) and inbounds(quasar['colc'][pointer] - chunk_size - 5, quasar['rowc'][pointer] - chunk_size - 5):
        
        #image = scidata[int(quasar['rowc'][pointer] - 10) : int(quasar['rowc'][pointer] + 11), int(quasar['colc'][pointer] - 10) : int(quasar['colc'][pointer] + 11)]

        #xc, yc = centroid_2dg(image, mask = None)
        #print("%f, %f, %f, %f" % (quasar['colc'][pointer], quasar['rowc'][pointer], xc, yc))
        #xc += int(quasar['colc'][pointer] - 10)
        #yc += int(quasar['rowc'][pointer] - 10)

        xc = quasar['colc'][pointer]
        yc = quasar['rowc'][pointer]
        print("%f, %f" % (quasar['colc'][pointer], quasar['rowc'][pointer]))
        image = scidata[int(yc - chunk_size - 5) : int(yc + chunk_size + 6), int(xc - chunk_size - 5) : int(xc + chunk_size + 6)]
        spline = interpolate.interp2d(np.arange(int(xc - chunk_size - 5), int(xc + chunk_size + 6)),
                                      np.arange(int(yc - chunk_size - 5), int(yc + chunk_size + 6)),
                                      image)
        xrang = np.arange(xc - chunk_size, xc + chunk_size + 1)
        yrang = np.arange(yc - chunk_size, yc + chunk_size + 1)

        if len(xrang) > 2 * chunk_size + 1:
            xrang = xrang[:-1]
        if len(yrang) > 2 * chunk_size + 1:
            yrang = yrang[:-1]

        shifted = spline(xrang, yrang)


        bkg_sigma = mad_std(scidata)
        mean1, median1, std1 = sigma_clipped_stats(shifted, sigma=3.0, iters=5)
        print("%f, %f, %f" % (mean1, median1, std1))
        #daofind = DAOStarFinder(fwhm = 2, threshold=10*bkg_sigma)
        #sources = daofind.find_stars(shifted)
        #print(sources)


        mag22 = 0
        if color == 'g':
            mag22 = 2500
        if color == 'r':
            mag22 = 1900
        if color == 'i':
            mag22 = 1400
        if color == 'z':
            mag22 = 300

        shifted = checkInner(shifted, obj_table, xc, yc, mean1, std1, pointer)

        #print("%f, %f" % (mean1, std1))
        shifted = perimeter(shifted, mean1, mean1 + 3 * std1)

        hdulist[0].header['XCOORD'] = xc
        hdulist[0].header['YCOORD'] = yc

        try:
            #shifted /= hdulist[0].header['NMGY']
            #fits.writeto('/data/marvels/billzhu/2175 Reference Quasar Cut/' + color + '/' + str(index) + '_REF.fit', shifted, hdulist[0].header, clobber = True)
            fits.writeto('/data/marvels/billzhu/2175 Quasar Cut/' + color + '/' + str(index) + '_DUST.fit', shifted, hdulist[0].header, clobber = True)
            #fits.writeto('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_MG.fit', shifted, hdulist[0].header, clobber = True)
            #fits.writeto('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/' + color + '/' + str(i) + '-' + str(mgi) + '-' + color + '_REF.fit', shifted, hdulist[0].header, clobber = True)
        except:
            print('Header is corrupt')
            return
    else:
        print('Can not cut in bounds')
        return



# Main method that can be wired for multiprocessing purposes using Pool

if __name__ == '__main__':
    #multiprocessing.set_start_method('spawn')
    dust_table = Table.read('final_catalog_full.fit')
    dr12_table = Table.read('DR12Q.fits')
    

    #dirs = os.listdir('/data/marvels/billzhu/2175 Dataset/')
    #dirs = os.listdir('/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/')

    #gdirs = os.listdir('/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/g/')
    #rdirs = os.listdir('/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/r/')
    #idirs = os.listdir('/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/i/')
    #zdirs = os.listdir('/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/z/')
    #udirs = os.listdir('/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/u/')

    #gdirs = os.listdir('/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/g/')
    #rdirs = os.listdir('/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/r/')
    #idirs = os.listdir('/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/i/')
    #zdirs = os.listdir('/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/z/')

    
    #gdirs = os.listdir('/data/marvels/billzhu/2175 Dataset/g/')
    rdirs = os.listdir('/data/marvels/billzhu/2175 Dataset/r/')
    idirs = os.listdir('/data/marvels/billzhu/2175 Dataset/i/')
    zdirs = os.listdir('/data/marvels/billzhu/2175 Dataset/z/')
    
    
    #gcheck_dirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/g/')
    #icheck_dirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/i/')
    #zcheck_dirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/z/')
    #print(len(check_dirs))
    rangelist = []
    #rangelist.append('42725-g')

    #rangelist.append(gdirs[0].split('.')[0])

    #rangelist.append('29066-g')


    
    #for d in gdirs:
    #    index = d.split('.')[0]
        #if d.split('.')[0] + '_MG.fit' not in gcheck_dirs:
    #    rangelist.append(index)
        #begin(index)

    #print(len(rangelist))
    #begin(rangelist[1])
    
        #if d.split('-')[0] + '_MG.fit' not in check_dirs:
        #rangelist.append(index + '-r')
        #rangelist.append(index + '-i')
        #rangelist.append(index + '-z')
        #rangelist.append(index + '-u')


    
    for d in rdirs:
        index = d.split('.')[0]
        #if d.split('-')[0] + '_MG.fit' not in rcheck_dirs:
        rangelist.append(index)
        #begin(index)
    
    for d in idirs:
        index = d.split('.')[0]
        #if d.split('-')[0] + '_MG.fit' not in icheck_dirs:
        rangelist.append(index)
        #begin(index)

    for d in zdirs:
        index = d.split('.')[0]
        #if d.split('-')[0] + '_MG.fit' not in zcheck_dirs:
        rangelist.append(index)
        #begin(index)
    
    
    
    #for d in udirs:
    #    index = d.split('.')[0]
    #    #if d.split('-')[0] + '_MG.fit' not in check_dirs:
    #    rangelist.append(index)
    

    print(len(rangelist))
    #try:
    #rlist = []
    #rlist.append('9807-g')
    pool = Pool(os.cpu_count())
    pool.map(begin, rangelist)
    #except:
    #print("Error: Unable to process file")

    print(True)
