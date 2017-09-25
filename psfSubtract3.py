#Last updated: 9/23/2017
#Description: Final version of PCA PSF subtraction that uses the fpObj files given, gauranteeing optimal normalization, background + sky subtractions, only bright stars are used, etc. 

from astropy.io import fits
from astropy.table import Column
from astropy.table import Table
import scipy as sp
import scipy.optimize as opt
from scipy import interpolate
import numpy as np
np.set_printoptions(threshold=np.inf)
from numpy import *
import functools
import string
import decimal
import matplotlib.pyplot as plt
import fileinput
import math
import os
import photutils
from photutils import DAOStarFinder
from photutils import IRAFStarFinder
from photutils import subtract_psf
from photutils import centroid_2dg
from astropy.stats import mad_std
from astropy.stats import sigma_clipped_stats
import multiprocessing
from multiprocessing import Pool
import linecache
from sklearn.decomposition import NMF, PCA, IncrementalPCA
import astropy.units as u
import astropy.coordinates as coord
import sys
sys.setrecursionlimit(15000)



# Checks if the given coordinates are within the bounds of the image

def inbounds(x, y):
    if x > 0 and x < 2048 and y > 0 and y < 1489:
        return True
    return False



# Calculates the distance btween two given points

def distance(x, y, x1, y1):
    return math.sqrt((x - x1)**2 + (y - y1)**2)



# Method that checks if a particular residue is noise or not by drawing a ring of similar distance as the "bright" point around the centroid

def checkNoise(x, y, qX, qY, data):
    halo = []
    for i in range(42):
        for j in range(42):
            if abs(distance(i, j, qX, qY) - distance(x, y, qX, qY)) <= 2 and distance(i, j, x, y) >= 2:
                halo.append(data[j, i])
    mean, median, std = sigma_clipped_stats(halo, sigma=3.0, iters=5)

    if data[y, x] > mean + 3 * std:
        return True
    return False
                

    
# Method that calculates the total photon count within 3 sigma of the quasar centroid

def photonCount(xc, yc, sigma, data):
    count = 0
    for i in range(len(data)):
        for j in range(len(data)):
            if distance(i, j, xc, yc) <= sigma:
                count += data[i][j]
    return count



# Method that checks if a source > 1 sigma is outside of a certain radius, and if so, masks it by putting mean value

def checkOutter(data, mean, std, visited):
    count = 0
    for i in range(len(data)):
        for j in range(len(data[0])):
            if data[i][j] > mean + 3 * std and visited[i][j] == False:
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





def floodfill(data, x, y, mean, threshold, visited):
    if x >= 0 and x < len(data) and y >= 0 and y < len(data) and data[y][x] >= threshold and visited[y][x] == False:
        data[y][x] = mean
        visited[y][x] = True
    else:
        return data

    floodfill(data, x - 1, y, mean, threshold, visited)
    floodfill(data, x + 1, y, mean, threshold, visited)
    floodfill(data, x, y - 1, mean, threshold, visited)
    floodfill(data, x, y + 1, mean, threshold, visited)

    return data



def perimeter(data, mean, threshold):
    visited = np.zeros((len(data), len(data)), dtype=bool)
    for i in range(101):
        if data[len(data) - 1][i] > threshold:
            data = floodfill(data, i, len(data) - 1, mean, threshold, visited)
        if data[i][len(data) - 1] > threshold:
            data = floodfill(data, len(data) - 1, i, mean, threshold, visited)
        if data[0][i] > threshold:
            data = floodfill(data, i, 0, mean, threshold, visited)
        if data[i][0] > threshold:
            data = floodfill(data, 0, i, mean, threshold, visited)
            
    return data




def normalize(data, x, y, p1, p2, threshold, visited):
    if x >= 0 and x < len(data) and y >= 0 and y < len(data) and data[y][x] >= threshold and visited[y][x] == False:
        data[y][x] /= p1
        data[y][x] *= p2
        visited[y][x] = True
    else:
        return data

    normalize(data, x - 1, y, p1, p2, threshold, visited)
    normalize(data, x + 1, y, p1, p2, threshold, visited)
    normalize(data, x, y - 1, p1, p2, threshold, visited)
    normalize(data, x, y + 1, p1, p2, threshold, visited)

    return data



def connected(data, x, y, threshold, visited):
    if x >= 0 and x < len(data) and y >= 0 and y < len(data) and data[y][x] >= threshold and visited[y][x] == False:
        visited[y][x] = True
    else:
        return visited

    connected(data, x - 1, y, threshold, visited)
    connected(data, x + 1, y, threshold, visited)
    connected(data, x, y - 1, threshold, visited)
    connected(data, x, y + 1, threshold, visited)

    return visited



    
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




# Class that assists in choosing the optimal PSF sources by storing index and Z value

class tuplet:
    def __init__(self, i, z):
        self.i = i
        self.z = z
        
    def getSelf(z):
        return self

    def getIndex():
        return i

    def getZ():
        return z


#reader = open('Pixel Coordinates 50000.txt', 'r')
#reader = open('MG II Pixel Coordinates.txt', 'r')



# Main method for executing PSF Subtraction
# Only searches in a 600x600 box around the quasar for a PSF source
# If none are found then the file is passed on and not used

def begin(index):
    #i = int(index)
    i = int(index.split('-')[0])
    mgi = int(index.split('-')[1])
    color = index.split('-')[2]
    
    #try:
    print(index)
    #filename = 'Test Data Extract/' + str(i) + '.fit'
    #filename = str(i) + '-g.fit'


    filename = '/data/marvels/billzhu/2175 Dataset/' + color + '/' + str(index) + '.fit'
    #filename = '/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/' + color + '/' + str(index) + '.fit'
    #filename = '/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '.fit'
    hdulist = fits.open(filename)

    #qlist = fits.open('MG II Test Cut/' + str(i) + '_MG.fit')

    qlist = fits.open('/data/marvels/billzhu/2175 Quasar Cut/' + color + '/' + str(index) + '_DUST.fit')
    #qlist = fits.open('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/' + color + '/' + str(index) + '_REF.fit')
    #qlist = fits.open('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_MG.fit')

    qx = qlist[0].header['XCOORD']
    qy = qlist[0].header['YCOORD']
    obj_id = qlist[0].header['ID']
    #print("%f, %f" % (x, y))
    qlist.close()

    #except:
    #    print("No coordinates")
    #    return

    # Save some frickin time

        
    
    scidata = hdulist[0].data.astype(float)
    mean, median, std = sigma_clipped_stats(scidata, sigma=3.0, iters=5)
    #print(median)

    """
    if 'SKY' in hdulist[0].header.keys():
        scidata -= float(hdulist[0].header['SOFTBIAS'])
        scidata -= float(hdulist[0].header['SKY'])
    else:
        scidata -= median
        print(str(i) + ' No sky')
        #return
    """
    

    #print(sigma_clipped_stats(scidata, sigma=3.0, iters=5))

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


    bkg_sigma = mad_std(scidata)

    try:
        #print('/data/marvels/billzhu/MG II Obj/0.37 - 0.55/' + str(i) + '.fit')
        obj_table = Table.read('/data/marvels/billzhu/2175 Obj/' + str(i) + '-' + str(mgi) + '.fit', hdu=1)
        #obj_table = Table.read('/data/marvels/billzhu/MG II Obj/0.37 - 0.55/' + str(i) + '.fit', hdu=1)
        #obj_table = Table.read('/data/marvels/billzhu/Reference Obj/0.37 - 0.55/' + str(i) + '.fit', hdu=1)
    except:
        print(str(i) + ' No Table')
        return
    
    #line_data = linecache.getline('Full Data.txt', i).split()
    #line_data = linecache.getline('DR12 QSO.txt', i).split()
    #print(len(line_data))
    #obj_id = int(line_data[52])
    quasar = obj_table[obj_id - 1]

    gauss = 0
    #scidata /= hdulist[0].header['NMGY']
    print("%f, %f" % (qx, qy))
    data = scidata[int(qy) - 10 : int(qy) + 11, int(qx) - 10 : int(qx) + 11]
    #print(data)
    if(np.ma.count(data) >= 7):
        gauss = photutils.fit_2dgaussian(data, mask = None)

    fwhm_x = 0
    fwhm_y = 0

    #print(gauss.x_stddev)
    #print(gauss.y_stddev)
    if gauss != 0:
        fwhm_x = 2*np.sqrt(2*np.log(2))*np.sqrt(abs(gauss.x_stddev.value))
        fwhm_y = 2*np.sqrt(2*np.log(2))*np.sqrt(abs(gauss.y_stddev.value))
        #qsigma = np.sqrt(gauss.x_stddev**2 + gauss.y_stddev**2)

    print(fwhm_x)
    print(fwhm_y)
    if fwhm_x == math.nan:
        fwhm_x = 2.8
    if fwhm_y == math.nan:
        fwhm_y = 2.8
    
    # If no quasar is found, the field image is deemed corrupt and not used
      
    if quasar == 0:
        print(str(i) + ' No quasar')
        return


    # Calculate the 18 magnitude threshold

    mag18 = 0
    header = hdulist[0].header

    try:
        mag18 = header['FLUX20'] * 10 **(8. - 18/2.5)
    except:
        if color == 'g':
            mag18 = 1500
        if color == 'r':
            mag18 = 10500
        if color == 'i':
            mag18 = 8800
        if color == 'z':
            mag18 = 1900
        print(str(i) + ' MAG20 APPROX = 14000')


    
    qsocut = fits.open('/data/marvels/billzhu/2175 Quasar Cut/' + color + '/' + str(index) + '_DUST.fit')
    #qsocut = fits.open('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/' + color + '/' + str(index) + '_REF.fit')
    #qsocut = fits.open('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_MG.fit')
    qsodata = qsocut[0].data.astype(float)

    
    """
    if 'SKY' in hdulist[0].header.keys():
        qsodata -= float(hdulist[0].header['SOFTBIAS'])
        qsodata -= float(hdulist[0].header['SKY'])
    else:
        qsodata -= median
    """



    largearr = []
    stars = []
    chunk_size = 50
    diff_fwhm = 1000000


    psf_fwhm = 100
    qsovisited = connected(qsodata, 50, 50, mean + 3 * std, np.zeros((101, 101), dtype=bool))
    qmax = np.max(qsodata)

    
    for j in range(len(obj_table)):
        sx = obj_table['colc'][j][pointer]
        sy = obj_table['rowc'][j][pointer]
        if obj_table['objc_type'][j] == 6 and distance(sx, sy, qx, qy) > 5 and inbounds(sx + chunk_size + 6, sy + chunk_size + 6) and inbounds(sx - chunk_size - 5, sy - chunk_size - 5) and obj_table['psfCounts'][j][pointer] > mag18:
            #try:
            """
            preshift = scidata[int(sy - 10) : int(sy + 11), int(sx - 10) : int(sx + 11)]
            xc, yc = centroid_2dg(preshift, mask = None)
            xc += quasar['colc'][pointer] - 10
            yc += quasar['rowc'][pointer] - 10
            """

            xc = obj_table['colc'][j][pointer]
            yc = obj_table['rowc'][j][pointer]
            preshift = scidata[int(yc - chunk_size - 5) : int(yc + chunk_size + 6), int(xc - chunk_size - 5) : int(xc + chunk_size + 6)]
            spline = interpolate.interp2d(np.arange(int(xc - chunk_size - 5), int(xc + chunk_size + 6)),
                                          np.arange(int(yc - chunk_size - 5), int(yc + chunk_size + 6)),
                                          preshift)

            xrang = np.arange(xc - chunk_size, xc + chunk_size + 1)
            yrang = np.arange(yc - chunk_size, yc + chunk_size + 1)

            if len(xrang) > 2 * chunk_size + 1:
                xrang = xrang[:-1]
            if len(yrang) > 2 * chunk_size + 1:
                yrang = yrang[:-1]

            shifted1 = spline(xrang, yrang)

            bkg_sigma = mad_std(scidata)
            #mean1, median1, std1 = sigma_clipped_stats(shifted1, sigma=3.0, iters=5)
            #print("%f, %f, %f" % (mean1, median1, std1))
            #daofind = DAOStarFinder(fwhm = 2., threshold = 3.*bkg_sigma)
            #sources = daofind.find_stars(shifted1)


            """
            mag22 = 0
            if color == 'g':
                mag22 = 2500
            if color == 'r':
                mag22 = 1900
            if color == 'i':
                mag22 = 1400
            if color == 'z':
                mag22 = 300
            """

            
            shifted1 = checkInner(shifted1, obj_table, xc, yc, mean, std, pointer)
            shifted1 = perimeter(shifted1, mean, mean + 3 * std)
            #shifted1 = checkOutter(shifted1, mean1, std1)
            #print('REACHED')

            gauss1 = photutils.fit_2dgaussian(shifted1[40 : 61, 40 : 61], mask = None)
            fit_fwhm_x = 2*np.sqrt(2*np.log(2))*np.sqrt(abs(gauss1.x_stddev.value))
            fit_fwhm_y = 2*np.sqrt(2*np.log(2))*np.sqrt(abs(gauss1.y_stddev.value))
            #print("%f, %f, %f, %f" % (fit_fwhm_x, fit_fwhm_y, fwhm_x, fwhm_y))


            if abs(fit_fwhm_x - fwhm_x) < 0.25 and abs(fit_fwhm_y - fwhm_y) < 0.25:

                #shifted1 = normalize(shifted1, 50, 50, np.max(shifted1), np.max(qsodata), mean1 + 5 * std1, np.zeros((101, 101), dtype=bool))

                #visited = connected(shifted1, 50, 50, mean1 + 3 * std, np.zeros((101, 101), dtype=bool))
                smax = np.max(shifted1)
                for r in range(len(qsovisited)):
                    for c in range(len(qsovisited)):
                        if qsovisited[r][c] == True:
                            shifted1[r][c] /= obj_table['psfCounts'][j][pointer]
                            shifted1[r][c] *= obj_table['psfCounts'][obj_id - 1][pointer]
                
                #shifted1 /= np.max(shifted1)
                #shifted1 *= np.max(qsodata)
                
                
                largearr.append(np.reshape(shifted1, 10201))

                #stars.append(shifted1)
            #except:
            #    continue

    

    
    largearr = np.array(largearr)
    print(np.shape(largearr))


    # Set number of components in PCA, use incremental PCA (IPCA) due to high efficiency and speed
    
    numcomp = len(largearr)


    # Need a healthy number of sources to make the PSF fitting in order to decrease noise, setting at 5% threshold

    if len(largearr) < 10:
        print('No Sources')
        return

    print(numcomp)
    mean_vector = []
 
    #print(np.shape(largearr))

    try:
        for j in range(0, 10201):
            mean_vector.append(np.mean(largearr[:, j]))
    except:
        print("NO SOURCE FOUND")
        return

    largearr -= mean_vector
        
    ipca = IncrementalPCA(n_components=numcomp)
    ipca.fit(largearr)
    ipca_comp = ipca.components_
    #print(np.shape(ipca_comp))



    # Only use the components of the central portion of the quasar, since there may be overfitting due to strength of ipca
    
    new_comp = []
    for j in range(len(largearr)):
        temp = np.reshape(ipca_comp[j, :], (101, 101))
        new_comp.append(np.reshape(temp[47 : 54, 47 : 54], 49))

    new_comp = np.array(new_comp)
    new_comp = new_comp.T
    print(np.shape(new_comp))

    ipca_comp = ipca_comp.T

    
    #print(ipca_comp)

    #print(np.shape(largearr[0, :]))
    #print(np.shape(ipca_comp))
    total_res = 0
    max_median = 10000000




    take_final = 10


    # Final fitting of the first n components, as determined by take_final, into the quasar to build a PSF fit
    print(np.shape(ipca_comp))
    qsodata = np.reshape(qsodata, 10201)
    qsodata -= mean_vector
    #qsodata = np.reshape(qsodata, (101, 101))
    #coeff = np.dot(np.reshape(qsodata[47 : 54, 47 : 54], 49), new_comp)
    
    coeff = np.dot(qsodata, ipca_comp)

    final_fit = np.dot(ipca_comp[:, 0:take_final], coeff[0:take_final])
    final_fit += mean_vector
    final_fit = np.reshape(final_fit, (101, 101))
    #final_fit /= len(largearr)

    qsodata = np.reshape(qsodata, 10201)
    qsodata += mean_vector
    qsodata = np.reshape(qsodata, (101, 101))

    #final_fit /= final_fit[50, 50]
    #final_fit *= qsodata[50, 50]

    """
    fmax = np.max(final_fit)
    qmax = np.max(qsodata)
    for r in range(len(qsovisited)):
        for c in range(len(qsovisited)):
            if qsovisited[r][c] == True:
                final_fit[r][c] /= fmax
                final_fit[r][c] *= qmax

    print(np.min(final_fit))
    print(np.min(qsodata))
    """

    
    """
    fxc, fyc = centroid_2dg(final_fit[40 : 61, 40 : 61], mask = None)
    fxc += 40
    fyc += 40
    print("%f, %f" % (fxc, fyc))
    spline = interpolate.interp2d(np.arange(101), np.arange(101), final_fit)
    fxcr = np.arange(fxc - 50, fxc + 51)
    fycr = np.arange(fyc - 50, fyc + 51)

    if len(fxcr) > 101:
        fxcr = fxcr[1:]
    if len(fycr) > 101:
        fycr = fycr[1:]
        
    final_fit = spline(fxcr, fycr)
    """

           
    
    # Section to normalize the PSF fitting by photon count, but this is unneccesary since CORRECT PCA fits better

    
    gauss_fit = photutils.fit_2dgaussian(final_fit[40 : 61, 40 : 61], mask = None)
    fit_fwhm = 2*np.sqrt(2*np.log(2))*np.sqrt(gauss_fit.x_stddev)
    #print(fit_fwhm)

    
    print("%f, %f" % (fwhm_x, fit_fwhm))
    print("%f, %f" % (np.max(qsodata), np.max(final_fit)))
    ffwhm = max(fwhm_x, fit_fwhm)
    ffphoton_1sig = photonCount(50, 50, 4 * fit_fwhm, final_fit) 
    qsophoton_1sig = photonCount(50, 50, 4 * fit_fwhm, qsodata)

    print("%f, %f" % (qsophoton_1sig, ffphoton_1sig))
    
    
    for j in range(len(final_fit)):
        for k in range(len(final_fit)):
            if distance(50, 50, j, k) < 4 * fit_fwhm:
                final_fit[j][k] /= ffphoton_1sig
                final_fit[j][k] *= qsophoton_1sig
    
    

    print("%f, %f" % (np.max(qsodata), np.max(final_fit)))

    #final_fit /= ffphoton_1sig
    #final_fit *= qsophoton_1sig

    
    """
    line_data = linecache.getline('Full Data.txt', i).split()

    if color == 'g':
        mag = float(line_data[6])
    if color == 'r':
        mag = float(line_data[8])
    if color == 'i':
        mag = float(line_data[10])
    if color == 'z':
        mag = float(line_data[12])
    if color == 'u':
        mag = float(line_data[4])
        

    #try:
    multiplier = 10**(8 - mag / 2.5) * header['FLUX20']
    multiplier1 = quasar['psfCounts'][pointer]
    print(multiplier)
    print(str(multiplier1 - header['SKY']))
    final_fit /= ffphoton_1sig
    final_fit *= multiplier
    #except:
        #final_fit *= qsodata[50, 50]
        #return
    
    
    print("%f, %f" % (np.max(qsodata), np.max(final_fit)))
    
    """


    # Final residue from subtraction of PSF from QSO


    residue = qsodata - final_fit
    #mean, median, stddev = sigma_clipped_stats(residue[0 : 10, 0 : 10], sigma=3.0, iters=5)
    #residue -= median
    

    try:
        fits.writeto('/data/marvels/billzhu/2175 PSF Cut/' + color + '/'+ str(index) + '_PSF.fit', final_fit, hdulist[0].header, clobber = True)
        fits.writeto('/data/marvels/billzhu/2175 PSF Subtract/' + color + '/' + str(index) + '_SUB.fit', residue, hdulist[0].header, clobber = True)

        #fits.writeto('/data/marvels/billzhu/Reference PSF Cut/0.37 - 0.55/' + color + '/' + str(index) + '_PSF.fit', final_fit, hdulist[0].header, clobber = True)
        #fits.writeto('/data/marvels/billzhu/Reference PSF Subtract/0.37 - 0.55/' + color + '/' + str(index) + '_SUB.fit', residue, hdulist[0].header, clobber = True)
        #fits.writeto('/data/marvels/billzhu/MG II PSF Cut/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_PSF.fit', final_fit, hdulist[0].header, clobber = True)
        #fits.writeto('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_SUB.fit', residue, hdulist[0].header, clobber = True)

        """
        fits.writeto('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_SUB-1.fit', residue1, hdulist[0].header, clobber = True)
        fits.writeto('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_SUB-2.fit', residue2, hdulist[0].header, clobber = True)
        fits.writeto('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_SUB-3.fit', residue3, hdulist[0].header, clobber = True)
        fits.writeto('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_SUB-4.fit', residue4, hdulist[0].header, clobber = True)
        fits.writeto('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_SUB-5.fit', stars[0], hdulist[0].header, clobber = True)
        fits.writeto('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_SUB-6.fit', stars[1], hdulist[0].header, clobber = True)
        """
        
        #fits.writeto('Reference Subtract/' + str(i) + '_SUB.fit', residue, hdulist[0].header, clobber = True)
        #fits.writeto('Reference PSF Cut/' + str(i) + '_PSF.fit', final_fit, hdulist[0].header, clobber = True)
        print('\n')

        print("DONE TO BOTTOM")
    except:
        print('HEADER IS CORRUPT')


    
    
# Code that opens up a maximum of 8 processes for concurrent execution
    
if __name__ == '__main__':
    #multiprocessing.set_start_method('spawn')

    #dirs = os.listdir('/data/marvels/billzhu/2175 Quasar Cut/')
    #dirs = os.listdir('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/')

    #gdirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/g/')
    #rdirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/r/')
    #idirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/i/')
    #zdirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/z/')
    #udirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/u/')


    gdirs = os.listdir('/data/marvels/billzhu/2175 Quasar Cut/g/')
    #rdirs = os.listdir('/data/marvels/billzhu/2175 Reference Quasar Cut/r/')
    #idirs = os.listdir('/data/marvels/billzhu/2175 Reference Quasar Cut/i/')
    #zdirs = os.listdir('/data/marvels/billzhu/2175 Reference Quasar Cut/z/')
    
    #check_dirs = os.listdir('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/z/')
    rangelist = []
    #rangelist.append('99795-r')
    #rangelist.append('89-230830-g')

    #rangelist.append('100728-g')



    #rangelist.append(gdirs[0].split('_')[0])
    
    for d in gdirs:
        index = d.split('_')[0]
        #if str(index) + '_SUB.fit' not in check_dirs:
        rangelist.append(index)


    
    
    print(len(rangelist))
    #try:
    pool = Pool(os.cpu_count())
    pool.map(begin, rangelist)
    #except:
    #print("Error: Unable to process file")





