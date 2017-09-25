#Last updated: 9/23/2017
#Description: This code calculates the Surface Brightness (apparent magnitude / kpc^2) at anulli with geomtrically increasing radii and plots the data points 

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import numpy as np
import scipy as sp
import os
import linecache
import math
from scipy import interpolate
from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from decimal import Decimal


# Calculates the distance btween two given points

def distance(x, y, x1, y1):
    return math.sqrt((x - x1)**2 + (y - y1)**2)


# Finds the mean of all photon counts where the value is 3 sigma above the mean

def photoncount(scidata, radius1, radius2):
    flux = 0
    length = 0
    #print(np.shape(scidata))
                
    for i in range(len(scidata)):
        for j in range(len(scidata[0])):
            if distance(i, j, 500, 500) <= radius1 and distance(i, j, 500, 500) >= radius2:
                flux += scidata[i][j]
                length += 1

    return flux / length



# Find the number of points within 10 kpc

def getlength(scidata, radius1, radius2):
    length = 0
    mean, median, stddev = sigma_clipped_stats(scidata, sigma=3.0, iters=5)
    for i in range(len(scidata)):
        for j in range(len(scidata[0])):
            if distance(i, j, 500, 500) <= radius1 and distance(i, j, 500, 500) >= radius2:
                length += 1

    return length


def photoncount2(scidata, radius1, radius2):
    flux = 0
    length = 0
    #print(np.shape(scidata))
                
    for i in range(len(scidata)):
        for j in range(len(scidata[0])):
            if scidata[i][j] > 0 and distance(i, j, 500, 500) <= radius1 and distance(i, j, 500, 500) >= radius2:
                flux += scidata[i][j]
                length += 1

    return flux / length


#def getdist(z):
    
"""
dirs = os.listdir('/data/marvels/billzhu/stacked/')
writer = open('SBMagnitudes2.txt', 'w')

success = 0
fail = 0

for i in range(len(dirs)):
    filename = '/data/marvels/billzhu/stacked/' + dirs[i]
    print(filename)
    hdulist = fits.open(filename)
    scidata = hdulist[0].data.astype(float)
    SBarray = []
    scale = cosmo.kpc_proper_per_arcmin(float(dirs[i].split('_')[2])) * u.arcmin / u.kiloparsec * 0.396 / 60
    print(scale)
    
    for j in range(1, 7):
        count = photoncount(scidata, 10 * j / scale)
        if count == 0:
            continue
        f = count/(10**8 * 2000)
        mag = -2.5 * np.log10(f)
        print(mag)

        # Calculate the surface brightness by adding the amount of kpc, similar to traditional mg/arcsec^2 formula
        surface_brightness = mag + 2.5 * np.log10(math.pi/4 * (10 * j)**2)
        print(surface_brightness)
        SBarray.append(surface_brightness)

    writer.write(dirs[i] + '\t' + str(mag) + '\t' + str(SBarray) + '\n')
    success += 1



writer.write('\n\n\nSuccess: ' + str(success))
writer.write('\nFail: ' + str(fail))
writer.close()        
    # How to convert counts to mag in image where no reference stars are available
    # Need to convert reference stars to same redshift, then calculate instrumental magnitude
    # Need apparent magnitude of reference star as well
    # AstroImageJ?



    # Which flux20 to use in final image?
    # Rescaling seems ok based off of flux20 data

    # Bigger cutouts means masking array?
    # Their method is janky
"""


"""
x = np.array([0, 1, 2, 3])
y = np.array([-1, 0.2, 0.9, 2.1])
A = np.vstack([x, np.ones(len(x))]).T
m, c = np.linalg.lstsq(A, y)[0]
plt.plot(x, y, 'o', label='Original data', markersize=10)
plt.plot(x, m*x + c, 'r', label='Fitted line')
plt.legend()
plt.show()
"""



"""
hdulist = fits.open('TotComb13-z.fit')
scidata = hdulist[0].data.astype(float)
mean, median, stddev = sigma_clipped_stats(scidata, sigma=3.0, iters=5)
#scidata -= median
spline = interpolate.interp2d(np.arange(len(scidata)), np.arange(len(scidata)), scidata)
scidata = spline(np.arange(0, len(scidata), 0.1), np.arange(0, len(scidata), 0.1))
#scidata *= 1.15
SBarray = []
error_low = []
error_high = []
outter = []
scale = cosmo.kpc_proper_per_arcmin(0.48) * u.arcmin / u.kiloparsec * 0.396 / 60
print(scale)

for j in range(1, 9):
    #print(5 * j / scale)
    f = photoncount(scidata, 60 * (25 / 3)**((1.0 / 8)*j) / scale, 60 * (25 / 3)**((1.0 / 8)*(j - 1)) / scale)
        
    print("%f, %f" % (60 * (25 / 3)**((1.0 / 8)*j) / scale, 60 * (25 / 3)**((1.0 / 8)*(j - 1)) / scale))
    outter.append(6 * (25 / 3)**((1.0 / 8)*j))
    f /= 100
    f *= 2000 * 10**8
    error = math.sqrt(f / 4.8 + getlength(scidata, 60 * (25 / 3)**((1.0 / 8)*j) / scale, 60 * (25 / 3)**((1.0 / 8)*(j - 1)) / scale) / 100 * stddev) / 2000 / 10**8
    f /= 2000 * 10**8
    
    #print("%f, %f" % (f, error))
    mag = -2.5 * np.log10(f)
    el_mag = -2.5 * np.log10(f - error) + 2.5 * np.log10(math.pi * ((60 * (25 / 3)**((1.0 / 8) * j))**2 - (60 * (25 / 3)**((1.0 / 8) * (j - 1)))**2) / 100)
    eh_mag = -2.5 * np.log10(f + error) + 2.5 * np.log10(math.pi * ((60 * (25 / 3)**((1.0 / 8) * j))**2 - (60 * (25 / 3)**((1.0 / 8) * (j - 1)))**2) / 100)
    print("%f, %f" % (el_mag, eh_mag))

    # Calculate the surface brightness by adding the amount of kpc, similar to traditional mg/arcsec^2 formula
    surface_brightness = mag + 2.5 * np.log10(math.pi * ((60 * (25 / 3)**((1.0 / 8) * j))**2 - (60 * (25 / 3)**((1.0 / 8) * (j - 1)))**2) / 100)
    print(surface_brightness)
    SBarray.append(surface_brightness)
    error_low.append(el_mag)
    error_high.append(eh_mag)

print(SBarray)

# fit least-squares with an intercept
#print(np.array(outter))
#print(np.array(SBarray))
#w = np.linalg.lstsq(np.vstack([outter, np.ones(8)]).T, SBarray)[0]

# plot best-fit line
#plt.plot(outter, w[0]*np.array(outter) + w[1], 'r')

fig, ax = plt.subplots()
error = [np.array(SBarray) - np.array(error_high), np.array(error_low) - np.array(SBarray)]
#yerr = np.multiply(error, [[6, 6, 6, 6], [6, 6, 6, 6]])
yerr = error

print(error)

#ax.errorbar(outter, SBarray, yerr = yerr, fmt='bo')
(_, caps, _) = ax.errorbar(outter, SBarray, yerr=yerr, capsize=5, elinewidth=2, fmt='bo')

for cap in caps:
    cap.set_color('black')
    cap.set_markeredgewidth(1.5)

plt.plot(np.unique(outter), np.poly1d(np.polyfit(outter, SBarray, 3))(np.unique(outter)), 'b')
#plt.plot(outter, SBarray, 'bo')
"""






fit, ax = plt.subplots()

hdulist = fits.open('2175SUBComb22-g.fit')
scidata = hdulist[0].data.astype(float)
mean, median, stddev = sigma_clipped_stats(scidata, sigma=3.0, iters=5)
print(median)
#scidata -= median
spline = interpolate.interp2d(np.arange(len(scidata)), np.arange(len(scidata)), scidata)
scidata = spline(np.arange(0, len(scidata), 0.1), np.arange(0, len(scidata), 0.1))
#scidata *= 1.15

SBarray = []
outter = []
scale = cosmo.kpc_proper_per_arcmin(1.44) * u.arcmin / u.kiloparsec * 0.396 / 60
print(scale)
#print(cosmo.kpc_proper_per_arcmin(1.03))
#print(cosmo.kpc_proper_per_arcmin(5))

#print("%s" % ('%.4E' % Decimal(photoncount(scidata, 250, 200))))
#print(-2.5 * math.log10(photoncount(scidata, 250, 200)) + 2.5 * math.log10(scale**2))

#print(-2.5 * math.log10(2/(10**8 * 2000)) + 2.5 * np.log10(scale**2))
#print(scale)

for j in range(1, 12):
    #print(5 * j / scale)
    f = photoncount(scidata, 60 * (25 / 3)**((1.0 / 8)*j) / scale, 60 * (25 / 3)**((1.0 / 8)*(j - 1)) / scale)
        
    print("%f, %f" % (60 * (25 / 3)**((1.0 / 8)*j) / scale, 60 * (25 / 3)**((1.0 / 8)*(j - 1)) / scale))
    outter.append(6 * (25 / 3)**((1.0 / 8)*j))
    #f /= 100

    mag = -2.5 * np.log10(f) + 22.5
    if np.isnan(mag):
        mag = 38
    
    #print(mag)

    # Calculate the surface brightness by adding the amount of kpc, similar to traditional mg/arcsec^2 formula
    surface_brightness = mag + 2.5 * math.log10(scale**2)#np.log10(math.pi * ((60 * (25 / 3)**((1.0 / 8) * j))**2 - (60 * (25 / 3)**((1.0 / 8) * (j - 1)))**2) / 100)

    print(surface_brightness)
    SBarray.append(surface_brightness)
    

print(SBarray)

# fit least-squares with an intercept
#print(np.array(outter))
#print(np.array(SBarray))
#w = np.linalg.lstsq(np.vstack([outter, np.ones(8)]).T, SBarray)[0]

# plot best-fit line
#plt.plot(outter, w[0]*np.array(outter) + w[1], 'r')
plt.plot(np.unique(outter), np.poly1d(np.polyfit(outter, SBarray, 3))(np.unique(outter)), 'r')
plt.plot(outter, SBarray, 'ro')
#plt.plot(outter, error_low, 'go')
#plt.plot(outter, error_high, 'go')

plt.axis([6, 120, 38, 28])

majorLocator = MultipleLocator(2)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(0.5)
ax.yaxis.set_major_locator(majorLocator)
ax.yaxis.set_major_formatter(majorFormatter)
ax.yaxis.set_minor_locator(minorLocator)


#plt.axes().yaxis.set_tick_params(which='minor', right = 'off')
plt.xscale('log')
plt.xlabel('R [kpc]', fontsize = 12)
plt.ylabel('Total SB [mag/kpc^2]', fontsize = 12)
plt.title('z-band 0.37 <= Zabs < 0.55', fontsize=20)

MG_patch = mpatches.Patch(color='red', label='Mg II Absorber SB Profile')
NET_patch = mpatches.Patch(color='blue', label='Net SB Profile')
plt.legend(handles=[MG_patch, NET_patch])
plt.show()

    
