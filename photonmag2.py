#Last updated: 9/20/2017
#Description: This script calculates the SED for the 2175 A dust absorbers.

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
import numpy as np
import math
import string
import os
import linecache
from scipy import interpolate
from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from decimal import Decimal



# Calculates the distance btween two given points

def distance(x, y, x1, y1):
    return math.sqrt((x - x1)**2 + (y - y1)**2)



# Finds the mean of all photon counts where the value is 3 sigma above the mean

def photoncount(scidata, radius1, radius2):
    flux = 0
    length = 0
    #print(np.shape(scidata))
                
    mean, median, std = sigma_clipped_stats(scidata, sigma=3.0, iters=5)
    
    for i in range(len(scidata)):
        for j in range(len(scidata[0])):
            if scidata[i][j] > 0 and distance(i, j, 500, 500) <= radius1 and distance(i, j, 500, 500) >= radius2:
                flux += scidata[i][j]

    return flux


def sedfitter(flux, wavelength):
    flux = 10**(flux / -2.5)
    return flux * 3631 * 10**(-23) * 29979245800/(wavelength * 10**-8)**2 * 10**-8 * 2.44**2


def errorcalc(filename):
    band = fits.open(filename)
    scidata = band[0].data.astype(float)
    outside = []
    count = 0
    total_flux = 0
    
    for i in range(len(scidata)):
        for j in range(len(scidata)):
            if distance(i, j, 50, 50) >= 50 and distance(i, j, 50, 50) <= 60:
                outside.append(scidata[i][j])
                count += 1

    mean, median, stddev = sigma_clipped_stats(outside, sigma=3.0, iters=5)
    return stddev * math.sqrt(count)



total_flux = []

hdulist = fits.open('2175TotComb13-g.fit')
scidata = hdulist[0].data.astype(float)
nmgy_count = 0

for i in range(len(scidata)):
    for j in range(len(scidata[0])):
        if scidata[i][j] > 0 and distance(i, j, 50, 50) >= 3 and distance(i, j, 50, 50) <= 16:
            nmgy_count += scidata[i][j]

gmag = -2.5 * math.log10(nmgy_count) + 22.5
print(gmag)
total_flux.append(10**(gmag/-2.5) * 3631 * 10**-23 * 29979245800/(4770 * 10**-8)**2 * 10**-8 * 2.44**2)


hdulist = fits.open('2175TotComb13-r.fit')
scidata = hdulist[0].data.astype(float)
nmgy_count = 0

for i in range(len(scidata)):
    for j in range(len(scidata[0])):
        if scidata[i][j] > 0 and distance(i, j, 50, 50) >= 3 and distance(i, j, 50, 50) <= 16:
            nmgy_count += scidata[i][j]

rmag = -2.5 * math.log10(nmgy_count) + 22.5
print(rmag)
total_flux.append(10**(rmag/-2.5) * 3631 * 10**-23 * 29979245800/(6231 * 10**-8)**2 * 10**-8 * 2.44**2)


hdulist = fits.open('2175TotComb13-i.fit')
scidata = hdulist[0].data.astype(float)
nmgy_count = 0

for i in range(len(scidata)):
    for j in range(len(scidata[0])):
        if scidata[i][j] > 0 and distance(i, j, 50, 50) >= 3 and distance(i, j, 50, 50) <= 16:
            nmgy_count += scidata[i][j]

imag = -2.5 * math.log10(nmgy_count) + 22.5
print(imag)
total_flux.append(10**(imag/-2.5) * 3631 * 10**-23 * 29979245800/(7625 * 10**-8)**2 * 10**-8 * 2.44**2)


hdulist = fits.open('2175TotComb13-z.fit')
scidata = hdulist[0].data.astype(float)
nmgy_count = 0

for i in range(len(scidata)):
    for j in range(len(scidata[0])):
        if scidata[i][j] > 0 and distance(i, j, 50, 50) >= 3 and distance(i, j, 50, 50) <= 16:
            nmgy_count += scidata[i][j]

zmag = -2.5 * math.log10(nmgy_count) + 22.5
print(zmag)
total_flux.append(10**(zmag/-2.5) * 3631 * 10**-23 * 29979245800/(9134 * 10**-8)**2 * 10**-8 * 2.44**2)

gflux = gmag
rflux = rmag
iflux = imag
zflux = zmag

gsed = total_flux[0]
rsed = total_flux[1]
ised = total_flux[2]
zsed = total_flux[3]
gerror_low = -2.5 * math.log10(10**((-2.5 * math.log10(errorcalc('2175REFSUBComb13-g.fit')) + 22.5) / -2.5) + 10**(gflux / -2.5))
rerror_low = -2.5 * math.log10(10**((-2.5 * math.log10(errorcalc('2175REFSUBComb13-r.fit')) + 22.5) / -2.5) + 10**(rflux / -2.5))
ierror_low = -2.5 * math.log10(10**((-2.5 * math.log10(errorcalc('2175REFSUBComb13-i.fit')) + 22.5) / -2.5) + 10**(iflux / -2.5))
zerror_low = -2.5 * math.log10(10**((-2.5 * math.log10(errorcalc('2175REFSUBComb13-z.fit')) + 22.5) / -2.5) + 10**(zflux / -2.5))
gerror_high = -2.5 * math.log10(-10**((-2.5 * math.log10(errorcalc('2175REFSUBComb13-g.fit')) + 22.5) / -2.5) + 10**(gflux / -2.5))
rerror_high = -2.5 * math.log10(-10**((-2.5 * math.log10(errorcalc('2175REFSUBComb13-r.fit')) + 22.5) / -2.5) + 10**(rflux / -2.5))
ierror_high = -2.5 * math.log10(-10**((-2.5 * math.log10(errorcalc('2175REFSUBComb13-i.fit')) + 22.5) / -2.5) + 10**(iflux / -2.5))
zerror_high = -2.5 * math.log10(-10**((-2.5 * math.log10(errorcalc('2175REFSUBComb13-z.fit')) + 22.5) / -2.5) + 10**(zflux / -2.5))

print("%s: %f, %f, %f, %s" % ('g-band', gflux, gerror_low, gerror_high, '%.4E' % Decimal(gsed)))
print("%s: %f, %f, %f, %s" % ('r-band', rflux, rerror_low, rerror_high, '%.4E' % Decimal(rsed)))
print("%s: %f, %f, %f, %s" % ('i-band', iflux, ierror_low, ierror_high, '%.4E' % Decimal(ised)))
print("%s: %f, %f, %f, %s" % ('z-band', zflux, zerror_low, zerror_high, '%.4E' % Decimal(zsed)))



fig, ax = plt.subplots()
error = [[sedfitter(gerror_low, 4770) - gsed, sedfitter(rerror_low, 6231) - rsed, sedfitter(ierror_low, 7625) - ised, sedfitter(zerror_low, 9134) - zsed], [gsed - sedfitter(gerror_high, 4770), rsed - sedfitter(rerror_high, 6231), ised - sedfitter(ierror_high, 7625), zsed - sedfitter(zerror_high, 9134)]]
#yerr = np.multiply(error, [[6, 6, 6, 6], [6, 6, 6, 6]])
yerr = error

print(error)

ax.errorbar([4770, 6231, 7625, 9134], [gsed, rsed, ised, zsed], yerr = yerr, fmt='bo')
plt.axis([3000, 10000, 0, 7 * 10**-17])
plt.xlabel('Wavelength [$\AA$]')
plt.ylabel('Flux Density [erg sec^-1 cm^-2 $\AA$^-1]')
ax.set_title('SED of DUST Absorbers, 0.37 <= Zabs < 0.55')

plt.show()





"""
spline = interpolate.interp2d(np.arange(len(scidata)), np.arange(len(scidata)), scidata)
scidata = spline(np.arange(0, len(scidata), 0.1), np.arange(0, len(scidata), 0.1))
#scidata *= 1.15
SBarray = []
outter = []
scale = cosmo.kpc_proper_per_arcmin(1.00) * u.arcmin / u.kiloparsec * 0.396 / 60
print(scale)

for j in range(1, 9):
    #print(5 * j / scale)
    f = photoncount(scidata, 60 * (25 / 3)**((1.0 / 8)*j) / scale, 60 * (25 / 3)**((1.0 / 8)*(j - 1)) / scale)
        
    print("%f, %f" % (60 * (25 / 3)**((1.0 / 8)*j) / scale, 60 * (25 / 3)**((1.0 / 8)*(j - 1)) / scale))
    outter.append(6 * (25 / 3)**((1.0 / 8)*j))
    f /= 100
    mag = -2.5 * np.log10(f) + 22.5
    #print(mag)

    # Calculate the surface brightness by adding the amount of kpc, similar to traditional mg/arcsec^2 formula
    surface_brightness = mag + 2.5 * np.log10(math.pi * ((60 * (25 / 3)**((1.0 / 8) * j))**2 - (60 * (25 / 3)**((1.0 / 8) * (j - 1)))**2) / 100)
    print(surface_brightness)
    SBarray.append(surface_brightness)

print(SBarray)


plt.plot(np.unique(outter), np.poly1d(np.polyfit(outter, SBarray, 3))(np.unique(outter)))
plt.plot(outter, SBarray, 'ro')

plt.axis([6, 60, 38, 28])
plt.xscale('log')
plt.xlabel('R [kpc]')
plt.ylabel('Total SB [mag/kpc^2]')
plt.show()
"""
