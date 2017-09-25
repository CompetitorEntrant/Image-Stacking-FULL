from astropy.io import fits
import numpy as np
import math

"""
f1 = fits.open('RefCombine1.fit')
dat1 = f1[0].data
f2 = fits.open('RefCombine2.fit')
dat2 = f2[0].data
f3 = fits.open('RefCombine3.fit')
dat3 = f3[0].data
f4 = fits.open('RefCombine4.fit')
dat4 = f4[0].data


fits.writeto('S1.fit', dat1 - dat2, clobber = True)
fits.writeto('S2.fit', dat1 - dat3, clobber = True)
fits.writeto('S3.fit', dat1 - dat4, clobber = True)
fits.writeto('S4.fit', dat2 - dat3, clobber = True)
fits.writeto('S5.fit', dat2 - dat4, clobber = True)
fits.writeto('S6.fit', dat3 - dat4, clobber = True)
"""

mg = fits.open('2175SUBComb13-i.fit')
mgdata = mg[0].data
ref = fits.open('2175REFSUBComb13-i.fit')
refdata = ref[0].data
comb = mgdata

"""
for i in range(len(mgdata)):
    for j in range(len(mgdata[0])):
        if math.sqrt((i - 50)**2 + (j - 50)**2) <= 25:
            comb[j][i] -= refdata[j][i]
"""            
#comb = comb[8 : 37, 8 : 37]
fits.writeto('2175TotComb13-i.fit', mgdata - refdata, clobber = True)
