#Last updated: 9/14/2017
#Description: This script calculates the average luminosity-weighted impact parameters of the absorber host galaxies in the different color bands. This code is not complete.

from astropy.io import fits
import numpy as np
from astropy.stats import sigma_clipped_stats
import math


def distance(x1, y1, x2, y2):
    return math.sqrt((x1 - x2)**2 + (y1 - y2)**2)

    
def ipcalc(filename):
    band = fits.open(filename)
    scidata = band[0].data.astype(float)
    ip_num = 0
    ip_den = 0

    for i in range(len(scidata)):
        for j in range(len(scidata)):
            if distance(i, j, 50, 50) <= 50 and distance(i, j, 50, 50) >= 5 and scidata[i][j] > 0:
                ip_num += -2.5 * math.log10(scidata[i][j]) * distance(i, j, 50, 50)**2 * 4.2
                ip_den += -2.5 * math.log10(scidata[i][j]) * distance(i, j, 50, 50) * 2

    ip = ip_num / ip_den
    return ip

-2.5 * math.log10(scidata[i][j])
g_ip = ipcalc('TotComb13-g.fit')
r_ip = ipcalc('TotComb13-r.fit')
i_ip = ipcalc('TotComb13-i.fit')
z_ip = ipcalc('TotComb13-z.fit')

print(g_ip)
print(r_ip)
print(i_ip)
print(z_ip)
