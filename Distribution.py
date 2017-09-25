#Last updated: 9/24/2017
#Description: This script was used to create histograms of the redshift data for the Mg II and 2175 A dust absorbers.

from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

table_list = fits.open('Trimmed_SDSS_DR7_107.fits')
table_readMG = Table.read('Trimmed_SDSS_DR7_107.fits')
table_readDUST = Table.read('final_catalog_full.fit')

redshift_MG = []
redshift_DUST = []
average = []
for i in range(len(table_readMG)):
    redshift_MG.append(table_readMG['ZABS'][i])
    if table_readMG['ZABS'][i] >= 0.37 and table_readMG['ZABS'][i] < 0.55:
        average.append(table_readMG['ZABS'][i])

average = np.mean(average)
print(average)

average = []
for i in range(len(table_readDUST)):
    redshift_DUST.append(table_readDUST['ZABS'][i])
    average.append(table_readDUST['ZABS'][i])

average = np.mean(average)
print(average)
    
bins = np.arange(0.2, 2.5, 0.1)
arr_MG = plt.hist(redshift_MG, bins = bins, alpha = 0.1, color = 'b')
#arr_DUST = plt.hist(redshift_DUST, bins = bins, alpha = 0.1, color = 'g')


for i in range(len(bins) - 1):
    
    if arr_MG[0][i] > 0:
        plt.text(arr_MG[1][i], arr_MG[0][i], str(int(arr_MG[0][i])))

    #if arr_DUST[0][i] > 0:
    #    plt.text(arr_DUST[1][i], arr_DUST[0][i], str(int(arr_DUST[0][i])))


MG_patch = mpatches.Patch(color='blue', label='Mg II Absorber')
#DUST_patch = mpatches.Patch(color='green', label='2175 $\AA$ Dust Absorber')
plt.legend(handles=[MG_patch])
plt.xticks(np.arange(0.2, 2.5, 0.1))
plt.yticks(np.arange(0, 3600, 400))
plt.xlabel('Redshift', fontsize=12)
plt.ylabel('Number of Absorber QSOs', fontsize=12)
plt.title('Redshift Distribution of Mg II Absorbers', fontsize=20)
plt.show()
