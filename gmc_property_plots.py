#import libraries
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np

#load file
mytable = Table.read('m100.co10_props_cprops.fits')

#Radius vs line width best fit line
R = np.arange(10,1000)
S = np.power(np.pi,1/2)*R/3.4
sigma = np.power(S,0.5)

#Virial Mass best fit line
M_vir = 540*np.power(R,2)

#Luminous Mass best fit line
M_lum = 2000*np.power(np.power(M_vir/39,1.234567901)/130,0.8)

#Virial mass vs. luminous mass plot
figure = plt.figure(figsize=(4.5,4)) #figure size in inches
line1, line2 = plt.loglog(M_lum,M_vir,mytable['MASS_EXTRAP'],mytable['VIRMASS_EXTRAP_DECONV'])
line1.set_linestyle('-')
line2.set_linestyle('None')
line2.set_marker('.')
plt.xlabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$') 
plt.ylabel(r'$M_{\mathrm{vir}}\ (M_{\odot})$')
plt.tight_layout() 	
plt.savefig('MlumMvir_matplotlib.png')

#Radius vs. line width(velocity dispersion) plot
figure = plt.figure(figsize=(4.5,4)) #figure size in inches
line1, line2 = plt.loglog(R,sigma, mytable['RADRMS_EXTRAP_DECONV'],mytable['VRMS_EXTRAP_DECONV'])
line1.set_linestyle('-')
line2.set_linestyle('None')
line2.set_marker('.')
plt.ylabel(r'$\sigma\ (km\ s^{-1})$') 
plt.xlabel(r'$R\ (pc)$')
plt.tight_layout() 	
plt.savefig('LwRad_matplotlib.png')

#Luminous mass vs. radius plot
figure = plt.figure(figsize=(4.5,4)) #figure size in inches
line1, line2 = plt.loglog(R,M_lum,mytable['RADRMS_EXTRAP_DECONV'],mytable['MASS_EXTRAP'])
line1.set_linestyle('-')
line2.set_linestyle('None')
line2.set_marker('.')
plt.xlabel(r'$R\ (pc)$') 
plt.ylabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$')
plt.tight_layout() 	
plt.savefig('MlumRad_matplotlib.png')