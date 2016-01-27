#import libraries
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from galaxies import Galaxy
import astropy.units as u

#load file
mytable = Table.read('m100.co10_props_cprops.fits')


##VARIABLE DEFINITIONS##

#Radius vs line width best fit line
R = np.arange(10,1000)
S = np.power(np.pi,1/2)*R/3.4
sigmav = np.power(S,0.5)
one = np.arange(10000,10000000000,10000)

#Virial Mass best fit line
M_vir = 540*np.power(R,2)

#Luminous Mass best fit line
M_lum = 2000*np.power(np.power(M_vir/39,1.234567901)/130,0.8)

#Sigma_0, Mass Density, R_gal
sigma0 = mytable['VRMS_EXTRAP_DECONV']/np.sqrt(mytable['RADRMS_EXTRAP_DECONV'])  
M_den = mytable['VIRMASS_EXTRAP_DECONV']/(np.pi*np.power(mytable['RADRMS_EXTRAP_DECONV'],2)) 

mygalaxy = Galaxy("M100")
cpropstable = Table.read('m100.co10_props_cprops.fits')
rgal=mygalaxy.radius(ra = cpropstable['XPOS'], dec = cpropstable['YPOS'])

#Make array of indexes for low R_gal group
A = []
n = 0
test = rgal.value
for item in test:
	if item < 1200:
		A = A + [n]
	n = n +1
	
	
##PLOTS##

#Virial mass vs.luminous mass plot
figure = plt.figure(figsize=(4.5,4)) #figure size in inches
line1, line2 = plt.loglog(one,one, mytable['MASS_EXTRAP'],mytable['VIRMASS_EXTRAP_DECONV'])
line1.set_linestyle('-')
line1.set_color('k')
line2.set_linestyle('None')
line2.set_marker('.')
for item in A:
	plt.loglog(mytable['MASS_EXTRAP'][item],mytable['VIRMASS_EXTRAP_DECONV'][item],marker='.',c='b')
plt.loglog(mytable['MASS_EXTRAP'][92],mytable['VIRMASS_EXTRAP_DECONV'][92],marker='.',c='r')
plt.loglog(mytable['MASS_EXTRAP'][132],mytable['VIRMASS_EXTRAP_DECONV'][132],marker='.',c='m')
plt.xlabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$') 
plt.ylabel(r'$M_{\mathrm{vir}}\ (M_{\odot})$')
plt.tight_layout() 	
plt.savefig('MlumMvir_matplotlib.png')

#Radius vs. line width(velocity dispersion) plot
figure = plt.figure(figsize=(4.5,4)) #figure size in inches
line1, line2 = plt.loglog(R,sigmav, mytable['RADRMS_EXTRAP_DECONV'],mytable['VRMS_EXTRAP_DECONV'])
line1.set_linestyle('-')
line1.set_color('k')
line2.set_linestyle('None')
line2.set_marker('.')
plt.ylabel(r'$\sigma\ (km\ s^{-1})$') 
plt.xlabel(r'$R\ (pc)$')
for item in A:
	plt.loglog(mytable['RADRMS_EXTRAP_DECONV'][item],mytable['VRMS_EXTRAP_DECONV'][item],marker='.',c='b')
plt.loglog(mytable['RADRMS_EXTRAP_DECONV'][92],mytable['VRMS_EXTRAP_DECONV'][92],marker='.',c='r')
plt.loglog(mytable['RADRMS_EXTRAP_DECONV'][132],mytable['VRMS_EXTRAP_DECONV'][132],marker='.',c='m')
plt.tight_layout() 	
plt.savefig('LwRad_matplotlib.png')

#Luminous mass vs. radius plot
figure = plt.figure(figsize=(4.5,4)) #figure size in inches
line1, line2 = plt.loglog(R,M_lum,mytable['RADRMS_EXTRAP_DECONV'],mytable['MASS_EXTRAP'])
line1.set_linestyle('-')
line1.set_color('k')
line2.set_linestyle('None')
line2.set_marker('.')
plt.xlabel(r'$R\ (pc)$') 
plt.ylabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$')
for item in A:
	plt.loglog(mytable['RADRMS_EXTRAP_DECONV'][item],mytable['MASS_EXTRAP'][item],marker='.',c='b')
plt.loglog(mytable['RADRMS_EXTRAP_DECONV'][92],mytable['MASS_EXTRAP'][92],marker='.',c='r')
plt.loglog(mytable['RADRMS_EXTRAP_DECONV'][132],mytable['MASS_EXTRAP'][132],marker='.',c='m')
plt.tight_layout() 	
plt.savefig('MlumRad_matplotlib.png')

#Sigma_0 vs Mass Density
figure = plt.figure(figsize=(4.5,4))
plt.loglog(M_den,sigma0,marker='.',linestyle='None',c='g')
plt.xlabel('$M/\pi R^2\ ((M_{\odot})/pc^2)$')
plt.ylabel('$\sigma_0$')
for item in A:
	plt.loglog(M_den[item],sigma0[item],marker='.',c='b')
plt.loglog(M_den[92],sigma0[92],marker='.',c='r')
plt.loglog(M_den[132],sigma0[132],marker='.',c='m')
plt.tight_layout() 	
plt.savefig('Sigma0_Mden_matplotlib.png')

###Sigma_0 vs Galactocentric Radii
figure = plt.figure(figsize=(4.5,4))
plt.loglog(rgal.to(u.pc),sigma0,marker='.',linestyle='None',c='g')
plt.xlabel('$R_{gal} (pc)$')
plt.ylabel('$\sigma_0$')
for item in A:
	plt.loglog(rgal.to(u.pc)[item],sigma0[item],marker='.',c='b')
plt.loglog(rgal.to(u.pc)[92],sigma0[92],marker='.',c='r')
plt.loglog(rgal.to(u.pc)[132],sigma0[132],marker='.',c='m')
plt.tight_layout() 	
plt.savefig('sigma0_Rgal_matplotlib.png')





	
