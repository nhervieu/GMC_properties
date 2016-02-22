#import libraries
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from galaxies import Galaxy
import astropy.units as u
import powerlaw

#load file
mytable = Table.read('m100.co10.kkms_props_cprops.fits')

##VARIABLE DEFINITIONS##

# Pull out the mass variable into a numpy array.
mass = mytable['MASS_EXTRAP'].data 
#remove two outliers    
mass = np.delete(mass,132)
mass = np.delete(mass,92)
 
#Radius vs line width best fit line
R = np.arange(10,1000)
S = np.power(np.pi,1/2)*R/3.4
sigmav = np.sqrt(S)
one = np.arange(10000,10000000000,10000)

#Virial Mass best fit line
M_vir = 540*np.power(R,2)

#Luminous Mass best fit line
M_lum = 2000*np.power(np.power(M_vir/39,1.234567901)/130,0.8)

#Sigma_0, Mass Density, R_gal
sigma0 = mytable['VRMS_EXTRAP_DECONV']/np.sqrt(mytable['RADRMS_EXTRAP_DECONV'])  
M_den = mytable['MASS_EXTRAP']/(np.pi*np.power(mytable['RADRMS_EXTRAP_DECONV'],2)) 

mygalaxy = Galaxy("M100")
cpropstable = Table.read('m100.co10.kkms_props_cprops.fits')
rgal=mygalaxy.radius(ra = cpropstable['XPOS'], dec = cpropstable['YPOS'])

#indexes for low R_gal group 
index = np.where(rgal.value < 1100)	
	
##PLOTS##

#Virial mass vs.luminous mass plot
figure = plt.figure(figsize=(4.5,4)) #figure size in inches
line1, line2 = plt.loglog(one,one, mytable['MASS_EXTRAP'],mytable['VIRMASS_EXTRAP_DECONV'])
line1.set_linestyle('-')
line1.set_color('k')
line2.set_linestyle('None')
line2.set_marker('.')
plt.loglog(mytable['MASS_EXTRAP'][index],mytable['VIRMASS_EXTRAP_DECONV'][index],marker='.',c='b',linestyle='None')
plt.loglog(mytable['MASS_EXTRAP'][92],mytable['VIRMASS_EXTRAP_DECONV'][92],marker='.',c='w')
plt.loglog(mytable['MASS_EXTRAP'][132],mytable['VIRMASS_EXTRAP_DECONV'][132],marker='.',c='w')
plt.xlabel(r'$log(M_{\mathrm{lum}})\ (M_{\odot})$') 
plt.ylabel(r'$log(M_{\mathrm{vir}})\ (M_{\odot})$')
plt.tight_layout() 	
plt.savefig('MlumMvir_matplotlib.png')

#Radius vs. line width(velocity dispersion) plot
figure = plt.figure(figsize=(4.5,4)) #figure size in inches
line1, line2 = plt.loglog(R,sigmav, mytable['RADRMS_EXTRAP_DECONV'],mytable['VRMS_EXTRAP_DECONV'])
line1.set_linestyle('-')
line1.set_color('k')
line2.set_linestyle('None')
line2.set_marker('.')
plt.ylabel(r'$log(\sigma)\ (km\ s^{-1})$') 
plt.xlabel(r'$log(R)\ (pc)$')
plt.loglog(mytable['RADRMS_EXTRAP_DECONV'][index],mytable['VRMS_EXTRAP_DECONV'][index],marker='.',c='b',linestyle='None')
plt.loglog(mytable['RADRMS_EXTRAP_DECONV'][92],mytable['VRMS_EXTRAP_DECONV'][92],marker='.',c='w')
plt.loglog(mytable['RADRMS_EXTRAP_DECONV'][132],mytable['VRMS_EXTRAP_DECONV'][132],marker='.',c='w')
plt.tight_layout() 	
plt.savefig('LwRad_matplotlib.png')

#Luminous mass vs. radius plot
figure = plt.figure(figsize=(4.5,4)) #figure size in inches
line1, line2 = plt.loglog(R,M_lum,mytable['RADRMS_EXTRAP_DECONV'],mytable['MASS_EXTRAP'])
line1.set_linestyle('-')
line1.set_color('k')
line2.set_linestyle('None')
line2.set_marker('.')
plt.xlabel(r'$log(R)\ (pc)$') 
plt.ylabel(r'$log(M_{\mathrm{lum}})\ (M_{\odot})$')
plt.loglog(mytable['RADRMS_EXTRAP_DECONV'][index],mytable['MASS_EXTRAP'][index],marker='.',c='b',linestyle='None')
plt.loglog(mytable['RADRMS_EXTRAP_DECONV'][92],mytable['MASS_EXTRAP'][92],marker='.',c='w')
plt.loglog(mytable['RADRMS_EXTRAP_DECONV'][132],mytable['MASS_EXTRAP'][132],marker='.',c='w')
plt.tight_layout() 	
plt.savefig('MlumRad_matplotlib.png')

#Sigma_0 vs Mass Density
figure = plt.figure(figsize=(4.5,4))
plt.loglog(M_den,sigma0,marker='.',linestyle='None',c='g')
#plt.xlabel('$M/\pi R^2\ ((M_{\odot})/pc^2)$')
plt.xlabel('$log(\Sigma)\ ((M_{\odot})/pc^2)$')
plt.ylabel('$log(\sigma_0)\ (km\ s^{-1})$')
plt.loglog(M_den[index],sigma0[index],marker='.',c='b',linestyle='None')
plt.loglog(M_den[92],sigma0[92],marker='.',c='w')
plt.loglog(M_den[132],sigma0[132],marker='.',c='w')
plt.tight_layout() 	
plt.savefig('Sigma0_Mden_matplotlib.png')

#Sigma_0 vs Galactocentric Radii
figure = plt.figure(figsize=(4.5,4))
plt.loglog(rgal.to(u.pc),sigma0,marker='.',linestyle='None',c='g')
plt.xlabel('$log(R_{gal})\ (pc)$')
plt.ylabel('$log(\sigma_0)\ (km\ s^{-1})$')
plt.loglog(rgal.to(u.pc)[index],sigma0[index],marker='.',c='b',linestyle='None')
plt.loglog(rgal.to(u.pc)[92],sigma0[92],marker='.',c='w')
plt.loglog(rgal.to(u.pc)[132],sigma0[132],marker='.',c='w')
plt.tight_layout() 	
plt.savefig('sigma0_Rgal_matplotlib.png')

#X,Y position
figure = plt.figure(figsize=(4.5,4)) #figure size in inches
plt.plot(mytable['XPOS'],mytable['YPOS'],linestyle = 'None', marker = '.',c = 'g')
plt.xlabel('X Position') 
plt.ylabel('Y Position')
plt.plot(mytable['XPOS'][index],mytable['YPOS'][index],marker='.',c='b',linestyle='None')
plt.plot(mytable['XPOS'][92],mytable['YPOS'][92],marker='.',c='w')
plt.plot(mytable['XPOS'][132],mytable['YPOS'][132],marker='.',c='w')
plt.tight_layout() 	
plt.savefig('xypos_matplotlib.png')


## MASS DISTRIBUTIONS ##

#Mass Distribution for All Clouds
figure = plt.figure(figsize=(4.5,4)) 
#Fit data 
myfit = powerlaw.Fit(mass)
myfit.plot_ccdf(label='Fit')
R, p = myfit.distribution_compare('power_law','truncated_power_law')
#Plot
myfit.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
myfit.power_law.plot_ccdf(label='Power Law')
myfit.plot_ccdf(drawstyle='steps',label='Data')
plt.legend(loc ='lower left',prop={'size':8})
plt.ylabel(r'$N$')
plt.xlabel(r'$Mass\ (M_{\odot})$')
plt.tight_layout() 	
plt.savefig('powerlaw.png')


#Mass Distribution for Nuclear Clouds
figure = plt.figure(figsize=(4.5,4)) 
#Keep only clouds with rgal <1kpc
mass_nuc = mass[index]
myfit_nuc = powerlaw.Fit(mass_nuc,xmin = myfit.xmin)
R_nuc, p_nuc = myfit_nuc.distribution_compare('power_law','truncated_power_law')
#Plot
myfit_nuc.plot_ccdf(label='Fit')
myfit_nuc.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
myfit_nuc.power_law.plot_ccdf(label='Power Law')
myfit_nuc.plot_ccdf(drawstyle='steps',label='Data (Nuclear Clouds)')
plt.legend(loc ='lower left',prop={'size':8})
plt.ylabel(r'$N$')
plt.xlabel(r'$Mass\ (M_{\odot})$')
plt.tight_layout() 	
plt.savefig('powerlaw_nuc.png')


#Mass Distribution for Disk Clouds
figure = plt.figure(figsize=(4.5,4)) 
#Ignore all clouds with rgal >1kpc
mass_disk = np.delete(mass,index)
#Fit Data
myfit_disk = powerlaw.Fit(mass_disk,xmin = myfit.xmin)
R_disk, p_disk = myfit_disk.distribution_compare('power_law','truncated_power_law')
#Plot
myfit_disk.plot_ccdf(label='Fit')
myfit_disk.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
myfit_disk.power_law.plot_ccdf(label='Power Law')
myfit_disk.plot_ccdf(drawstyle='steps',label='Data (Disk Clouds)')
plt.legend(loc ='lower left',prop={'size':8})
plt.ylabel(r'$N$')
plt.xlabel(r'$Mass\ (M_{\odot})$')
plt.tight_layout() 	
plt.savefig('powerlaw_disk.png')


#print out table of alpha, R and p values for each mass distribution
tb = {'Clouds': ['All','Nuclear','Disk'],'Alpha': [myfit.alpha,myfit_nuc.alpha,myfit_disk.alpha],'R':[R,R_nuc,R_disk], 'p':[p,p_nuc,p_disk]}
print Table(tb,names =('Clouds','Alpha','R','p'))








	
