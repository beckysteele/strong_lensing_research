import numpy as np
import random
from astropy.coordinates import SkyCoord
from astropy import units as u
from copy import deepcopy
from astropy import constants as const
from astropy import units as u
from astropy import cosmology

# Becky comments

datadir = '/disks/strw9/BO_XUV/Catalogue/'

Mag_sun_w1 = 3.254
Mag_sun_error_w1 = np.sqrt(0.008**2+0.02**2) #extra 0!!
Mag_sun_kt = 3.302
Mag_sun_error_kt = np.sqrt(0.008**2+0.02**2)
lum_sun = 3.828e26 #in W from NASA fact sheet
c = const.c.to('km/s').value #c in km/s
nu_3 = (c*1000)/(11.5608e-6) #in s^-1 Jarrett11
nu_3_error = 0.0446e-6*(c*1000)/(nu_3**2) #in s^-1 Jarrett11
nu_4 = (c*1000)/(22.0883e-6) #in s^-1 Jarrett11
nu_4_error = 0.1184e-6*(c*1000)/(nu_4**2) #in s^-1 Jarrett11
zp_3 = 31.674 #in Jy Jarrett11
zp_3_error = 0.450 #in Jy Jarrett11
zp_4 = 8.363 #in Jy Jarrett11
zp_4_error = 0.124 #in Jy Jarrett11
a_3 = [1.13, 10.24] #Cluver14
a_4 = [0.82, 7.3] #Cluver14

h0 = 73.8 #h0 in km/(Mpc*s) Riess2011
h0_error = 2.4 # km/(Mpc*s)
parsec = 3.08567758e16 #m
jansky = 1e-26 #W/m^2/Hz

# In order to use this function, need to calculate cz and cz_error and pass as args
# Use redshift as a parameter directly instead of doing cz stuff.
def velocity(cz, cz_error):
	z = cz/c
	v = c * ((1+z)**2 -1)/((1+z)**2+1) #v in km/s
	v = abs(v)
	z_error = cz_error/c
	v_error = z_error*c*4*(1+z)/(((1+z)**2+1)**2) #v_error in km/s
	#print 'v: ', v
	#print 'verr: ', v_error
	return v, v_error #km/s

# In order to use this function, need to execute velocity() first	
def distance(v, v_error):
	D = v/h0 #D in Mpc
	D = D * 1e6 #D in pc
	D_error = np.sqrt((v_error/h0)**2+(h0_error*(-v)/(h0**2))**2) *1e6 #D_error in pc
	#print 'D: ', D
	#print 'Derr: ', D_error
	return D, D_error #pc

# Absolute magnitude calculation. In order to use this function, need to do distance() and velocity() first.
def Mag(D, D_error, m, m_error):
	M = m-5*np.log10(D) + 5
	M_error = np.sqrt((m_error)**2+(D_error*(-5)/(D*np.log(10)))**2)
	#print 'M: ', M
	#print 'Merr: ', M_error
	return M, M_error

# Luminosity calculation. Need to run Mag(), distance(), velocity() first.
def luminosity(M, M_error, Mag_sun, Mag_sun_error):
	L = 10**(0.4*(Mag_sun - M)) #L in solar luminosities
	a = np.log(10**0.4)
	L_error = np.sqrt((L*a*M_error)**2+(Mag_sun_error*L*a)**2) #L_error in solar luminosities
	#print 'L: ', L
	#print 'Lerr: ', L_error
	return L, L_error

# Stellar mass estimation based on 0.6*L.	
def mass(L, L_error):
	M = L*0.6 #M in solar mass
	M_error = 0.6*L_error
	#print 'Mass: ', M
	#print 'Masserr: ', M_error
	return M, M_error

# Calculated mass based on calculated absolute and apparent mag
def calculate_mass(cz, cz_error, m, m_error, Mag_sun, Mag_sun_error):
	v, v_error = velocity(cz, cz_error)
	D, D_error = distance(v, v_error)
	M, M_error = Mag(D, D_error, m, m_error)
	L, L_error = luminosity(M, M_error, Mag_sun, Mag_sun_error)
	Mass, Mass_error = mass(L, L_error)
	#M, M_error = mass(luminosity(Mag((distance(*velocity(cz))), m, m_error)))
	return Mass, Mass_error

# Star formation rate based on cz, m, nu, zp, a
def SFR(cz, cz_error, m, m_error, nu, nu_error, zp, zp_error, a):
	v, v_error = velocity(cz, cz_error) #km/s
	
	D, D_error = distance(v, v_error) #pc
	D = D * parsec #convert D to meters
	D_error = D_error * parsec #convert D_error to meters
	
	M, M_error = Mag(D, D_error, m, m_error)
	
	fluxdens = zp*10**(-0.4*m) #in Jy = 10^-26 W/m^2/Hz
	fluxdens_error = np.sqrt((m_error * np.log(10**-0.4) * fluxdens)**2 + (fluxdens/zp *zp_error)**2) #in Jy
	fluxdens = fluxdens * jansky #convert fluxdens to W/m^2/Hz
	fluxdens_error = fluxdens_error * jansky #convert fluxdens_error to W/m^2/Hz
	
	lumdens = 4*np.pi*D**2*fluxdens*nu #W
	lumdens_error = np.sqrt((D_error*2*lumdens/D)**2 + (fluxdens_error*lumdens/fluxdens)**2 + (nu_error * lumdens/nu)**2) #W
	lumdens = lumdens / lum_sun #in solar luminosity
	lumdens_error = lumdens_error / lum_sun #in solar luminosity
	
	SFR = 10**(a[0]*np.log10(lumdens)-a[1]) #in M_sun/year
	SFR_error = lumdens_error*a[0]*10**(-a[1])*lumdens**(a[0]-1)
	return SFR, SFR_error

# Estimated SFR based on two WISE bands, wa - wb
def SFR_barbaric(wa, wa_error, wb, wb_error):
	wa_wb = wa - wb
	wa_wb_error = np.sqrt(wa_error**2 + wb_error**2)
	return wa_wb, wa_wb_error

# Radius estimate	
def radius(D, angle): 
	rad_angle = (angle/60)*(np.pi/180) #arcmin to radians
	return D*np.tan(rad_angle/2) #D in pc -> radius in pc #no error

# This returns the SFR per square area and its error	
def sig(dist, major, minor, SFR, e_SFR):
	opp = np.pi*radius(dist, major)*radius(dist, minor) #pc^2
	#print radius(dist,major)*3.26163344 #to ly
	sig = SFR/opp #Msun/pc^2/yr
	e_sig = 1/opp*e_SFR
	print sig
	return sig, e_sig	

# Writes calcs to files
total = np.loadtxt(datadir + 'total_2mass_wise_huchra_calc_new2.txt', skiprows = 7, dtype='string')
lem = np.loadtxt(datadir + 'total_lemonias_cor_calc_new3.txt', skiprows = 8, dtype='string')
thi = np.loadtxt(datadir + 'total_thilker_cor_calc_new3.txt', skiprows = 8, dtype='string')
meu = np.loadtxt(datadir + 'total_meurer_cor_calc_new3.txt', skiprows = 8, dtype='string')
tot_randsamp = np.loadtxt(datadir + 'total_2mass_wise_huchra_calc_new3_randsamp_backup.txt', skiprows = 7, dtype='string')
#np.savetxt(datadir+'/namesradii/lem_names.txt', lem[:,0], fmt='%s')
#np.savetxt(datadir+'/namesradii/thi_names.txt', thi[:,0], fmt='%s')
#np.savetxt(datadir+'/namesradii/meu_names.txt', meu[:,0], fmt='%s')

#Creating csv files with comma's as delimiters:
np.savetxt(datadir+'total_2mass_wise_huchra_calc_new2_csv.csv', total, fmt='%s', delimiter=',')
np.savetxt(datadir+'total_lemonias_cor_calc_new3_csv.csv', lem, fmt='%s', delimiter=',')
np.savetxt(datadir+'total_thilker_cor_calc_new3_csv.csv', thi, fmt='%s', delimiter=',')
np.savetxt(datadir+'total_meurer_cor_calc_new3_csv.csv', meu, fmt='%s', delimiter=',')
np.savetxt(datadir+'total_2mass_wise_huchra_calc_new3_randsamp_backup_csv.csv', tot_randsamp, fmt='%s', delimiter=',')

#Create random sample of 2MASS galaxies, 2500 galaxies! (250 per keer)
"""
randomarray = np.zeros(2500)
for i in range(len(randomarray)):
	randomarray[i] = random.randint(0,len(total[:,0]))
randomarray = np.sort(randomarray)
print 'random array:', randomarray
randomtwomass = np.empty((2500, len(total[0,:])), dtype='S60')
for i in range(len(randomarray)):
	randomtwomass[i,range(0,len(total[0,:]))] = total[randomarray[i], range(0,len(total[0,:]))]
print randomtwomass
randomnames = np.empty(len(randomtwomass[:,0]), dtype='S50')
for i in range(len(randomnames)):
	name = randomtwomass[i,0]
	name1 = name[0:6]
	name2 = name[6:15]
	name3 = name[15]
	randomnames[i] = "2MASXJ"+name1+'.'+name2+'.'+name3
np.savetxt(datadir+'/namesradii/total_names.txt', randomnames, fmt='%s')
np.savetxt(datadir+'total_2mass_wise_huchra_calc_new2_randsamp.txt', randomtwomass, fmt='%s', delimiter = '	')
"""

lem_radii = np.loadtxt(datadir+'namesradii/lem_radii.txt', usecols=(1,2), skiprows=3) #1st col: major, 2nd col: minor
thi_radii = np.loadtxt(datadir+'namesradii/thi_radii.txt', usecols=(1,2), skiprows=3) #1st col: major, 2nd col: minor
meu_radii = np.loadtxt(datadir+'namesradii/meu_radii.txt', usecols=(1,2), skiprows=3) #1st col: major, 2nd col: minor
tot_radii = np.loadtxt(datadir+'namesradii/total_radii.txt', usecols=(1,2), skiprows=3) #1st col: major, 2nd col: minor
#print tot_radii
#sig(dist, major, minor, SFR, e_SFR)
sig_lem, e_sig_lem = sig(lem[:,45].astype('float'), lem_radii[:,0], lem_radii[:,1], lem[:,57].astype('float'), lem[:,58].astype('float'))
sig_thi, e_sig_thi = sig(thi[:,45].astype('float'), thi_radii[:,0], thi_radii[:,1], thi[:,57].astype('float'), thi[:,58].astype('float'))
sig_meu, e_sig_meu = sig(meu[:,50].astype('float'), meu_radii[:,0], meu_radii[:,1], meu[:,62].astype('float'), meu[:,63].astype('float'))
sig_tot, e_sig_tot = sig(tot_randsamp[:,39].astype('float'), tot_radii[:,0], tot_radii[:,1], tot_randsamp[:,51].astype('float'), tot_randsamp[:,52].astype('float'))

#Adding the star formation rate surface density
"""
lem_plus = np.empty((len(lem[:,0]),len(lem[0,:])+2), dtype='S60')
lem_plus[:,range(0,len(lem[0,:]))] = lem[:, range(0,len(lem[0,:]))]
lem_plus[:,len(lem[0,:])] = sig_lem
lem_plus[:,len(lem[0,:])+1] = e_sig_lem
#np.savetxt(datadir+'total_lemonias_cor_calc_new3.txt', lem_plus, fmt='%s', delimiter = '    ')
thi_plus = np.empty((len(thi[:,0]),len(thi[0,:])+2), dtype='S60')
thi_plus[:,range(0,len(thi[0,:]))] = thi[:, range(0,len(thi[0,:]))]
thi_plus[:,len(thi[0,:])] = sig_thi
thi_plus[:,len(thi[0,:])+1] = e_sig_thi
#np.savetxt(datadir+'total_thilker_cor_calc_new3.txt', thi_plus, fmt='%s', delimiter = '    ')
meu_plus = np.empty((len(meu[:,0]),len(meu[0,:])+2), dtype='S60')
meu_plus[:,range(0,len(meu[0,:]))] = meu[:, range(0,len(meu[0,:]))]
meu_plus[:,len(meu[0,:])] = sig_meu
meu_plus[:,len(meu[0,:])+1] = e_sig_meu
#np.savetxt(datadir+'total_meurer_cor_calc_new3.txt', meu_plus, fmt='%s', delimiter = '    ')
tot_randsamp_plus = np.empty((len(tot_randsamp[:,0]),len(tot_randsamp[0,:])+2), dtype='S60')
tot_randsamp_plus[:,range(0,len(tot_randsamp[0,:]))] = tot_randsamp[:, range(0,len(tot_randsamp[0,:]))]
tot_randsamp_plus[:,len(tot_randsamp[0,:])] = sig_tot
tot_randsamp_plus[:,len(tot_randsamp[0,:])+1] = e_sig_tot
#np.savetxt(datadir+'total_2mass_wise_huchra_calc_new3_randsamp.txt', tot_randsamp_plus, fmt='%s', delimiter = '    ')
"""



"""
cz = total[:,15].astype('float')
cz_error = total[:,16].astype('float')
m_kt = total[:,6].astype('float')
m_kt_error = total[:,12].astype('float')
w1 = total[:,21].astype('float')
w1_error = total[:,22].astype('float')
w2 = total[:,25].astype('float')
w2_error = total[:,26].astype('float')
w3 = total[:,29].astype('float')
w3_error = total[:,30].astype('float')
w4 = total[:,33].astype('float')
w4_error = total[:,34].astype('float')

print cz
v, v_error = velocity(cz, cz_error)
dist, dist_error = distance(v, v_error)
mass_kt, mass_kt_error = calculate_mass(cz, cz_error, m_kt, m_kt_error, Mag_sun_kt, Mag_sun_error_kt)
mass_w1, mass_w1_error = calculate_mass(cz, cz_error, w1, w1_error, Mag_sun_w1, Mag_sun_error_w1)
w1_w3, w1_w3_error = SFR_barbaric(w1, w1_error, w3, w3_error)
w1_w2, w1_w2_error = SFR_barbaric(w1, w1_error, w2, w2_error)
w2_w3, w2_w3_error = SFR_barbaric(w2, w2_error, w3, w3_error)
SFR_w3, SFR_w3_error = SFR(cz, cz_error, w3, w3_error, nu_3, nu_3_error, zp_3, zp_3_error, a_3)
SFR_w4, SFR_w4_error = SFR(cz, cz_error, w4, w4_error, nu_4, nu_4_error, zp_4, zp_4_error, a_4)


total_plus = np.empty((len(total[:,0]),55), dtype='S60')
total_plus[:,range(0,37)] = total[:,range(0,37)]
total_plus[:,37] = v
total_plus[:,38] = v_error
total_plus[:,39] = dist
total_plus[:,40] = dist_error
total_plus[:,41] = mass_kt
total_plus[:,42] = mass_kt_error
total_plus[:,43] = mass_w1
total_plus[:,44] = mass_w1_error
total_plus[:,45] = w1_w3
total_plus[:,46] = w1_w3_error
total_plus[:,47] = w1_w2
total_plus[:,48] = w1_w2_error
total_plus[:,49] = w2_w3
total_plus[:,50] = w2_w3_error
total_plus[:,51] = SFR_w3
total_plus[:,52] = SFR_w3_error
total_plus[:,53] = SFR_w4
total_plus[:,54] = SFR_w4_error

print total_plus[0,:]

np.savetxt(datadir+'total_2mass_wise_huchra_calc_new2.txt', total_plus, fmt='%s', delimiter = '    ')

cz = thi[:,21].astype('float')
cz_error = thi[:,22].astype('float')
m_kt = thi[:,12].astype('float')
m_kt_error = thi[:,18].astype('float')
w1 = thi[:,27].astype('float')
w1_error = thi[:,28].astype('float')
w2 = thi[:,31].astype('float')
w2_error = thi[:,32].astype('float')
w3 = thi[:,35].astype('float')
w3_error = thi[:,36].astype('float')
w4 = thi[:,39].astype('float')
w4_error = thi[:,40].astype('float')

print cz
v, v_error = velocity(cz, cz_error)
dist, dist_error = distance(v, v_error)
mass_kt, mass_kt_error = calculate_mass(cz, cz_error, m_kt, m_kt_error, Mag_sun_kt, Mag_sun_error_kt)
mass_w1, mass_w1_error = calculate_mass(cz, cz_error, w1, w1_error, Mag_sun_w1, Mag_sun_error_w1)
w1_w3, w1_w3_error = SFR_barbaric(w1, w1_error, w3, w3_error)
w1_w2, w1_w2_error = SFR_barbaric(w1, w1_error, w2, w2_error)
w2_w3, w2_w3_error = SFR_barbaric(w2, w2_error, w3, w3_error)
SFR_w3, SFR_w3_error = SFR(cz, cz_error, w3, w3_error, nu_3, nu_3_error, zp_3, zp_3_error, a_3)
SFR_w4, SFR_w4_error = SFR(cz, cz_error, w4, w4_error, nu_4, nu_4_error, zp_4, zp_4_error, a_4)


total_plus = np.empty((len(thi[:,0]),61), dtype='S60')
total_plus[:,range(0,43)] = thi[:,range(0,43)]
total_plus[:,43] = v
total_plus[:,44] = v_error
total_plus[:,45] = dist
total_plus[:,46] = dist_error
total_plus[:,47] = mass_kt
total_plus[:,48] = mass_kt_error
total_plus[:,49] = mass_w1
total_plus[:,50] = mass_w1_error
total_plus[:,51] = w1_w3
total_plus[:,52] = w1_w3_error
total_plus[:,53] = w1_w2
total_plus[:,54] = w1_w2_error
total_plus[:,55] = w2_w3
total_plus[:,56] = w2_w3_error
total_plus[:,57] = SFR_w3
total_plus[:,58] = SFR_w3_error
total_plus[:,59] = SFR_w4
total_plus[:,60] = SFR_w4_error
"""


#print total_plus[0,:]

#np.savetxt(datadir+'total_thilker_cor_calc_new2.txt', total_plus, fmt='%s', delimiter = '    ')

print 'done'

