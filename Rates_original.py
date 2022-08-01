
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.signal import savgol_filter
from time_conversion import read_units
from astropy import units as u
from astropy.cosmology import Planck15, z_at_value

#---- Important functions

def inspiral_time_peters(a0,e0,m1,m2,af=0):
    """
    Computes the inspiral time, in Gyr, for a binary
    a0 in Au, and masses in solar masses
    
    if different af is given, computes the time from a0,e0
    to that af
    
    for af=0, just returns inspiral time
    for af!=0, returns (t_insp,af,ef)
    """
    coef = 6.086768e-11 #G^3 / c^5 in au, gigayear, solar mass units
    beta = (64./5.) * coef * m1 * m2 * (m1+m2)
    
    if e0 == 0:
        print(e0,a0)
        if not af == 0:
            print("ERROR: doesn't work for circular binaries")
            return 0
        return a0**4 / (4*beta)
    
    c0 = a0 * (1.-e0**2.) * e0**(-12./19.) * (1.+(121./304.)*e0**2.)**(-870./2299.)
    
    if af == 0:
        eFinal = 0.
    else:
        r = ode(deda_peters)
        r.set_integrator('lsoda')
        r.set_initial_value(e0,a0)
        r.integrate(af)
        if not r.successful():
            print("ERROR, Integrator failed!")
        else:
            eFinal = r.y[0]      
    
    time_integrand = lambda e: e**(29./19.)*(1.+(121./304.)*e**2.)**(1181./2299.) / (1.-e**2.)**1.5
    integral,abserr = integrate.quad(time_integrand,eFinal,e0)
    
    if af==0:
        return integral * (12./19.) * c0**4. / beta
    else:
        return (integral * (12./19.) * c0**4. / beta), af, eFinal

def comovingDistance(z):
    h = 0.679
    omegaM = 0.306
    omegaK = 0.
    omegaL = 1 - 0.306
    dh = 3000. / h
    e = lambda zp: 1./np.sqrt(omegaM*(1+zp)**3 + omegaL)
    return dh*integrate.quad(e,0,z)[0]


def lookbackTime(z):
    h = 0.679
    omegaM = 0.306
    omegaK = 0.
    omegaL = 1 - 0.306
    th = 9.78/h
    e = lambda zp: 1./(np.sqrt(omegaM*(1+zp)**3 + omegaL)*(1+zp))
    return th*integrate.quad(e,0,z)[0]

def zAtLookbackTime(t):
    zero = lambda z: lookbackTime(z) - t
    return brentq(zero,0,10)




#---This function gets the merger times for all In-cluster mergers (adapated to my runs, modify as needed)
def get_t_merge(m_min, m_max):
    
    fbh = ['0.50', '0.75', '1.0']
    m_times = []
    m_m1 = []
    m_m2 = []
    for h in range(len(fbh)):
        end = 1
        # end = 15
        if fbh[h] == '0.75':
            end = 25
        for i in range(0,end):
            model = 'FBH_HS1.0_rv1rg20z0.002N8e5_primMS_alpha1_fb0.05_hbf' + fbh[h]+ '_' + str(i) + '_long/'
            path = '/projects/b1095/egp8636/elena_cmc/rundir/triples/' + model

            units=read_units(path+'initial')
            time_units = units[0]['t_myr'] 

            print(model)
            
            flag = 'inactive'
            # if os.path.isfile(path +  'restarted.bhmerger.dat'):
            #     flag, rt, _ = restart_time(path)
               

            if flag == 'inactive':

                bh_merger = np.genfromtxt(path + 'initial.bhmerger.dat')

                if (len(bh_merger))>0:
                    bh_merger_time = bh_merger[:,0]
                    bh_merger_m1 = bh_merger[:,5]
                    bh_merger_m2 = bh_merger[:,6]
                    
                    for j in range(len(bh_merger_time)):
                        if m_max == 40.5:
                            if (bh_merger_m1[j] >= m_min and bh_merger_m1[j] < m_max) and (bh_merger_m2[j] >= m_min and bh_merger_m2[j] < m_max):
                                m_times.append(bh_merger_time[j]*time_units)
                                m_m1.append(bh_merger_m1[j])
                                m_m2.append(bh_merger_m2[j])
                                
                        else:#to consider cases with only one component 
                            
                            if (bh_merger_m1[j] >= m_min and bh_merger_m1[j] < m_max) or (bh_merger_m2[j] >= m_min and bh_merger_m2[j] < m_max):
                                m_times.append(bh_merger_time[j]*time_units)
                                m_m1.append(bh_merger_m1[j])
                                m_m2.append(bh_merger_m2[j])   
                                    
            # elif flag == 'active':
            #     bh_merger_ini = np.genfromtxt(path + 'initial.bhmerger.dat') 
            #     bh_merger_res = np.genfromtxt(path + 'restarted.bhmerger.dat') 
            #     if (len(bh_merger_ini))>0:
            #         bh_merger_time = bh_merger_ini[:,0]
            #         bh_merger_m1 = bh_merger_ini[:,5]
            #         bh_merger_m2 = bh_merger_ini[:,6]
            #         for j in range(len(bh_merger_time)):
            #             if m_max == 40.5:
            #                 if bh_merger_time[j] < rt :
            #                     if (bh_merger_m1[j] >= m_min and bh_merger_m1[j] < m_max) and (bh_merger_m2[j] >= m_min and bh_merger_m2[j] < m_max):
            #                         m_times.append(bh_merger_time[j]*time_units)
            #                         m_m1.append(bh_merger_m1[j])
            #                         m_m2.append(bh_merger_m2[j])
            #             else:#to consider cases with only one component 
            #                 if bh_merger_time[j] < rt :
            #                     if ((bh_merger_m1[j] >= m_min and bh_merger_m1[j] < m_max) or (bh_merger_m2[j] >= m_min and bh_merger_m2[j] < m_max)):
            #                         m_times.append(bh_merger_time[j]*time_units)
            #                         m_m1.append(bh_merger_m1[j])
            #                         m_m2.append(bh_merger_m2[j])

            #     if (len(bh_merger_res))>0:
            #         bh_merger_time = bh_merger_res[:,0]
            #         bh_merger_m1 = bh_merger_res[:,5]
            #         bh_merger_m2 = bh_merger_res[:,6]
            #         for j in range(len(bh_merger_time)):
            #             if m_max == 40.5:
            #                 if bh_merger_time[j]  >= rt: 
            #                     if (bh_merger_m1[j] >= m_min and bh_merger_m1[j] < m_max) and (bh_merger_m2[j] >= m_min and bh_merger_m2[j] < m_max):
            #                         m_times.append(bh_merger_time[j]*time_units)
            #                         m_m1.append(bh_merger_m1[j])
            #                         m_m2.append(bh_merger_m2[j])
            #             else:
            #                 if bh_merger_time[j] >= rt:
            #                     if ((bh_merger_m1[j] >= m_min and bh_merger_m1[j] < m_max) or (bh_merger_m2[j] >= m_min and bh_merger_m2[j] < m_max)):
            #                         m_times.append(bh_merger_time[j]*time_units)
            #                         m_m1.append(bh_merger_m1[j])
            #                         m_m2.append(bh_merger_m2[j])
            
    return m_times, m_m1, m_m2
            





#---This function gets the merger times for all out-cluster mergers (adapated to my runs, modify as needed)

def out_t_merge(m_min, m_max):
    from ipynb.fs.full.LISA_calculations import inspiral_time_peters
    times = []
    m1 = []
    m2 = []
    fbh = ['0.50', '0.75', '1.0']

    for hmbf in fbh:
        end = 1
        # end = 15
        if hmbf == '0.75':
            end = 25
        for e in range(0,end):
            model = 'FBH_HS1.0_rv1rg20z0.002N8e5_primMS_alpha1_fb0.05_hbf' + hmbf + '_' + str(e) 
            path = '/projects/b1095/egp8636/elena_cmc/rundir/triples/'+model+'_long/'
            print(model)
            units=read_units(path+'initial')
            time_units = units[0]['t_myr'] 
            

            flag ='inactive'
            # if os.path.isfile(path +  'restarted.esc.dat'):
            #         flag, rt, _ = restart_time(path)

            if flag == 'inactive':

                data = np.genfromtxt(path + 'initial.esc.dat')
                binflag = data[:,14]

                for j in range(len(binflag)):
                    if binflag[j] == 1:
                        t = data[j,1]
                        t = t
                        type1= int(data[j,22])
                        type2 = int(data[j,23])
                        mass1 =  (data[j,15])
                        mass2 =  (data[j,16])
                        a = (data[j,19])
                        e = (data[j,20])
                        if  (type1 == 14 and type2 == 14):
                            if  m_max == 40.5:
                                if (mass1 >= m_min and mass1 < m_max) and (mass2 >= m_min and mass2 < m_max):
                                    
                                    inspiral_time = inspiral_time_peters(a,e,mass1,mass2)
                                    t_tot = t*time_units + inspiral_time*10**3
                                    print(t,inspiral_time*10**3,t_tot)
                                    if t_tot  <= 13700.0:
                        
                                            times.append(t_tot)
                            else:
                                if (mass1 >= m_min and mass1 < m_max) or (mass2 >= m_min and mass2 < m_max):
                                
                                    inspiral_time = inspiral_time_peters(a,e,mass1,mass2)
                                    t_tot = t*time_units + inspiral_time*10**3
                                    if t_tot<= 13700.0:
                                            times.append(t_tot)
            # if flag == 'active':

            #     data_ini = np.genfromtxt(path + 'initial.esc.dat')
            #     binflag = data_ini[:,14]

            #     for j in range(len(binflag)):
            #         if binflag[j] == 1:
            #             t = data_ini[j,1]
            #             type1= int(data_ini[j,22])
            #             type2 = int(data_ini[j,23])
            #             mass1 =  (data_ini[j,15])
            #             mass2 =  (data_ini[j,16])
            #             a = (data_ini[j,19])
            #             e = (data_ini[j,20])
            #             if  (type1 == 14 and type2 == 14):
                            
            #                 if  m_max == 40.5:
            #                     if t < rt :
                                    
            #                         if (mass1 >= m_min and mass1 < m_max) and (mass2 >= m_min and mass2 < m_max):
            #                             inspiral_time = inspiral_time_peters(a,e,mass1,mass2)
            #                             t_tot = t*time_units + inspiral_time*10**3
            #                             if t_tot <= 13700.0:
            #                                 times.append(t_tot)
            #                 else:

            #                     if t < rt :
            #                         if (mass1 >= m_min and mass1 < m_max) or (mass2 >= m_min and mass2 < m_max):
                                
            #                             inspiral_time = inspiral_time_peters(a,e,mass1,mass2)
            #                             t_tot = t*time_units + inspiral_time*10**3
            #                             if t_tot<= 13700.0:
            #                                 times.append(t_tot)
                                
            #     data_res = np.genfromtxt(path + 'restarted.esc.dat')
            #     binflag = data_res[:,14]

            #     for j in range(len(binflag)):
            #         if binflag[j] == 1:
                        
            #             t = data_res[j,1]
            #             t = t
            #             type1= int(data_res[j,22])
            #             type2 = int(data_res[j,23])
            #             mass1 =  (data_res[j,15])
            #             mass2 =  (data_res[j,16])
            #             a = (data_res[j,19])
            #             e = (data_res[j,20])
            #             if  (type1 == 14 and type2 == 14):
            #                 if  m_max == 40.5:
            #                     if t >= rt:
            #                         if (mass1 >= m_min and mass1 < m_max) and (mass2 >= m_min and mass2 < m_max):

            #                             inspiral_time = inspiral_time_peters(a,e,mass1,mass2)
            #                             t_tot = t*time_units + inspiral_time*10**3
            #                             if t_tot <= 13700.0:
            #                                 times.append(t_tot)
            #                 else:
            #                     if t  >= rt:
            #                         if (mass1 >= m_min and mass1 < m_max) or (mass2 >= m_min and mass2 < m_max):

            #                             inspiral_time = inspiral_time_peters(a,e,mass1,mass2)
            #                             t_tot = t*time_units + inspiral_time*10**3
            #                             if t_tot <= 13700.0:
                                            
                                            
            #                                 times.append(t_tot)

    return times


in_cluster_merger_times=open("TEST_in_cluster_merger_times.txt",'w')
m_times, m_m1, m_m2 = get_t_merge(0, 100000000)
for m in m_times:
    m = str(m)
    in_cluster_merger_times.write(f'{m}\n')
in_cluster_merger_times.close()

# Write out-of-cluster merger times to a separate .txt file                    ## NEED TO FIGURE OUT STRUCTURE OF VAR 'times'
out_cluster_merger_times=open("TEST_out_cluster_merger_times.txt",'w')
times = get_t_merge(0, 100000000)

temp_count = 0
print(f'type(times) = {type(times)}')
print(f'len(times) = {len(times)}')

for t in times:
    # if temp_count == 0:
    print(f'type(t) = {type(t)}')
    print(f'len(t) = {len(t)}')
    # temp_count += 1
    t = str(t)
    out_cluster_merger_times.write(f'{t}\n')
out_cluster_merger_times.close()

files = ["TEST_in_cluster_merger_times.txt", "TEST_out_cluster_merger_times.txt"]    
    


#--- I saved all of the merger times to data files. I deleted this part not to confuse you.

#--- This gets the cluster t_ages from an El-Badry paper

gc_ages = np.genfromtxt('Mvir1e14.txt')
gc_age  = gc_ages[:,0]
gc_met  = gc_ages[:,1]

#get only age values with metallicities between -0.8 and -1.2 
gc_age_met = [gc_age[i]*10**3 for i in range(len(gc_age)) if gc_met[i] >= -1.3 and gc_met[i] <= -0.7]

#create files of effective times
for file in files:
    t_mergers = np.genfromtxt('merger_times/' + file)
    f1 = open('effective_times/eff_' + file, 'w')
    for t_merger in t_mergers:
        
        #for each merger, draw 100 cluster ages 

        t_ages = np.random.choice(gc_age_met, size = 100)

        #compute the effective merger time
        
        for tage in t_ages:
            t_eff = 13.7 * 10**3 - tage + t_merger #in Myr
            f1.write(str(t_eff))
            f1.write('\n')
    f1.close()




#--- This function calculates the cumulative merger rate and volumetric merger rates.

def rate(t_effs):
    
    '''This is used to calculate the merger rate
    
    PARAMETERS
    ---------------

    t_effs : list of effective merger times in Myr
    
    
    OUTPUT
    ---------------
    rates: rates in units of merger/Gpc^3/yr
    '''

    n = 55 #number of models (personal to my project)
    m = 100
    rates = []
    d_mins = []
    d_maxs = []
    zs = []
    

    #divide the list of effective times by into separate redshift bins

    t_bins = np.linspace(10,13200, num = 300) #300 bins in Myr
    for j in range(1,len(t_bins)-1):
        
        t_lower = t_bins[j]
        t_upper = t_bins[j+1]

        count  = 0 
        for teff in t_effs:
            if teff >= t_lower and teff < t_upper:
                count += 1
                
        
        delta_t = (t_upper-t_lower)*10**6 #in years
        rho = 2.31*10**9 #Gpc^-3
        weight = (m*n)
        R = ((count*rho)/weight)/(delta_t)*4.65
        rates.append(R)
        
        t_lb1 = (13.3*10**9- t_lower*10**6)*u.yr #lookback time
        t_lb2 = (13.3*10**9- t_upper*10**6)*u.yr
        t_mid = (13.3*10**9- (t_upper*10**6+t_lower*10**6)/2)*u.yr
        
        zmin = z_at_value(Planck15.lookback_time,t_lb1)
        zmax = z_at_value(Planck15.lookback_time,t_lb2)
        zmid = z_at_value(Planck15.lookback_time,t_mid)
        
        d_min = (Planck15.comoving_distance(zmin)).to(u.Gpc)
        d_max = (Planck15.comoving_distance(zmax)).to(u.Gpc)
        
        zs.append(zmid.value)
        d_mins.append(d_min.value)
        d_maxs.append(d_max.value)
        
    #calculating cumulative rates
    rc = cumulative_rate(rates, zs, d_maxs,d_mins)
    return rates, rc, zs



def cumulative_rate(rates,zmid, d_max,d_min):
    cumulative_rates = []

    r_c = 0
    
    for i in range(len(rates)-1,-1,-1):
        r = rates[i]
        
        r_c += r * (1 + zmid[i])**(-1) * (4/3) * np.pi * abs((d_max[i]**3 - d_min[i]**3))

        cumulative_rates.append(r_c)


    return cumulative_rates


