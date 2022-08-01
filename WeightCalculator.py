import kick as kc
import detection as dt
import numpy as np
import os
import sys
from scipy.interpolate import interp1d
from scipy.stats import lognorm
import astropy.units as u
from astropy.cosmology import Planck15, z_at_value

import random
from datetime import datetime
random.seed(datetime.now())

# ---------- constants (units AU=1, Msun=1, G=1)

pi = 3.14159
pc = 2.0*10**5      # pc in AU
vsc = 30.0          # scale of velocity 30 km/s
c = 10 ** 4           # velocity of light

dyear = 365.24  	  # from yr to day

#---------- calculations

#appr = 'IMRPhenomC'
appr = 'IMRPhenomPv2'
psd = 'aLIGOZeroDetHighPower'
deltaf = 1./40.
flow = 10.
rhothr = 8.0

output = open('results_all.txt', 'w')

nsim = 100

data = np.loadtxt('all.txt')

idc = data[:, 0]
delt = data[:, 1]
m1 = data[:, 2]
m2 = data[:, 3]
e12 = data[:, 4]
iclus = data[:, 5]
mcl = data[:, 6]
z = data[:, 7]
bhf = data[:, 8]
nsf = data[:, 9]
km = data[:, 10]
scalr = data[:, 11]
scalt = data[:, 12]
scalv = data[:, 13]

zeta = [0.0001, 0.0002, 0.001, 0.005, 0.01, 0.02]
ecc = 0.0
dmini = 0.01

z_sf = np.arange(0., 10., 0.01)
max_zsf = max(z_sf)
min_zsf = min(z_sf)

sf = []
for i in range(len(z_sf)):
    sf.append(dt.star_formation(z_sf[i]))
arr_sf = np.array(sf)
max_sf = max(arr_sf)
min_sf = min(arr_sf)

massesc = [10000, 20000, 30000, 50000, 75000, 100000, 200000]
w_mc = []
sum_wmc = 0
for imc in range(len(massesc)):
    w_mc.append(1.0/massesc[imc] ** 2.0)
    sum_wmc = sum_wmc + 1.0/massesc[imc] ** 2.0
sum_wmc = np.float64(sum_wmc)
weight_mc = np.array(w_mc / sum_wmc)

for i in range(len(idc)):

    mb1 = m1[i]
    mb2 = m2[i]

    zeta_str = str(z[i])
    bhf_str = str(int(bhf[i]))
    nsf_str = str(int(nsf[i]))

    # sampling spins from pre-computed sse tracks

    ff = 'bse/sse_res_zeta_'+zeta_str+'_bhf_'+bhf_str+'_nsf_'+nsf_str
    fdata = np.loadtxt(ff)

    mbh = fdata[:, 1]
    abh = fdata[:, 3]

    if mb1 < 45.0:
        ibh1 = []
        dm = dmini
        while len(ibh1) <= 2:
            (ibh1,) = np.where((mbh < mb1 + dm) & (mbh > mb1 - dm))
            dm = dm * 1.5
        abh1max = max(abh[ibh1])
        abh1min = min(abh[ibh1])
        s1 = random.random() * (abh1max - abh1min) + abh1min
    else:
        s1 = 0.7
    #print(dm, mbh[ibh1], abh[ibh1])

    if mb2 < 45.0:
        ibh2 = []
        dm = dmini
        while len(ibh2) <= 2:
            (ibh2,) = np.where((mbh < mb2 + dm) & (mbh > mb2 - dm))
            dm = dm * 1.5
        abh2max = max(abh[ibh2])
        abh2min = min(abh[ibh2])
        s2 = random.random() * (abh2max - abh2min) + abh2min
    else:
        s2 = 0.7
    #print(dm, mbh[ibh2], abh[ibh2])
    #print(mb1, mb2, s1, s2)

    # weight cluster mass

    mclg = np.int64(mcl[i])
    prob_mc = weight_mc[np.where(massesc == mclg)][0]

    # nsim random realizations of spin orientations

    for j in range(0, nsim):

        print(i, j)

        # sample star formation

        yval = 10.0  # compute star formation redshift
        yfunc = 1.0
        while yval/yfunc > 1.0:
            xval = random.random() * (max_zsf - min_zsf) + min_zsf
            yfunc = dt.star_formation(xval)
            yval = random.random() * (max_sf - min_sf) + min_sf
        zzsf = xval

        tformation = Planck15.age(zzsf).value  # t formation of cluster
        tmerger = tformation + delt[i] / 1000.  # t merger (tformation + delay)
        if(tmerger < Planck15.age(0).value):
            zmerger = z_at_value(Planck15.age, tmerger * u.Gyr)

            # compute vkick chieff mfin chifin

            vkick, chieff = kc.kick(mb1, mb2, s1, s2, ecc)
            mfin = kc.mrem(mb1, mb2, s1, s2)  # compute final mass
            chifin = kc.spinrem(mb1, mb2, s1, s2)  # compute final spin

            # weight metallicity

            met_zzsf = dt.mean_metallicity(zzsf) * np.log(10.0)
            sigma_met_zzsf = 10 ** 0.5 * np.log(10.0)
            # distribution of metallicity at zzsf
            dist = lognorm(s=sigma_met_zzsf, loc=0.0, scale=np.exp(met_zzsf))
            w_zeta = []
            sum_wzeta = 0
            for iz in range(len(zeta)):
                w_zeta.append(dist.pdf(zeta[iz]))
                sum_wzeta = sum_wzeta + dist.pdf(zeta[iz])
            weight_zeta = np.array(w_zeta / sum_wzeta)
            prob_met = weight_zeta[np.where(zeta == z[i])][0]
            ran_met = random.random()

            # weight cosmology

            prob_cosmo = Planck15.differential_comoving_volume(zmerger).value
            prob_cosmo = 4.0 * np.pi * prob_cosmo / (1.0 + zmerger)

            # compute snr and detection probability

            snr = dt.snr(mb1, mb2, zmerger, appr, psd,
                         deltaf, flow)  # signal-to-noise
            prob_snr = dt.detec_prob_1d(rhothr/snr)  # detection probability
            ran_det = random.random()

            # overall weight

            prob_all = prob_mc * prob_met * prob_cosmo * prob_snr

            # write

            print(vkick, chieff, mfin, chifin, mb1, mb2, s1, s2, zzsf, tformation, zmerger, tmerger, met_zzsf, prob_met, prob_mc, prob_cosmo, snr, prob_snr,
                  prob_all, idc[i], delt[i], m1[i], m2[i], e12[i], iclus[i], mcl[i], z[i], bhf[i], nsf[i], km[i], scalr[i], scalt[i], scalv[i], file=output)

output.close()