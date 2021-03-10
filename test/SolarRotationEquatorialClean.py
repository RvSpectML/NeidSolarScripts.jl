# Python code by Shubham Kanodia March 2021 (with input from Jason Wright and Eric Ford)

from astropy.coordinates import EarthLocation
from astropy.coordinates import get_body_barycentric_posvel, get_body_barycentric
from astropy.time import Time
import astropy.constants as ac

import os
import numpy as np
import pandas as pd

import barycorrpy
from barycorrpy import PINT_erfautils as PINT
from barycorrpy.utils import CalculatePositionVector
from barycorrpy.PhysicalConstants import *

from scipy.spatial.transform import Rotation as R

# For full precision, do not use obsname='KPNO', but use actual WIYN coords
obsname = 'KPNO'
ephemeris = 'de430'

obsname = None
longi = -17.8893
lat = 28.7541
alt = 2370


if obsname:
	loc = EarthLocation.of_site(obsname)
	lat = loc.lat.value
	longi = loc.lon.value
	alt = loc.height.value
else:
	loc = EarthLocation.from_geodetic(longi, lat, height=alt)


HARPSN_df = pd.read_csv("../data/Sun_harpsn_qualflag.csv")

jd = np.array(HARPSN_df['JD']) + 2400000
F_obs = np.array(HARPSN_df['FWHMobs'])
F_sid = np.array(HARPSN_df['FWHMsid'])
HARPS_Delta = F_obs**2 - F_sid**2

JDUTC_master = Time(jd[::50] , format='jd', scale='utc')


def CalculateFWHMDifference_SolarRotation_Equatorial(loc, JDUTC):
	"""
	Calculate the difference between the Observed Solar FWHM and Sidereal Solar FWHM
	Based on Colier Cameron et al. (2019)

	INPUTS:
		loc: Astropy Earth Location object. https://docs.astropy.org/en/stable/api/astropy.coordinates.EarthLocation.html
		JDUTC: Astropy Time Object.

	OUTPUT:
		Delta:  F_obs**2 - F_sid**2 [(km/s)^2]
	"""
	# Convert times to obtain TDB and TT
	JDTDB = JDUTC.tdb
	JDTT = JDUTC.tt

	################################
	######EARTH EPHEMERIS ############
	################################
	##### NUTATION, PRECESSION, ETC. #####

	# Observatory position wrt Geocenter
	r_pint, v_pint = PINT.gcrs_posvel_from_itrf(loc, JDUTC, JDTT)

	r_eci = r_pint[0]  # [m]
	v_eci = v_pint[0]  # [m/s]

	##### EPHEMERIDES #####

	earth_geo = get_body_barycentric_posvel('earth', JDTDB, ephemeris=ephemeris) # [km]
	r_geo = np.reshape(earth_geo[0].xyz.value*1000., 3) # [m]
	v_geo = np.reshape(earth_geo[1].xyz.value*1000./86400., 3)  # [m/s]

	PosVector_EarthSSB = r_eci + r_geo # [m]

	# Relativistic Addition of Velocities
	VelVector_EarthSSB = (v_eci+v_geo) / (1.+ np.sum(v_eci*v_geo)/c**2) # [m/s]

	################################
	######SOLAR EPHEMERIS ############
	################################

	solar_ephem = get_body_barycentric_posvel('sun', JDTDB, ephemeris=ephemeris)

	PosVector_SolSSB = np.reshape(solar_ephem[0].xyz.value*1000., 3) #[m]
	VelVector_SolSSB = np.reshape(solar_ephem[1].xyz.value*1000./86400., 3)  # [m/s]


	################################
	####EQUATORIAL COORD VECTORS ####
	################################

	PosVector_EarthSol, PosMag_EarthSol, PosHat_EarthSol = CalculatePositionVector(r1=PosVector_EarthSSB, r2=PosVector_SolSSB)
	VelVector_EarthSol = (VelVector_EarthSSB - VelVector_SolSSB) / (1. + np.sum(VelVector_SolSSB*VelVector_EarthSSB)/c**2)

	# omega = (x.vy - y.vx)/|(x,y)|^2
	# Omega = (PosVector_EarthSol[0]*VelVector_EarthSol[1] - PosVector_EarthSol[1]*VelVector_EarthSol[0]) / (PosMag_EarthSol**2)
	OmegaVector = np.cross(PosVector_EarthSol, VelVector_EarthSol) / (PosMag_EarthSol**2)


	################################
	################################


	# Rotation Axis of the Sun with respect to Equatorial System
	# https://www2.mps.mpg.de/homes/fraenz/systems/systems3art.pdf

	#Eqn 13
	alpha = 286.13 # Rotate along the z axis
	dec = 63.87 # Rotate along the x axis

	Theta = alpha * np.pi/180
	Phi = dec * np.pi/180

	EclipticEpsilon = 23.44

	SolarInclination = 7.25
	SolarLongitude = 75.76

	# Theta = SolarLongitude * 180/np.pi
	# Phi = SolarInclination * 180/np.pi

	# Need to perform Extrinsic Euler Rotations

	# Transform to solar rotation axis
	REquatorial = R.from_euler("zyx", [alpha, 0, dec], degrees=True).as_matrix()

	REquatorial = np.array([np.cos(Theta)*np.cos(Phi), np.sin(Theta)*np.cos(Phi), np.sin(Phi)])


	"""
	################################
	####ROTATED COORD VECTORS ####
	################################

	RotatedPositionVector = np.matmul(Rzyx, PosVector_EarthSol)
	RotatedPositionVector = PosVector_EarthSol
	RotatedPositionHat = RotatedPositionVector / np.linalg.norm(RotatedPositionVector)

	OmegaEarthVectorRotated = np.matmul(Rzyx, OmegaVector)

	################################
	################################
	"""

	OmegaSolVector = np.array([0, 0,  2.972 *1e-6])
	OmegaSolHat = np.array([0,0,1])

	# Rotated Solar rotation vector to ecliptic plane

	# OmegaSolVector = np.matmul(REquatorial, OmegaSolVector)
	OmegaSolVector = REquatorial * np.linalg.norm(OmegaSolVector)

	OmegaSolHat = OmegaSolVector / np.linalg.norm(OmegaSolVector)

	sini = np.sqrt(1 - np.matmul(OmegaSolHat, PosHat_EarthSol)**2)
	Gamma = 1.04
	DeltaOmega = OmegaSolVector - OmegaVector


	Delta = ((Gamma* ac.R_sun.to(u.km).value)**2) * (np.matmul(DeltaOmega, DeltaOmega)*sini*sini - np.matmul(OmegaSolVector, OmegaSolVector))

	return Delta


Delta = np.array([CalculateFWHMDifference_SolarRotation(loc, JDUTC) for JDUTC in JDUTC_master])

#########
plt.plot(JDUTC_master.jd, HARPS_Delta[::50]-Delta)
plt.ylabel("HARPS - BCPy (km/s)^2")
plt.show(block=False)
