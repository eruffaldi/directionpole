#See: http://geographiclib.sourceforge.net/
from geocoord import *

ELLIPSOIDS = {
	# model           major (km)   minor (km)     flattening
	'WGS-84':        (6378.137,    6356.7523142,  1 / 298.257223563),
	'GRS-80':        (6378.137,    6356.7523141,  1 / 298.257222101),
	'Airy (1830)':   (6377.563396, 6356.256909,   1 / 299.3249646),
	'Intl 1924':     (6378.388,    6356.911946,   1 / 297.0),
	'Clarke (1880)': (6378.249145, 6356.51486955, 1 / 293.465),
	'GRS-67':        (6378.1600,   6356.774719,   1 / 298.25)
}

import numpy
from math import *
from numpy import *

# taken from here:
# http://gis.stackexchange.com/questions/66/trilateration-using-3-latitude-and-longitude-points-and-3-distances
def laterization(LatA,LonA,DistA,LatB,LonB,DistB,LatC,LonC,DistC):
	#assuming elevation = 0
	earthR = 6371

	#using authalic sphere
	#if using an ellipsoid this step is slightly different
	#Convert geodetic Lat/Long to ECEF xyz
	#   1. Convert Lat/Long to radians
	#   2. Convert Lat/Long(radians) to ECEF
	xA = earthR *(math.cos(math.radians(LatA)) * math.cos(math.radians(LonA)))
	yA = earthR *(math.cos(math.radians(LatA)) * math.sin(math.radians(LonA)))
	zA = earthR *(math.sin(math.radians(LatA)))

	xB = earthR *(math.cos(math.radians(LatB)) * math.cos(math.radians(LonB)))
	yB = earthR *(math.cos(math.radians(LatB)) * math.sin(math.radians(LonB)))
	zB = earthR *(math.sin(math.radians(LatB)))

	xC = earthR *(math.cos(math.radians(LatC)) * math.cos(math.radians(LonC)))
	yC = earthR *(math.cos(math.radians(LatC)) * math.sin(math.radians(LonC)))
	zC = earthR *(math.sin(math.radians(LatC)))

	P1 = array([xA, yA, zA])
	P2 = array([xB, yB, zB])
	P3 = array([xC, yC, zC])

	#from wikipedia
	#transform to get circle 1 at origin
	#transform to get circle 2 on x axis
	ex = (P2 - P1)/(numpy.linalg.norm(P2 - P1))
	i = dot(ex, P3 - P1)
	ey = (P3 - P1 - i*ex)/(numpy.linalg.norm(P3 - P1 - i*ex))
	ez = numpy.cross(ex,ey)
	d = numpy.linalg.norm(P2 - P1)
	j = dot(ey, P3 - P1)

	#from wikipedia
	#plug and chug using above values
	x = (pow(DistA,2) - pow(DistB,2) + pow(d,2))/(2*d)
	y = ((pow(DistA,2) - pow(DistC,2) + pow(i,2) + pow(j,2))/(2*j)) - ((i/j)*x)

	# only one case shown here
	a = pow(DistA,2) - pow(x,2) - pow(y,2)
	if a < 0:
		return None
	else:
		z = sqrt(a)

	#triPt is an array with ECEF x,y,z of trilateration point
	triPt = P1 + x*ex + y*ey + z*ez

	#convert back to lat/long from ECEF
	#convert to degrees
	lat = math.degrees(math.asin(triPt[2] / earthR))
	lon = math.degrees(math.atan2(triPt[1],triPt[0]))
	return (lat,lon)




#from: https://github.com/geopy/geopy
class xdistance: #geopy.distance.distance):
	ellipsoid_key = None
	ELLIPSOID = None
	def __init__(self,**kwargs):
		self.set_ellipsoid(kwargs.pop('ellipsoid', 'WGS-84'))
		major, minor, f = self.ELLIPSOID # pylint: disable=W0612
	def set_ellipsoid(self, ellipsoid):
			"""
			Change the ellipsoid used in the calculation.
			"""
			if not isinstance(ellipsoid, (list, tuple)):
				try:
					self.ELLIPSOID = ELLIPSOIDS[ellipsoid]
					self.ellipsoid_key = ellipsoid
				except KeyError:
					raise Exception(
						"Invalid ellipsoid. See geopy.distance.ELIPSOIDS"
					)
			else:
				self.ELLIPSOID = ellipsoid
				self.ellipsoid_key = None
			return

	def measure(self, a, b):
		assert isinstance(a,GeoPoint)
		assert isinstance(b,GeoPoint)
		from math import radians
		#a, b = Point(a), Point(b)
		lat1, lng1 = radians(a.lat.value), radians(a.lon.value)
		lat2, lng2 = radians(b.lat.value), radians(b.lon.value)

		if isinstance(self.ELLIPSOID, str):
			major, minor, f = ELLIPSOIDS[self.ELLIPSOID] #geopy.distance.ELLIPSOIDS[self.ELLIPSOID]
		else:
		    major, minor, f = self.ELLIPSOID
		from math import atan,tan,sin,cos,pi,sqrt,atan2
		delta_lng = lng2 - lng1

	

		reduced_lat1 = atan((1 - f) * tan(lat1))
		reduced_lat2 = atan((1 - f) * tan(lat2))

		sin_reduced1, cos_reduced1 = sin(reduced_lat1), cos(reduced_lat1)
		sin_reduced2, cos_reduced2 = sin(reduced_lat2), cos(reduced_lat2)

		lambda_lng = delta_lng
		lambda_prime = 2 * pi

		iter_limit = 20

		while abs(lambda_lng - lambda_prime) > 10e-12 and iter_limit > 0:
			sin_lambda_lng, cos_lambda_lng = sin(lambda_lng), cos(lambda_lng)

			sin_sigma = sqrt(
				(cos_reduced2 * sin_lambda_lng) ** 2 +
				(cos_reduced1 * sin_reduced2 -
				 sin_reduced1 * cos_reduced2 * cos_lambda_lng) ** 2
			)

			if sin_sigma == 0:
				return dict(distanceKm=0,startbearing=0,endbearing=0) # Coincident points

			cos_sigma = (
				sin_reduced1 * sin_reduced2 +
				cos_reduced1 * cos_reduced2 * cos_lambda_lng
			)

			sigma = atan2(sin_sigma, cos_sigma)

			sin_alpha = (
				cos_reduced1 * cos_reduced2 * sin_lambda_lng / sin_sigma
			)
			cos_sq_alpha = 1 - sin_alpha ** 2

			if cos_sq_alpha != 0:
				cos2_sigma_m = cos_sigma - 2 * (
					sin_reduced1 * sin_reduced2 / cos_sq_alpha
				)
			else:
				cos2_sigma_m = 0.0 # Equatorial line

			C = f / 16. * cos_sq_alpha * (4 + f * (4 - 3 * cos_sq_alpha))

			lambda_prime = lambda_lng
			lambda_lng = (
				delta_lng + (1 - C) * f * sin_alpha * (
					sigma + C * sin_sigma * (
						cos2_sigma_m + C * cos_sigma * (
							-1 + 2 * cos2_sigma_m ** 2
						)
					)
				)
			)
			iter_limit -= 1

		if iter_limit == 0:
			raise ValueError("Vincenty formula failed to converge!")

		u_sq = cos_sq_alpha * (major ** 2 - minor ** 2) / minor ** 2

		A = 1 + u_sq / 16384. * (
			4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq))
		)

		B = u_sq / 1024. * (256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq)))

		delta_sigma = (
			B * sin_sigma * (
				cos2_sigma_m + B / 4. * (
					cos_sigma * (
						-1 + 2 * cos2_sigma_m ** 2
					) - B / 6. * cos2_sigma_m * (
						-3 + 4 * sin_sigma ** 2
					) * (
						-3 + 4 * cos2_sigma_m ** 2
					)
				)
			)
		)
		#http://en.wikipedia.org/wiki/Vincenty%27s_formulae
		#U1 = reduced_lat1
		#U2 = reduced_lat2
		cl = cos_lambda_lng 
		sl = sin_lambda_lng
		c1 = cos_reduced1
		c2 = cos_reduced2
		s1 = sin_reduced1
		s2 = sin_reduced2
		a1 = math.atan2(c2* sl, c1*s2 - s1*c2*cl)
		a2 = math.atan2(c1* sl, -s1*c2 + c1*s2*cl)
		startbearing = radians2compassdegrees(a1)
		endbearing = radians2compassdegrees(a2)

		s = minor * A * (sigma - delta_sigma)
		return dict(distanceKm=s,startbearing=startbearing,endbearing=endbearing)        
	def destination(self, point, bearingDeg, distanceMi): # pylint: disable=W0621
		"""
		TODO docs.
		"""
		assert isinstance(a,GeoPoint)
		lat1 = math.radians(point.lat.value)
		lng1 = math.radians(point.lon.value)
		bearing = math.radians(bearing)

		distance = distanceMi*1.852

		ellipsoid = self.ELLIPSOID
		if isinstance(ellipsoid, string_compare):
			ellipsoid = ELLIPSOIDS[ellipsoid]

		major, minor, f = ellipsoid

		tan_reduced1 = (1 - f) * tan(lat1)
		cos_reduced1 = 1 / sqrt(1 + tan_reduced1 ** 2)
		sin_reduced1 = tan_reduced1 * cos_reduced1
		sin_bearing, cos_bearing = sin(bearing), cos(bearing)
		sigma1 = atan2(tan_reduced1, cos_bearing)
		sin_alpha = cos_reduced1 * sin_bearing
		cos_sq_alpha = 1 - sin_alpha ** 2
		u_sq = cos_sq_alpha * (major ** 2 - minor ** 2) / minor ** 2

		A = 1 + u_sq / 16384. * (
			4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq))
		)
		B = u_sq / 1024. * (256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq)))

		sigma = distance / (minor * A)
		sigma_prime = 2 * pi

		while abs(sigma - sigma_prime) > 10e-12:
			cos2_sigma_m = cos(2 * sigma1 + sigma)
			sin_sigma, cos_sigma = sin(sigma), cos(sigma)
			delta_sigma = B * sin_sigma * (
				cos2_sigma_m + B / 4. * (
					cos_sigma * (
						-1 + 2 * cos2_sigma_m
					) - B / 6. * cos2_sigma_m * (
						-3 + 4 * sin_sigma ** 2) * (
						-3 + 4 * cos2_sigma_m ** 2
					)
				)
			)
			sigma_prime = sigma
			sigma = distance / (minor * A) + delta_sigma

		sin_sigma, cos_sigma = sin(sigma), cos(sigma)

		lat2 = atan2(
			sin_reduced1 * cos_sigma + cos_reduced1 * sin_sigma * cos_bearing,
			(1 - f) * sqrt(
				sin_alpha ** 2 + (
					sin_reduced1 * sin_sigma -
					cos_reduced1 * cos_sigma * cos_bearing
				) ** 2
			)
		)

		lambda_lng = atan2(
			sin_sigma * sin_bearing,
			cos_reduced1 * cos_sigma - sin_reduced1 * sin_sigma * cos_bearing
		)

		C = f / 16. * cos_sq_alpha * (4 + f * (4 - 3 * cos_sq_alpha))

		delta_lng = (
			lambda_lng - (1 - C) * f * sin_alpha * (
				sigma + C * sin_sigma * (
					cos2_sigma_m + C * cos_sigma * (
						-1 + 2 * cos2_sigma_m ** 2
					)
				)
			)
		)

		lng2 = lng1 + delta_lng

		return GeoPoint(math.degrees(lat2),math.degrees(lng2))
