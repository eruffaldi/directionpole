# Exercise Verifier for Sailing
#
# Emanuele Ruffaldi 2013
#
# Objective:
# Using Geodesic Distance tools: http://code.google.com/p/geopy/
# but using my own library for angles as exercise
# See also possible libraries:
# http://pint.readthedocs.org/en/latest/
# http://code.google.com/p/geopy/source/browse/trunk/geopy/distance.py?r=105
# http://code.google.com/p/geopy/source/browse/trunk/geopy/point.py?r=105
# http://www.movable-type.co.uk/scripts/latlong.html
#
# TODO: parse pure numbers
#
# p2-p1 -> bearing
# p1 + bearing + distance -> p2
#
# TODO: bearing
# http://stackoverflow.com/questions/17624310/geopy-calculating-gps-heading-bearing
#
# TODO ROUND LAT LON as follow
"""
if abs(latitude) > 90:
            latitude = ((latitude + 90) % 180) - 90

        longitude = float(longitude or 0.0)
        if abs(longitude) > 180:
            longitude = ((longitude + 180) % 360) - 180
"""
import math,re,datetime
#import geopy.distance,geopy.units
#from geopy.point import Point
#from geopy.compat import string_compare

# Point("39.3 N 76.4 W")
# p = Point("41.12345;-81.98765") 

def radians2compassdegrees(r):
	d = math.degrees(r)
	# radians are anti clockwise with zero at EAST and limited to 
	# compass is clockwise with zero at NORTH
	d = 180+d	
	if d < 0:
		d = 360+d
	while d >= 360:
		d = 360-d
	return d

class AngCoord:
	Other = 0
	Lat = 1
	Long = 2
	
	# ns can be: "E" "W" "N" "S", or empty
	def __init__(self,deg,primes,ns):
		self.value = deg
		self.family = AngCoord.Other
		if primes is not None:
			self.value += primes / 60.0
		if type(ns) is str:
			if ns == "N":
				self.family = AngCoord.Lat
			elif ns == "S":
				self.family = AngCoord.Lat
				self.value = -self.value
			elif ns == "W":
				self.family = AngCoord.Long
				self.value = -self.value
			elif ns == "E":
				self.family = AngCoord.Long
			elif ns == "":
				self.family = AngCoord.Other
			else:
				raise Exception("Unsupported ns suffix: " + ns)	
		elif type(self.family) is int:
			self.family = ns
		elif ns is not None:
			raise Exception("Unsupported family value" + str(ns))
	def clone(self):
		return AngCoord(self.value,None,self.family)	
	def __div__(self,x):
		return AngCoord(self.value/x,None,self.family)	
	@staticmethod
	def parselow(s):
		ins =s
		degs = 0
		s = s.strip()
		k = s.find("\xB0")			
		if k < 0:
			k = s.find("d")			
		if k >= 0:
			degs = float(s[0:k])
			s = s[k+1:].strip()
		if s == "":
			return (degs,0,"")
		s = s.replace("''","\"")
		g = re.match("([NSEW])?",s.strip())
		if g:
			ns = g.group(1)
			tmp  = degs;
			degs = math.floor(tmp)
			primes = (tmp-degs)*60.0
			return (degs,primes,ns)
		else:
			g = re.match("([0-9.,]+)'([.,][0-9]+)?( [0-9.]+\")? ([NSEW])?",s)
			if g:
				# returns (degree,primes,direction letter)
				hasFracPrimes = g.group(1).find(".") >= 0 or g.group(1).find(",") >= 0 
				if g.group(2): # seconds
					if hasFracPrimes:
						raise Exception("Floating double specified in primes")
					hasFracPrimes = True
					primes = int(g.group(1)) + float(g.group(2)[1:])/60.0
				else: # just primes
					primes = float(g.group(1).replace(",","."))
				if g.group(3):
					if hasFracPrimes:
						raise Exception("Floating double specified in primes. Cannot specify seconds")
					primes = primes + float(g.group(3)[0:-1])/60.0
				if g.group(4):
					ns = g.group(4)
				else:
					ns = ""
				return (degs,primes,ns)
			else:
				raise Exception("AngCoord Format mismatch on: <%s> of <%s>" % (s,ins))
				return None
	@staticmethod
	def parse(s):
		r = AngCoord.parselow(s)
		return AngCoord(r[0],r[1],r[2])
	def __add__(self,other):
		if other.family != self.family:
			raise Exception("Family mismatch in add")
		r = self.clone()
		r.value = self.value + other.value
		return r
	def __sub__(self,other):
		if other.family != self.family:
			raise Exception("Family mismatch in sub")
		r = self.clone()
		r.value = self.value - other.value
		return r
	def __mul__(self,other):
		assert type(other) == int or type(other) == float
		r = self.clone()
		r.value = self.value * other
		return r
	def geopystr(self):
		isneg = self.value < 0
		av = math.fabs(self.value)
		if self.family == AngCoord.Lat:
			return "%f %s" % (av,isneg and "S" or "N")
		else:
			return "%f %s"  % (av,isneg and "W" or "E")
	def __str__(self):
		isneg = self.value < 0
		av = math.fabs(self.value)
		ip = int(math.floor(av))
		fp = (av-ip)*60
		if self.family == AngCoord.Lat:
			return "%03dd %02.2f' %s" % (ip,fp,isneg and "S" or "N")
		elif self.family == AngCoord.Long:
			return "%02dd %02.2f' %s" % (ip,fp,isneg and "W" or "E")
		else:
			return "%s%02dd %02.2f'" % (isneg and "-" or "+",ip,fp)

class LatCoord(AngCoord):
	def __init__(self,deg,primes=None,ns=None):
		AngCoord.__init__(self,deg,primes,ns)
		if ns is None:
			self.family = AngCoord.Lat
		elif self.family != AngCoord.Lat:
			raise Exception("Mismatch coordinate family")
	@staticmethod
	def parse(s):
		r = AngCoord.parselow(s)
		return LatCoord(r[0],r[1],r[2])
	def clone(self):
		return LatCoord(self.value,None,AngCoord.Lat)


class LongCoord(AngCoord):
	def __init__(self,deg,primes=None,ns=None):
		AngCoord.__init__(self,deg,primes,ns)
		if ns is None:
			self.family = AngCoord.Lat
		elif self.family != AngCoord.Long:
			raise Exception("Mismatch coordinate family")
	@staticmethod
	def parse(s):
		r = AngCoord.parselow(s)
		return LongCoord(r[0],r[1],r[2])
	def clone(self):
		return LongCoord(self.value,None,AngCoord.Long)

class Declination:
	def __init__(self,value,year,increment):
		assert isinstance(value,AngCoord)
		assert isinstance(increment,AngCoord)
		self.value = value
		self.year = year
		self.increment = increment
	def actualize(self,year):		
		return Declination(self.value + self.increment * (year-self.year),year,self.increment)
	def __str__(self):
		return "Declination(%s,%d,%s)" % (self.value,self.year,self.increment)

class GeoVector:
	"""Pure Geographical Vector where direction is in Degrees counting clockwise, and distance is in nautical miles"""
	def __init__(self,directionDeg,distanceMi=1):
		self.directionDeg = directionDeg
		self.distanceMi = distanceMi
	def correct(self,declination):
		assert isinstance(declination,Declination)
		self.directionDeg += declination.value
	def convert(self,declination):
		assert isinstance(declination,Declination)
		self.directionDeg -= declination.value
	def eta(self,starttime,etaKnots):
		assert isinstance(starttime,datetime.datetime)
		return starttime + datetime.timedelta(hours=self.distanceMi/etaKnots)
	def intersect(self,v):
		assert isinstance(v,GeoVector)
		raise "NotImpl"
	def distanceToPoint(self,pt):
		assert isinstance(pt,GeoPoint)
		raise "NotImpl"
	def normalize(self):
		return GeoVector(self.directionDeg,1)
	def __add__(self,other):
		raise "NotImpl"
		pass
	def traverso(self,rightSide):
		return
		raise "NotImpl"
	def rotate(self,deg):
		raise "NotImpl"
		return
	def anglebetween(self):
		raise "NotImpl"
		return
	def __str__(self):
		return "Vector(%s,%f mi)" % (str(self.directionDeg),self.distanceMi)

class GeoPoint:
	"""Pure Geographical Point expressed by latitude and longitude"""
	def __init__(self,lat,lon):
		assert isinstance(lat,LatCoord) or type(lat) is str or type(lat) is float
		assert isinstance(lon,LongCoord) or type(lon) is str or type(lon) is float
		if type(lat) is str:
			lat = LatCoord.parse(lat)
		elif type(lat) is float:
			lat = LatCoord(lat,None,None)
		if type(lon) is str:
			lon = LongCoord.parse(lon)
		elif type(lon) is float:
			lon = LongCoord(lon,None,None)
		self.lat = lat
		self.lon = lon
	def cartesian(self):
		la = self.lat.value*math.pi/180
		lo = self.lon.value*math.pi/180
		cla = math.cos(la)
		clo = math.cos(lo)
		sla = math.sin(la)
		slo = math.sin(lo)
		return (cla*clo,cla*slo,sla)
	@staticmethod
	def fromcartesian(ca):
		x,y,z = ca
		lon = math.atan2(y,x)
		hyp = math.sqrt(x*x+y*y)
		lat = math.atan2(z,hyp)
		return GeoPoint(lat*180/math.pi,lon*180/math.pi)
	@staticmethod
	def fromtuple(tup):
		return GeoPoint(tup[0],tup[1])		
	def simplemidpoint(self,other):
		return GeoPoint((self.lat.value+other.lat.value)*0.5,(self.lon.value+other.lon.value)*0.5)
		#return GeoPoint.centroid([self,other])
	def astuple(self):
		return (self.lat.value,self.lon.value)
	def asgeopy(self):
		return self.lat.geopystr() + " " + self.lon.geopystr()
	def distance(self,other):
		if isinstance(other,GeoPoint):
			raise "Not Implemented"
		elif isinstance(other,GeoVector) or isinstance(other,GeoPointVector):
			return other.distance(self)
	@staticmethod
	def weightedcentroid(pts,weights):
		assert len(pts) == len(weights)
		raise "Not Implemented"
	@staticmethod
	def centroid(pts):
		return GeoPoint.weightedcentroid(pts,[1 for i in range(0,len(pts))])
	def __add__(self,other):
		assert isinstance(other,GeoVector)
		import geodistance
		d = geodistance.xdistance()
		return d.destination(self,other.directionDeg,other.distanceMi)
	def coordDelta(self,otherorigin):
		assert isinstance(otherorigin,GeoPoint)
		return GeoPoint(self.lat-otherorigin.lat,self.lon-otherorigin.lon)
	def coordDirection(self):
		"""computes the angle between latitude and longitude vector as if being 2D vector. Use in planar"""
		raise "Not Implemented"
	def minModulus(self):
		"""modulus in minutes. Used for planar navigation"""
		return math.sqrt((self.lat.value*60.0)**2+(self.lon.value*60.0)**2)
	def __sub__(self,other):
		# TODO: bearing is ...?
		import geodistance
		d = geodistance.xdistance().measure(self,other)		
		return GeoVector(d["startbearing"],d["distanceKm"]/1.852)
	def __str__(self):	
		return "Lat: %s  Lon: %s" % (str(self.lat),str(self.lon))
	def __repr__(self):
		return self.__str__()
	@staticmethod
	def fromgeopy(gp):
		return GeoPoint(LatCoord(gp.latitude,None,None),LongCoord(gp.longitude,None,None))


class GeoPointVector:
	"""Pure Geographical Origin with Direction, used for Surveys"""
	def __init__(self,p,v):
		"""Construct by a point and a vector"""
		assert isinstance(p,GeoPoint) and isinstance(v,GeoVector)
		self.point = p
		self.dir = v
	def pointAtMiles(self,distanceMi):
		return self.point + GeoVector(self.dir.directionDeg,distanceMi)
	def distance(self,pt):
		# specialized because it looks at the correct side
		if isinstance(pt,GeoPointVector):
			pass
		elif isinstance(pt,GeoPoint):
			pass
		else:
			raise ValueError("op not supported")
		raise "Not Implemented"
	def translate(self,v):
		assert isinstance(v,GeoVector)
		# translate point component by amount
		raise "Not Implemented"

	def intersect(self,ptv):
		raise "Not Implemented"
		if isinstance(ptv,GeoVector):
			pass
		elif isinstance(ptv,GeoPointVector):
			pass

def triangleinfo(pa,pb,pc):
	#compute the lengths of the opposite edges
	a = pb.distance(pc)
	b = pc.distance(pa)
	c = pa.distance(pb)
	P = a+b+c
	return dict(incenter=weightedcentroid([pa,pb,pc],(a/P,b/P,c/P)),area=1/4*math.sqrt(P*(a-b+c)*(b-c+a)*(c-a+b)),perimeter=P,sides=(a,b,c))


# given a direction (estimated) and a distance in Miles between two surveys of the same or different points
#
# note that all directions are ABSOLUTE
#
# Returns the position and the area of the triangle
def locateTwoSurveys(direction,survey1,survey2,distanceMi):
	assert isinstance(direction,GeoVector) and isinstance(survey1,GeoPointVector) and isinstance(survey2,GeoPointVector)
	if distanceMi == 0:
		return (survey1.intersect(survey2),0)
	else:
		# correct the second based on the traveled distance along the direction
		# then intersect all the directions and obtain a triangle
		survey1moved = survey1.translate(GeoVector(direction.directionDeg,distanceMi))
		p1 = survey1moved.intersect(survey2)

		# Following was old code for triangle intersection
		#p2 = survey2.intersect(survey1moved)
		#p3 = survey2.intersect(direction)
		#(d["incenter"],d["area"])
		# assume in the incenter of triangle
		#d = triangleinfo(p1,p2,p3)
		return p1 

if __name__ == "__main__":
	
	# Varrone pag 20 #4
	cases = [("4d 42' W",2005,"10' W",2010),("5d 20' E",2000,"10' E",2010),("1d 15' W",2002,"9' E",2010),("0d 09' E",2002,"8' W",2010)]
	for d,y,i,ty in cases:
		print "Case ",d,y,i,ty
		print "\tparsed:",AngCoord.parselow(d),"value:",AngCoord.parse(d).value,"formatted:",AngCoord.parse(d),"reparsed:",AngCoord.parse(str(AngCoord.parse(d)))
		print "\tincrement:",AngCoord.parse(i)
		print "\tresult:",Declination(AngCoord.parse(d),y,AngCoord.parse(i)).actualize(ty)
	
	# Varrone pag 34 #22
	cases2 = [("3d 10' W",1988,"9' W",2004),("0d 46' E",1999,"10' W",2006),("4d 50' W",2000,"6' E",2003),("9d 30' E",1990,"10' W",2000)]
	for d,y,i,ty in cases2:
		print "Case ",d,"@",y,"by",i,"year ->",ty
		print "\tparsed:",AngCoord.parselow(d),"value:",AngCoord.parse(d).value,"formatted:",AngCoord.parse(d),"reparsed:",AngCoord.parse(str(AngCoord.parse(d)))
		print "\tincrement:",AngCoord.parse(i)
		print "\tresult:",Declination(AngCoord.parse(d),y,AngCoord.parse(i)).actualize(ty)

 	p1 = GeoPoint("42d 33' N","010d 24.5' E")
	p2 = GeoPoint("42d 40' N","010d 38' E")
	t0 = datetime.datetime.now()
	print "Now",t0
	print "From",p1,"To",p2
	d = p2-p1
	print "\tVector:",d
	print "\tEta at 10knots:",d.eta(t0,10)
