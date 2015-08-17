# Where is my Direction Pole?
# Emanuele Ruffaldi 2014
#
# Requirements: pykml scipy

from geodistance import *
import itertools
from scipy.optimize import minimize,leastsq
import numpy as np
import sys

# all from Google Search: name latitude longitude
locations = dict(
	photo="42d 50' 57.822'' N, 11d 41' 3.846'' E",
	locpole="42.775842d N, 11.695500d E",
	locpole2="42.773775d N, 11.697626d E",
	village="42.7667d N, 11.7000d E",
	BuenosAires="34.6033d S, 58.3817d W",
	London="51.5072d N, 0.1275d W",
	Cophenhagen="55.6761d N, 12.5683d E",
	SanPaolo="23.5500d S, 46.6333d W",
	Vienna="48.2000d N, 16.3667d E",
	Marrakech ="31.6300d N, 8.0089d W",
	Moscow="55.7500d N, 37.6167d E",
	Sofia="42.7000d N, 23.3333d E",
	Chisenau="47.0000d N, 28.9167d E",
	Varsavia="52.2333d N, 21.0167d E",
	Minsk="53.9000d N, 27.5667d E",
	Kiev="50.4500d N, 30.5233d E",
	Bucharest="44.4325d N, 26.1039d E",
	Skopje="42.0000d N, 21.4333d E",
	NewYork="40.7127d N, 74.0059d W",
	Paris="48.8567d N, 2.3508d E",
	Bern="46.9500d N, 7.4500d E",
	Berlin="52.5167d N, 13.3833d E")

# Photo 1 cities
distances1 = dict(
	BuenosAires=11175,
	Varsavia=1262,
	Minsk=1700,
	Kiev=1667,
	Bucharest=1174,
	Skopje=805,
	NewYork=6794
	)

# Photo 2 cities
distances2 = dict(
	London=1399,
	Cophenhagen=1440,
	SanPaolo=9499,
	Vienna=706,
	Marrakech=2135,
	Moscow=2350,
	Chisenau=1429,
	Sofia=950
	)

# PoleB Photo 1 cities
distancesB1 = dict(
	Paris=991,
	Bern=573,
	Berlin=1093,
	BuenosAires=11175,
	Varsavia=1262,
	Minsk=1700,
	Kiev=1667,
	Skopje=950
	)

# PoleB Photo 2 cities
distancesB2 = dict(
	London=1399,
	Cophenhagen=1440,
	SanPaolo=9499,
	Vienna=706,
	Marrakech=2135,
	Moscow=2350,
	Chisenau=1429
	)

# finds best using trilaterization
def findbesttrilat(target,ltarget):
	allp = []
	mind = None
	mindset = None
	minp = None
	failedcount = 0
	items = target.items()
	print "*",items[0]
	for iset in itertools.combinations(range(0,len(target)),3):

		dss = [items[x][1][0] for x in iset]
		pts = [items[x][1][1] for x in iset]
		curset = [items[x][0] for x in iset]
		#print dss,pts,curset

		try:
			r = laterization(pts[0].lat.value,pts[0].lon.value,dss[0],pts[1].lat.value,pts[1].lon.value,dss[1],pts[2].lat.value,pts[2].lon.value,dss[2])
		except:
			print iset,curset,"failed"
			failedcount += 1
			continue
		if r == None or math.isnan(r[0]) or math.isnan(r[1]):
			print iset,curset,"failed"
			failedcount += 1
		else:
			p = GeoPoint(r[0],r[1])
			dd = xd.measure(p,ltarget)["distanceKm"]
			allp.append((curset,p,dd))
			#print iset,dd
			if mind is None or mind > dd:
				mind = dd
				mindset = curset
				minp = p
	# now trilateration of multilateration given subsets of 2 elements of each distance pole
	print "best",mind,mindset
	print "failed",failedcount
	return allp,(minp,mindset,mind)

#computes center of gravity of points
def cog(points):
	pp = None
	n = 0
	for p in points:
		if pp is None:
			pp = p.cartesian()
		else:
			x = p.cartesian()
			pp = (pp[0]+x[0],pp[1]+x[1],pp[2]+x[2])
			n += 1
	pp = (pp[0]/n,pp[1]/n,pp[2]/n)
	return GeoPoint.fromcartesian(pp)



kmloutput = None


def kmlpush(pt,label,stylename="sn_shaded_dot"):
	global kmloutput,outkml
	from pykml.factory import KML_ElementMaker as KML
	if kmloutput is None:
		kmloutput = KML.kml(
		    KML.Document(
		        KML.Name("Sun Position"),
		        KML.Style(
			      KML.IconStyle(
			        KML.scale(1.0),
			        KML.Icon(
			          KML.href("http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png"),
			        ),
			        id="mystyle"
			      ),
			      id="pushpin"
			    ),
		        KML.Style(
			      KML.IconStyle(
			        KML.scale(1.0),
			        KML.Icon(
			          KML.href("http://maps.google.com/mapfiles/kml/pushpin/red-pushpin.png"),
			        ),
			        id="redmystyle"
			      ),
			      id="redpushpin"
			    ),
		        KML.Style(
		            KML.IconStyle(
		                KML.scale(1.2),
		                KML.Icon(
		                    KML.href("http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png")
		                ),
		            ),
		            id="sn_shaded_dot",
		        )
		    )
		)
	if pt is None:
		from lxml import etree
		open(outkml,"wb").write(etree.tostring(kmloutput, pretty_print=True))
		return
	pt = KML.Placemark(
	    KML.name(label),
	    KML.styleUrl('#{0}'.format(stylename)),
	    KML.Point(
	        KML.extrude(True),
	        KML.altitudeMode('relativeToGround'),
	        KML.coordinates("{0},{1}".format(pt.lon.value,pt.lat.value)),
	    ),
	)
	kmloutput.Document.append(pt)



#compute point that is the average of the trilat
def faveraginglat(c,pts):
	pp = GeoPoint(LatCoord(c[0],None,None),LongCoord(c[1],None,None))
	x = np.array([(xd.measure(p,pp)["distanceKm"]) for p in pts])
	return np.sum(x**2)


def tryminimizeND(pinitial,xset):
	def f(c,pts):
		pp = GeoPoint(LatCoord(c[0],None,None),LongCoord(c[1],None,None))
		x = np.array([(xd.measure(p,pp)["distanceKm"]) for p in pts])
		return np.sum(x**2)
	xset = [p[1] for p in xset.values()]
	initial_guess = np.array(pinitial.astuple())
	#res = leastsq(f,initial_guess,args = (sourcep,sourced,))
	res = minimize(f,initial_guess,args=(xset,),options=dict(maxiter=20000))
	pres = GeoPoint(*[float(x) for x in res.x])
	return dict(pos=pres,distanceKm=xd.measure(pres,lt)["distanceKm"],)

def tryminimizePD(pinitial,xset):
	def fPD(c,pts):
		pp = GeoPoint(LatCoord(c[0],None,None),LongCoord(c[1],None,None))
		x = np.array([(xd.measure(p[1],pp)["distanceKm"]-p[0]) for p in pts])
		return np.sum(x**2)

	xset = [p for p in xset.values()]
	initial_guess = np.array(pinitial.astuple())
	#res = leastsq(f,initial_guess,args = (sourcep,sourced,))
	res = minimize(fPD,initial_guess,args=(xset,),options=dict(maxiter=20000))
	pres = GeoPoint(*[float(x) for x in res.x])
	return dict(pos=pres,distanceKm=xd.measure(pres,lt)["distanceKm"],)

def tryminimizeRD(pinitial,xset):
	def fRD(c,pts):
		pp = GeoPoint(LatCoord(c[0],None,None),LongCoord(c[1],None,None))
		x = np.array([(xd.measure(p[1],pp)["distanceKm"]-p[0]) for p in pts])
		return np.sum(x**2)

	xset = [(xd.measure(p[1],lt)["distanceKm"],p[1]) for p in xset.values()]
	print xset
	initial_guess = np.array(pinitial.astuple())
	#res = leastsq(f,initial_guess,args = (sourcep,sourced,))
	res = minimize(fRD,initial_guess,args=(xset,))
	pres = GeoPoint(*[float(x) for x in res.x])
	return dict(pos=pres,distanceKm=xd.measure(pres,lt)["distanceKm"],)



if __name__ == "__main__":
	# TODO pole2 set
	nolondon = False
	secondset = False
	kml = False
	initnearest = False

	for x in sys.argv[1:]:
		if x == "--help":
			print "options: nolondon secondset kml initnearest"
			sys.exit(0)
		elif x == "nolondon":
			nolondon = True
			print "removing london"
		elif x == "secondset":
			secondset = True
			print "second set"
		elif x == "kml":
			kml = True
			print "produce kml"
		elif x == "initnearest":
			initnearest = True

	# second pole set
	if secondset:
		distances1 = distancesB1
		distances2 = distancesB2
		locations["locpole"] = locations["locpole2"]
		outkml = "outB.kml"
	else:
		outkml = "out.kml"

	if nolondon:
		del distances2["London"]

	xd = xdistance()

	# Dump Locations
	glocations ={}
	for k,v in locations.items():
		a = v.split(",")
		glocations[k] = GeoPoint(a[0].strip(),a[1].strip())
		print k,glocations[k]

	# Pick references
	l1 = glocations["village"]
	l2 = glocations["photo"]
	l3 = glocations["locpole"]

	if kml:
		kmlpush(l1,"Village","pushpin")
		kmlpush(l2,"Photo GPS")
		kmlpush(l1,"Pole Location","pushpin")

	lt = l3 # use the pole as target

	# Info
	print "google:",l1,"photo:",l2,"distanceKm",xd.measure(l1,l2)["distanceKm"]
	print "google:",l1,"locpole:",l3,"distanceKm",xd.measure(l1,l3)["distanceKm"]
	print "------listing-------"
	print "location,reported,measured zara,measured photo"


	print "----first photo-------"
	w = 0
	for x,d0 in distances1.items():
		l = glocations[x]
		d = xd.measure(l,lt)["distanceKm"]
		print d0,d,d0-d,(d0-d)/d*100,l.lat.value,l.lat.value
		w = w + math.fabs(d0-d)

	print "----second photo------"
	for x,d0 in distances2.items():
		l = glocations[x]
		d = xd.measure(l,lt)["distanceKm"]
		print d0,d,d0-d,(d0-d)*100.0/d,l.lat.value,l.lat.value
		w = w + math.fabs(d0-d)

	print "mean",w/(len(distances1)+len(distances2))
	#print "set1 lookup"

	print "set1 is set from first photo"
	print "set2 is set from second photo"
	pts1 = dict([(k,(v,glocations[k])) for k,v in distances1.items()])
	pts2 = dict([(k,(v,glocations[k])) for k,v in distances2.items()])
	ptsa = dict(pts1.items()+pts2.items())
	if kml:
		for k,v in ptsa.items():
			kmlpush(v[1],k)

	print "----use trilateration to find optimal---"
	allp,bestinfo = findbesttrilat(ptsa,lt)

	# Optional write to KML
	if kml:
		for curset,p,d in allp:			
			kmlpush(p,"trilat " + " ".join(curset))


	if bestinfo[0] is not None:
		if kml:
			kmlpush(bestinfo[0],"Best of trilateration")
		print bestinfo
		print "best all distanceKm: ",xd.measure(bestinfo[0],lt)["distanceKm"]

	cogp = cog([x[1] for x in allp])
	print "cog distanceKm ",xd.measure(cogp,lt)["distanceKm"]

	if kml:
		kmlpush(cogp,"Cog of trilateration")

	# long
	if True:
		pinitial_guess = bestinfo[0]
		initial_guess = np.array(pinitial_guess.astuple())
		alltrilat = [x[1] for x in allp]
		print "alltrilat",len(alltrilat),alltrilat[0],"from",initial_guess

		res = minimize(faveraginglat,initial_guess,args=(alltrilat,))#,options=dict(maxiter=20000))
		pres = GeoPoint(*[float(x) for x in res.x])
		print "average of trilat result",pres," distanceKm",xd.measure(pres,lt)["distanceKm"]

		if kml:
			kmlpush(pres,"Average point of trilateration")

	if initnearest:
		print "***using nearest place for initialization per set***"
		pinitial = (
			ptsa[min([k for k in ptsa.keys()],key=lambda x: ptsa[k][0])][1],
			pts1[min([k for k in pts1.keys()],key=lambda x: pts1[k][0])][1],
			pts2[min([k for k in pts2.keys()],key=lambda x: pts2[k][0])][1],
		)
	else:
		print "***using best of trilaterationf as initi"
		pinitial = (bestinfo[0],bestinfo[0],bestinfo[0])

	print pinitial

	print "------Using position only"

	r1 = tryminimizeND(pinitial[1],pts1)
	r2 = tryminimizeND(pinitial[2],pts2)
	ra = tryminimizeND(pinitial[0],ptsa)

	print "\tinitial",pinitial,"distanceKm is",xd.measure(pinitial[0],lt)["distanceKm"]
	print "\tset1",r1
	print "\tset2",r2
	print "\tall",ra

	if kml:
		kmlpush(r1[0],"Multilateration without Distance - Set 1")
		kmlpush(r2[0],"Multilateration without Distance - Set 2")
		kmlpush(ra[0],"Multilateration without Distance - All")


	print "------Using pole distance"

	r1 = tryminimizePD(pinitial[1],pts1)
	r2 = tryminimizePD(pinitial[2],pts2)
	ra = tryminimizePD(pinitial[0],ptsa)

	print "\tinitial",pinitial,"distanceKm is",xd.measure(pinitial[0],l1)["distanceKm"]
	print "\tset1",r1
	print "\tset2",r2
	print "\tall",ra

	if kml:
		kmlpush(r1[0],"Multilateration with Distance - Set 1","redpushpin")
		kmlpush(r2[0],"Multilateration with Distance - Set 2","redpushpin")
		kmlpush(ra[0],"Multilateration with Distance - All","redpushpin")


	print "------Using real distance"

	r1 = tryminimizeRD(pinitial[1],pts1)
	r2 = tryminimizeRD(pinitial[2],pts2)
	ra = tryminimizeRD(pinitial[0],ptsa)


	print "\tinitial",pinitial,"distanceKm is",xd.measure(pinitial[0],lt)["distanceKm"]
	print "\tset1",r1
	print "\tset2",r2
	print "\tall",ra

	# Flush the KML
	if kml:	
		kmlpush(None,None)