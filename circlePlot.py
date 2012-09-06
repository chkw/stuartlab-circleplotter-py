#!/usr/bin/env python
"""circlePlot.py: 

Usage:
  circlePlot.py [options] outputDir inputFile [inputFile ...]

Options:
  -s str		list file containing samples to include
  -f str		list file containing features to include
  -o str		feature;file[,file ...] or feature
  -c str		file to use as center colors
  -l			print the feature identifier in the circle or not (default: FALSE)
  -q			run quietly
"""
## Written By: Steve Benz and Zack Sanborn
## Modified By: Sam Ng
## Last Updated: 10/10/2011
## Modified by: chrisw DEC 2011
import getopt, math, os, sys, re
from matplotlib import use

# specify a backend to use the Anti-Grain Geometry C++ library to make a raster (pixel) image
use('Agg')

from pylab import plt, axes, axis, fill, text, xlim, ylim, savefig, close
from random import random
import mData

import hashlib

from time import time

verbose = True

# named html colors to hex codes
htmlColorNamesDict = dict()
htmlColorNamesDict['aliceblue'] = '#F0F8FF'
htmlColorNamesDict['antiquewhite'] = '#FAEBD7'
htmlColorNamesDict['aqua'] = '#00FFFF'
htmlColorNamesDict['aquamarine'] = '#7FFFD4'
htmlColorNamesDict['azure'] = '#F0FFFF'
htmlColorNamesDict['beige'] = '#F5F5DC'
htmlColorNamesDict['bisque'] = '#FFE4C4'
htmlColorNamesDict['black'] = '#000000'
htmlColorNamesDict['blanchedalmond'] = '#FFEBCD'
htmlColorNamesDict['blue'] = '#0000FF'
htmlColorNamesDict['blueviolet'] = '#8A2BE2'
htmlColorNamesDict['brown'] = '#A52A2A'
htmlColorNamesDict['burlywood'] = '#DEB887'
htmlColorNamesDict['cadetblue'] = '#5F9EA0'
htmlColorNamesDict['chartreuse'] = '#7FFF00'
htmlColorNamesDict['chocolate'] = '#D2691E'
htmlColorNamesDict['coral'] = '#FF7F50'
htmlColorNamesDict['cornflowerblue'] = '#6495ED'
htmlColorNamesDict['cornsilk'] = '#FFF8DC'
htmlColorNamesDict['crimson'] = '#DC143C'
htmlColorNamesDict['cyan'] = '#00FFFF'
htmlColorNamesDict['darkblue'] = '#00008B'
htmlColorNamesDict['darkcyan'] = '#008B8B'
htmlColorNamesDict['darkgoldenrod'] = '#B8860B'
htmlColorNamesDict['darkgray'] = '#A9A9A9'
htmlColorNamesDict['darkgrey'] = '#A9A9A9'
htmlColorNamesDict['darkgreen'] = '#006400'
htmlColorNamesDict['darkkhaki'] = '#BDB76B'
htmlColorNamesDict['darkmagenta'] = '#8B008B'
htmlColorNamesDict['darkolivegreen'] = '#556B2F'
htmlColorNamesDict['darkorange'] = '#FF8C00'
htmlColorNamesDict['darkorchid'] = '#9932CC'
htmlColorNamesDict['darkred'] = '#8B0000'
htmlColorNamesDict['darksalmon'] = '#E9967A'
htmlColorNamesDict['darkseagreen'] = '#8FBC8F'
htmlColorNamesDict['darkslateblue'] = '#483D8B'
htmlColorNamesDict['darkslategray'] = '#2F4F4F'
htmlColorNamesDict['darkslategrey'] = '#2F4F4F'
htmlColorNamesDict['darkturquoise'] = '#00CED1'
htmlColorNamesDict['darkviolet'] = '#9400D3'
htmlColorNamesDict['deeppink'] = '#FF1493'
htmlColorNamesDict['deepskyblue'] = '#00BFFF'
htmlColorNamesDict['dimgray'] = '#696969'
htmlColorNamesDict['dimgrey'] = '#696969'
htmlColorNamesDict['dodgerblue'] = '#1E90FF'
htmlColorNamesDict['firebrick'] = '#B22222'
htmlColorNamesDict['floralwhite'] = '#FFFAF0'
htmlColorNamesDict['forestgreen'] = '#228B22'
htmlColorNamesDict['fuchsia'] = '#FF00FF'
htmlColorNamesDict['gainsboro'] = '#DCDCDC'
htmlColorNamesDict['ghostwhite'] = '#F8F8FF'
htmlColorNamesDict['gold'] = '#FFD700'
htmlColorNamesDict['goldenrod'] = '#DAA520'
htmlColorNamesDict['gray'] = '#808080'
htmlColorNamesDict['grey'] = '#808080'
htmlColorNamesDict['green'] = '#008000'
htmlColorNamesDict['greenyellow'] = '#ADFF2F'
htmlColorNamesDict['honeydew'] = '#F0FFF0'
htmlColorNamesDict['hotpink'] = '#FF69B4'
htmlColorNamesDict['indianred'] = '#CD5C5C'
htmlColorNamesDict['indigo'] = '#4B0082'
htmlColorNamesDict['ivory'] = '#FFFFF0'
htmlColorNamesDict['khaki'] = '#F0E68C'
htmlColorNamesDict['lavender'] = '#E6E6FA'
htmlColorNamesDict['lavenderblush'] = '#FFF0F5'
htmlColorNamesDict['lawngreen'] = '#7CFC00'
htmlColorNamesDict['lemonchiffon'] = '#FFFACD'
htmlColorNamesDict['lightblue'] = '#ADD8E6'
htmlColorNamesDict['lightcoral'] = '#F08080'
htmlColorNamesDict['lightcyan'] = '#E0FFFF'
htmlColorNamesDict['lightgoldenrodyellow'] = '#FAFAD2'
htmlColorNamesDict['lightgray'] = '#D3D3D3'
htmlColorNamesDict['lightgrey'] = '#D3D3D3'
htmlColorNamesDict['lightgreen'] = '#90EE90'
htmlColorNamesDict['lightpink'] = '#FFB6C1'
htmlColorNamesDict['lightsalmon'] = '#FFA07A'
htmlColorNamesDict['lightseagreen'] = '#20B2AA'
htmlColorNamesDict['lightskyblue'] = '#87CEFA'
htmlColorNamesDict['lightslategray'] = '#778899'
htmlColorNamesDict['lightslategrey'] = '#778899'
htmlColorNamesDict['lightsteelblue'] = '#B0C4DE'
htmlColorNamesDict['lightyellow'] = '#FFFFE0'
htmlColorNamesDict['lime'] = '#00FF00'
htmlColorNamesDict['limegreen'] = '#32CD32'
htmlColorNamesDict['linen'] = '#FAF0E6'
htmlColorNamesDict['magenta'] = '#FF00FF'
htmlColorNamesDict['maroon'] = '#800000'
htmlColorNamesDict['mediumaquamarine'] = '#66CDAA'
htmlColorNamesDict['mediumblue'] = '#0000CD'
htmlColorNamesDict['mediumorchid'] = '#BA55D3'
htmlColorNamesDict['mediumpurple'] = '#9370D8'
htmlColorNamesDict['mediumseagreen'] = '#3CB371'
htmlColorNamesDict['mediumslateblue'] = '#7B68EE'
htmlColorNamesDict['mediumspringgreen'] = '#00FA9A'
htmlColorNamesDict['mediumturquoise'] = '#48D1CC'
htmlColorNamesDict['mediumvioletred'] = '#C71585'
htmlColorNamesDict['midnightblue'] = '#191970'
htmlColorNamesDict['mintcream'] = '#F5FFFA'
htmlColorNamesDict['mistyrose'] = '#FFE4E1'
htmlColorNamesDict['moccasin'] = '#FFE4B5'
htmlColorNamesDict['navajowhite'] = '#FFDEAD'
htmlColorNamesDict['navy'] = '#000080'
htmlColorNamesDict['oldlace'] = '#FDF5E6'
htmlColorNamesDict['olive'] = '#808000'
htmlColorNamesDict['olivedrab'] = '#6B8E23'
htmlColorNamesDict['orange'] = '#FFA500'
htmlColorNamesDict['orangered'] = '#FF4500'
htmlColorNamesDict['orchid'] = '#DA70D6'
htmlColorNamesDict['palegoldenrod'] = '#EEE8AA'
htmlColorNamesDict['palegreen'] = '#98FB98'
htmlColorNamesDict['paleturquoise'] = '#AFEEEE'
htmlColorNamesDict['palevioletred'] = '#D87093'
htmlColorNamesDict['papayawhip'] = '#FFEFD5'
htmlColorNamesDict['peachpuff'] = '#FFDAB9'
htmlColorNamesDict['peru'] = '#CD853F'
htmlColorNamesDict['pink'] = '#FFC0CB'
htmlColorNamesDict['plum'] = '#DDA0DD'
htmlColorNamesDict['powderblue'] = '#B0E0E6'
htmlColorNamesDict['purple'] = '#800080'
htmlColorNamesDict['red'] = '#FF0000'
htmlColorNamesDict['rosybrown'] = '#BC8F8F'
htmlColorNamesDict['royalblue'] = '#4169E1'
htmlColorNamesDict['saddlebrown'] = '#8B4513'
htmlColorNamesDict['salmon'] = '#FA8072'
htmlColorNamesDict['sandybrown'] = '#F4A460'
htmlColorNamesDict['seagreen'] = '#2E8B57'
htmlColorNamesDict['seashell'] = '#FFF5EE'
htmlColorNamesDict['sienna'] = '#A0522D'
htmlColorNamesDict['silver'] = '#C0C0C0'
htmlColorNamesDict['skyblue'] = '#87CEEB'
htmlColorNamesDict['slateblue'] = '#6A5ACD'
htmlColorNamesDict['slategray'] = '#708090'
htmlColorNamesDict['slategrey'] = '#708090'
htmlColorNamesDict['snow'] = '#FFFAFA'
htmlColorNamesDict['springgreen'] = '#00FF7F'
htmlColorNamesDict['steelblue'] = '#4682B4'
htmlColorNamesDict['tan'] = '#D2B48C'
htmlColorNamesDict['teal'] = '#008080'
htmlColorNamesDict['thistle'] = '#D8BFD8'
htmlColorNamesDict['tomato'] = '#FF6347'
htmlColorNamesDict['turquoise'] = '#40E0D0'
htmlColorNamesDict['violet'] = '#EE82EE'
htmlColorNamesDict['wheat'] = '#F5DEB3'
htmlColorNamesDict['white'] = '#FFFFFF'
htmlColorNamesDict['whitesmoke'] = '#F5F5F5'
htmlColorNamesDict['yellow'] = '#FFFF00'
htmlColorNamesDict['yellowgreen'] = '#9ACD32'

class rgb:
	"""Object that represents an RGB color code."""
	def __init__(self, r, g, b):
		"""Initialize an object, ensuring that RGB values are within acceptable the acceptable range of 0 to 255."""
		self.r = int(round(r))
		self.g = int(round(g))
		self.b = int(round(b))

		if self.r > 255:
			self.r = 255
		elif self.r < 0:
			self.r = 0
		if self.g > 255:
			self.g = 255
		elif self.g < 0:
			self.g = 0
		if self.b > 255:
			self.b = 255
		elif self.b < 0:
			self.b = 0

	def tohex(self):
		"""Convert RGB values to a hex color code."""
		r = self.r
		g = self.g
		b = self.b
		hexchars = "0123456789ABCDEF"
		return "#" + hexchars[r / 16] + hexchars[r % 16] + hexchars[g / 16] + hexchars[g % 16] + hexchars[b / 16] + hexchars[b % 16]

def usage(code=0):
	"""Print docs."""
	print __doc__
	if code != None: sys.exit(code)

def log(msg, die=False):
	"""Perform logging to sys.stderr."""
	if verbose:
		sys.stderr.write(msg)
	if die:
		sys.exit(1)

def syscmd(cmd):
	"""Execute a system command via os.system(cmd)."""
	log("running:\n\t" + cmd + "\n")
	exitstatus = os.system(cmd)
	if exitstatus != 0:
		print "Failed with exit status %i" % exitstatus
		sys.exit(10)
	log("... done\n")

def scmp(a, b, feature, dataList):
	"""Comparison function for sorting the samples of a specified feature."""
	dataFeature = feature
	if (a not in dataList[0]) & (b in dataList[0]):
		return(1)
	elif (a in dataList[0]) & (b not in dataList[0]):
		return(-1)
	elif (b not in dataList[0]) & (a not in dataList[0]):
		return(0)
	if dataFeature not in dataList[0][a]:
		if "*" in dataList[0][a]:
			dataFeature = "*"
		else:
			return(0)
	val = cmp(dataList[0][a][dataFeature], dataList[0][b][dataFeature])
	if val == 0:
		if len(dataList) > 1:
			# recursively check subsequent datasets until some ordering is determined
			val = scmp(a, b, feature, dataList[1:])
		else:
			return(0)
	return(val)

def polar(r, val):
	"""Convert polar coordinates to cartesian coordinates where 12-o-clock is 0 degrees and increases in the clockwise direction."""
	theta = -2.0 * math.pi * val + math.pi / 2.0
	x = r * math.cos(theta)
	y = r * math.sin(theta)
	return x, y

def getColor(val, minVal, maxVal, minColor=rgb(0, 0, 255), zeroColor=rgb(255, 255, 255), maxColor=rgb(255, 0, 0), purple0Hack = False):
	"""Find the hex color code for a value via linear interpolation between minVal and maxVal."""
	# check if val is a number
	fval = None
	try:
		fval = float(val)
		if fval != fval:
			raise ValueError
	except ValueError:
		col = greyRGB;
		return col.tohex()

	# perform linear interpolation, guaranteeing the result in [-1,1]
	if fval < 0.0:
		col = minColor
		if fval < minVal:
			fval = -1.0
		else:
			fval = fval / minVal
	# PURPLE HACK
	elif purple0Hack == True and fval == 0.0 :
		return violetRGB.tohex()
	else:
		col = maxColor
		if fval > maxVal:
			fval = 1.0
		else:
			fval = fval / maxVal

	r = fval * float(col.r - zeroColor.r) + zeroColor.r
	g = fval * float(col.g - zeroColor.g) + zeroColor.g
	b = fval * float(col.b - zeroColor.b) + zeroColor.b
	col = rgb(r, g, b)

	return col.tohex()

def plotScale(imgFile, minVal, maxVal):
	imgSize = (2, 4)
	fig = plt.figure(figsize=imgSize, dpi=100, frameon=True, facecolor='w')
	for i in xrange(10):
		val = minVal + i * (maxVal - minVal) / 10
		col = getColor(val, minVal, maxVal)
		X = [float(i) / 10, float(i + 1) / 10, float(i + 1) / 10, float(i) / 10, float(i) / 10]
		Y = [1, 1, 0, 0, 1]
		fill(X, Y, col, lw=1, ec=col)
	savefig(imgFile)
	close()

def plotCircle(imgFile, label, centerColHex=rgb(255, 255, 255).tohex(), circleCols=[[rgb(200, 200, 200).tohex()]], innerRadTotal=0.2, outerRadTotal=0.5, width=5, tstep=0.01):
	"""Plot and save a circlePlot image using matplotlib module."""
	## image settings
	imgSize = (width, width)
	fig = plt.figure(figsize=imgSize, dpi=100, frameon=True, facecolor='w')
	axes([0, 0, 1, 1], frameon=True, axisbg='w')
	axis('off')
	circleWid = (outerRadTotal - innerRadTotal) / float(len(circleCols))

	## color center
	outerRadCenter = innerRadTotal
	outerRadCenter -= .01
	X = []
	Y = []
	x, y = polar(outerRadCenter, 0)
	X.append(x)
	Y.append(y)
	ti = 0
	while ti < 1:
		x, y = polar(outerRadCenter, ti)
		X.append(x)
		Y.append(y)
		ti += tstep
		if ti > 1:
			break
	x, y = polar(outerRadCenter, 1)
	X.append(x)
	Y.append(y)

	if centerColHex == None:
		# transparent center
		fill(X, Y, rgb(255, 255, 255).tohex(), lw=1, ec='none', fill=False)
	else:
		# color-filled center
		fill(X, Y, centerColHex, lw=1, ec=centerColHex)

	time0 = time()

	## color rings
	# this part is slow ~0.6 sec for one dataset (536 samples)
	for i in xrange(len(circleCols)):
		innerRadRing = (i * circleWid) + innerRadTotal
		outerRadRing = ((i + 1) * circleWid) + innerRadTotal - .01
		for j in xrange(len(circleCols[i])):
			t0 = float(j) / len(circleCols[i])
			t1 = float(j + 1) / len(circleCols[i])
			X = []
			Y = []
			x, y = polar(innerRadRing, t0)
			X.append(x)
			Y.append(y)
			ti = t0
			while ti < t1:
				x, y = polar(outerRadRing, ti)
				X.append(x)
				Y.append(y)
				ti += tstep
				if ti > t1:
					break
			x, y = polar(outerRadRing, t1)
			X.append(x)
			Y.append(y)
			ti = t1
			while ti > t0:
				x, y = polar(innerRadRing, ti)
				X.append(x)
				Y.append(y)
				ti -= tstep
				if ti < t0:
					break
			x, y = polar(innerRadRing, t0)
			X.append(x)
			Y.append(y)
			fill(X, Y, circleCols[i][j], lw=1, ec=circleCols[i][j])

#	log("%s to get ring colors\n" % (time() - time0))
	time0 = time()
	
	## save image
	text(0, 0, label, ha='center', va='center', size='xx-large')
	xlim(-0.5, 0.5)
	ylim(-0.5, 0.5)
	savefig(imgFile, transparent=True)
	close()

#	log("%s to savefig\n" % (time() - time0))
	time0 = time()

def getCohortMinMaxValues(featureList, sampleList, circleData):
	"""Get the minVal and maxVal of sample scores among the specified featureList for the ring/dataset."""
	minValList = []
	maxValList = []

	for ring in xrange(len(circleData)):
		ringVals = []

		# get ring values in effort to find min/max values for each *ring*
		for sample in sampleList:
			if sample in circleData[ring]:
				for feature in featureList:
					if feature in circleData[ring][sample]:
						ringVals.append(circleData[ring][sample][feature])
					elif "*" in circleData[ring][sample]:
						ringVals.append(circleData[ring][sample]["*"])
	
		# find the min & max sample scores for this ring in this feature
		floatList = mData.floatList(ringVals)

		minValList.append(min([-0.01] + floatList))
		maxValList.append(max([0.01] + floatList))

	return (minValList, maxValList)
	
def drawCircleImageForFeature(feature, samples, label, imgFile, circleData, circleColors, centerColHex=None, width=5, minValList=None, maxValList=None, purple0Hack=False):
	"""Draw a circle map image and write it to a file."""
	# feature - feature to draw image for. This is some kind of concept: for example, a gene.
	# samples - sample names of data
	# label - label to use in image
	# imgFile - file object to which image will be written
	# circleData - data struct containing sample data for features.  It is a list of dict[col][row]=score .
	# circleColors - a list of (minColor),(zeroColor),(maxColor)
	# centerColHex - hex code for center color fill. If none, then make transparent center.

	# centerCol is the color of the center of the circleImage
#	centerCol = whiteRGB.tohex()
	
	# circleCols is a list.  Each member of the list represents a list of colors in a ring.
	circleCols = []

	# iterate through rings of data
	for ring in xrange(len(circleData)):
		ringCols = []

		# get minVal and maxVal
		minVal = None
		maxVal = None
		if minValList == None or maxValList == None:
			ringVals = []
			
			# get ring values in effort to find min/max values for each *ring*
			for sample in samples:
				if sample in circleData[ring]:
					if feature in circleData[ring][sample]:
						ringVals.append(circleData[ring][sample][feature])
					elif "*" in circleData[ring][sample]:
						ringVals.append(circleData[ring][sample]["*"])

			# find the min & max sample scores for this ring in this feature
			floatList = mData.floatList(ringVals)
			minVal = min([-0.01] + floatList)
			maxVal = max([0.01] + floatList)
		else:
			minVal = minValList[ring]
			maxVal = maxValList[ring]

		# convert scores into colors
		for sample in samples:
			if sample in circleData[ring]:
				if feature in circleData[ring][sample]:
					ringCols.append(getColor(circleData[ring][sample][feature], minVal, maxVal, minColor=circleColors[ring][0], zeroColor=circleColors[ring][1], maxColor=circleColors[ring][2], purple0Hack=purple0Hack))
				elif "*" in circleData[ring][sample]:
					ringCols.append(getColor(circleData[ring][sample]["*"], minVal, maxVal, minColor=circleColors[ring][0], zeroColor=circleColors[ring][1], maxColor=circleColors[ring][2], purple0Hack=purple0Hack))
				else:
					# sample exists, but no score for the feature
					ringCols.append(greyRGB.tohex())
			else:
				# this sample not found in the sample data
				ringCols.append(greyRGB.tohex())

		# add the ring
		circleCols.append(ringCols)

	# plot the image
	plotCircle(imgFile, label=label, centerColHex=centerColHex, circleCols=circleCols, innerRadTotal=0.2, outerRadTotal=0.5, width=width)

def getHashedImageFileName(prefix="", feature=""):
	"""Get a hashed image filename the is created from the parameters of the image plotting."""
	hasher = hashlib.md5()
	hasher.update(" ".join([prefix, feature]))
	return hasher.hexdigest()

def cgi_routine(outputDir, dataMatrixNameList, circleData, samples, features, orderFeature, matrixPriorityList=None, printLabel=True, cohortMinMax=False):
	"""Routine for program execution via cgi. Returns a dictionary describing the fileName for each feature."""
	# circleData is a list of dict[col][row]=score from each circleFile
	# The ordering of data in circlData must by synced with their names in the list
	# dataMatrixNameList is a list of rings in order of display
	# matrixPriorityList is a list of rings in order of sorting priority
	
	if matrixPriorityList == None:
		matrixPriorityList = dataMatrixNameList

	# circleColors is a list of (minColor),(zeroColor),(maxColor)
	# one item in the list per dataset
	circleColors = []

	for i in xrange(len(dataMatrixNameList)):
		minCol = blueRGB
		zerCol = whiteRGB
		maxCol = redRGB

		circleColors.append((minCol, zerCol, maxCol))
	# end section for reading circleFiles
		
	## sort
	# reorder circleData to reflect priorities in matrixPriorityList
	orderData = list()
	for priority in xrange(len(matrixPriorityList)):
		matrixName = matrixPriorityList[priority]
		for ring in xrange(len(dataMatrixNameList)):
			if dataMatrixNameList[ring] == matrixName:
				orderData.append(circleData[ring])
				continue

	# sort samples based on sample score in orderData
	samples.sort(lambda x, y: scmp(x, y, orderFeature, orderData))
	# end section for sample ordering

	# get min/max values for datasets
	if cohortMinMax:
		(minValList, maxValList) = getCohortMinMaxValues(features, samples, circleData)
	else:
		(minValList, maxValList) = (None, None)

	## plot images
	outputFilesDict = dict()

	prefix = " ".join([" ".join(dataMatrixNameList), " ".join(samples), " ".join(matrixPriorityList), " ".join(str(cohortMinMax))])

	for feature in features:
		sanitizedFeatureName = getHashedImageFileName(prefix=prefix, feature=feature)
#		sanitizedFeatureName = re.sub("[/:]", "_", feature)
		imgFileName = "%s.png" % (sanitizedFeatureName)
		imgFilePath = "%s/%s" % (outputDir, imgFileName)

		if os.path.exists(imgFilePath):
			# don't plot if file exists
			outputFilesDict[feature] = imgFileName
			continue

		imgFile = "%s" % (imgFilePath)
		label = ""
		if printLabel:
			label = feature
		drawCircleImageForFeature(feature, samples, label, imgFile, circleData, circleColors, centerColHex=None, width=1.5, minValList=minValList, maxValList=maxValList)

		outputFilesDict[feature] = imgFileName

	return (outputFilesDict, samples)

def cli_routine(outputDir, circleFiles, orderFiles, sampleFile, featureFile, orderFeature, centerFile, printLabel, verbose, cohortMinMax=False, purpleHack = True):
	"""Routine for program execution via command-line."""
	# I've tried not to touch this method as much as possible.
	# I don't want to break the way it was working for Sam Ng.
	# chrisw
	
	## execute
	samples = []
	features = []
	if sampleFile != None:
		samples = mData.rList(sampleFile)
	if featureFile != None:
		features = mData.rList(featureFile)
	# end section for getting lists of samples and features
	
	## read circleFiles
	# circleData is a list of dict[col][row]=score from each circleFile
	circleData = []
	# circleColorsPalette is a list of (minColor),(zeroColor),(maxColor)
	circleColorsPalette = []

	for i in xrange(len(circleFiles)):
		# get data, samples, and features from each circleFile
		# data is a dict[col][row]=score
		# cols is a list of sample names
		# features is a list of feature names
		(data, cols, rows) = mData.rCRSData(circleFiles[i], retFeatures=True)
		circleData.append(data)
		minCol = lightBlueRGB
		zerCol = whiteRGB
		maxCol = redRGB

		# special cases for -meth and -mut
#		if circleFiles[i].endswith("meth"):
#			maxCol = blueRGB
#			minCol = redRGB
#			log("Color: meth\n")
#		elif circleFiles[i].endswith("mut"):
#			maxCol = blackRGB
#			minCol = whiteRGB
#			log("Color: mut\n")

		circleColorsPalette.append((minCol, zerCol, maxCol))

		# if no sampleFile/featureFile, default to using samples/features from circleFiles
		if sampleFile == None:
			samples = list(set(cols) | set(samples))
		if featureFile == None:
			features = list(set(rows) | set(features))
	# end section for reading circleFiles

	## read centerFile
	centerData = None
	if centerFile != None:
		centerData = mData.r2Col(centerFile, header=True)

	## sort
	if orderFeature != None:
		if len(orderFiles) > 0:
			orderData = []
			orderColors = []
			for i in xrange(len(orderFiles)):
				orderData.append(mData.rCRSData(orderFiles[i]))
				minCol = whiteRGB
				zerCol = whiteRGB
				maxCol = blackRGB
				orderColors.append((minCol, zerCol, maxCol))
		else:
			orderData = circleData

		# sort samples based on sample score in orderData
		# priority of sorting determined by orderFiles parameter
		samples.sort(lambda x, y: scmp(x, y, orderFeature, orderData))

		## cohort png
		# cgi will probably not use orderFiles
		if len(orderFiles) > 0:
			imgFile = "%s/Cohort.png" % (outputDir)
			label = "Cohort"
			centerCol = whiteRGB.tohex()
			cohortCircleCols = []
			for i in xrange(len(orderData)):
				ringCols = []
				ringVals = []
				for sample in samples:
					if sample in orderData[i]:
						if orderFeature in orderData[i][sample]:
							ringVals.append(orderData[i][sample][orderFeature])
						elif "*" in orderData[i][sample]:
							ringVals.append(orderData[i][sample]["*"])
				minVal = min([-0.01] + mData.floatList(ringVals))
				maxVal = max([0.01] + mData.floatList(ringVals))
				for sample in samples:
					if sample in orderData[i]:
						if orderFeature in orderData[i][sample]:
							ringCols.append(getColor(orderData[i][sample][orderFeature], minVal, maxVal, minColor=orderColors[i][0], zeroColor=orderColors[i][1], maxColor=orderColors[i][2]))
						elif "*" in orderData[i][sample]:
							ringCols.append(getColor(orderData[i][sample]["*"], minVal, maxVal, minColor=orderColors[i][0], zeroColor=orderColors[i][1], maxColor=orderColors[i][2]))
						else:
							ringCols.append(greyRGB.tohex())
					else:
						ringCols.append(greyRGB.tohex())
				cohortCircleCols.append(ringCols)
			plotCircle(imgFile, label=label, centerCol=centerCol, circleCols=cohortCircleCols, innerRadTotal=0.2, outerRadTotal=0.5, width=5)
	# end section for sample ordering

	## plot images
	if centerData != None:
		centerDataFloatList = mData.floatList(centerData.values())
		centerDataMinVal = min([-0.01] + centerDataFloatList)
		centerDataMaxVal = max([0.01] + centerDataFloatList)

	# get min/max values for datasets
	if cohortMinMax:
		(minValList, maxValList) = getCohortMinMaxValues(features, samples, circleData)
	else:
		(minValList, maxValList) = (None, None)

	for feature in features:
		log("Drawing %s\n" % (feature))
		centerColHex = None
		if centerData != None:
			if feature in centerData:
				centerColHex = getColor(centerData[feature], centerDataMinVal, centerDataMaxVal, minColor=lightBlueRGB, zeroColor=whiteRGB, purple0Hack=purpleHack)
				
		imgFile = "%s/%s.png" % (outputDir, re.sub("[/:]", "_", feature))

		label = ""
		if printLabel:
			label = feature

		image_width = 5.0

		drawCircleImageForFeature(feature, samples, label, imgFile, circleData, circleColorsPalette, width=image_width, centerColHex=centerColHex, minValList=minValList, maxValList=maxValList, purple0Hack=purpleHack)

	for sample in samples:
		log("ordered samples: %s\n" % (sample))

def main(args):
	# I've tried not to touch this method as much as possible.
	# I don't want to break the way it was working for Sam Ng.
	# chrisw

	## parse arguments
	try:
		opts, args = getopt.getopt(args, "s:f:o:c:lqm")
	except getopt.GetoptError, err:
		print str(err)
		usage(2)
	if len(args) < 2:
		usage(2)
	
	outputDir = args[0].rstrip("/")
	circleFiles = args[1:]
	
	sampleFile = None
	featureFile = None
	orderFeature = None
	centerFile = None
	printLabel = False
	cohortMinMax = False

	global verbose
	for o, a in opts:
		if o == "-s":
			sampleFile = a
		elif o == "-f":
			featureFile = a
		elif o == "-o":
			sa = re.split(";", a)
			if len(sa) == 1:
				orderFeature = sa[0]
				orderFiles = []
			else:
				orderFeature = sa[0]
				orderFiles = re.split(",", sa[1])
		elif o == "-c":
			centerFile = a
		elif o == "-l":
			printLabel = True
		elif o == "-q":
			verbose = False
		elif o == "-m":
			cohortMinMax = True
	# end section for parsing arguments
	
	# execute the routine for command-line usage
	cli_routine(outputDir, circleFiles, orderFiles, sampleFile, featureFile, orderFeature, centerFile, printLabel, verbose, cohortMinMax=cohortMinMax)

# some rgb colors
blueRGB = rgb(0, 0, 255)
lightBlueRGB = rgb(0, 100, 255)
whiteRGB = rgb(255, 255, 255)
redRGB = rgb(255, 0, 0)
blackRGB = rgb(0, 0, 0)
greyRGB = rgb(200, 200, 200)
purpleRGB = rgb(128, 0, 128)
violetRGB = rgb(143, 0, 255)

if __name__ == "__main__":
	main(sys.argv[1:])
