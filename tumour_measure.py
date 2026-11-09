# Take user input of a TIFF file or directory, find and parse Leica LAS AF widefield data organisation.
# Segment nuclei and GFP +ve area, measure other channels and make +/-ve call based on z-score.
# For a file input, finds associated channels, runs analysis and shows objects as Overlay. For a directory input, recursively searches for TIFF files and analyses all positions found.
# Takes 1st and 2nd parent directory names as experimental conditions.
#	by Richard Butler, Gurdon Institute Imaging Facility

import math as maths
import os, re
import xml.etree.ElementTree as ET

from ij import IJ, ImagePlus, ImageStack, Prefs

from ij.gui import Roi, ShapeRoi, Overlay
from ij.measure import Calibration, ResultsTable
from ij.process import ImageProcessor, FloatProcessor, ByteProcessor, AutoThresholder, Blitter, FloodFiller
from ij.plugin.filter import MaximumFinder, ThresholdToSelection, EDM

from java.awt import Color
from java.io import File
from javax.swing import JFileChooser


#ch00: Hoechst
#ch01: GFP
#ch02: HA
#ch03: EdU

posRegex = re.compile("(?<=(P\\s))\\d+(?=_)")
chanRegex = re.compile("(?<=(_ch))\\d{2}(?=.tif)")

minGFPMean = 1000
maxGFPskew = 10

eduZcall = 1.0
haZcall = 2.0


class Position:
	def __init__(self, p, c, str1, str2, path):
		self.p = p
		self.c = c
		self.pid = (str(str1),str(str2),p)
		self.path = path


def getChannels(positions, idd):
	clist = []
	for pos in positions:
		if pos.pid == idd:
			clist.append(pos)
	clist.sort(key=lambda pos:pos.c)
	return clist


def getPositions(root, positions):
	for fil in root.listFiles():
		if fil.isDirectory():
			getPositions(fil, positions)
		elif fil.getName().endswith(".tif"):
			p = parsePosition(fil)
			if p is not None:
				positions.append(p)


def parsePosition(fil):
	name = fil.getName()
	pmatch = re.search(posRegex, name)
	p = int(pmatch.group()) if pmatch else None
	cmatch = re.search(chanRegex, name)
	c = int(cmatch.group()) if cmatch else None
	if p is None or c is None:
		return None
	parent1 = fil.getParentFile()
	parent2 = parent1.getParentFile()
	pos = Position(p,c, parent2.getName(), parent1.getName(), fil.getAbsolutePath())
	return pos


def getCalibration(imgpath):
	f = File(imgpath)
	pdir = f.getParentFile()
	xmlpath = os.path.join(pdir.getAbsolutePath(), "MetaData", "P 1.xml")
	doc = ET.parse(xmlpath)
	root = doc.getroot() # no namespace
	camerasettings = root.findall(".//ATLCameraSettingDefinition") #"ATLCameraSettingDefinition"
	pixelW = None
	pixelH = None
	for cameratag in camerasettings: # find a populated camera settings tag
		if "TheoCamSensorPixelSizeX" in cameratag.attrib:
			pixelW = float(cameratag.attrib["TheoCamSensorPixelSizeX"])
		if "TheoCamSensorPixelSizeY" in cameratag.attrib:
			pixelH = float(cameratag.attrib["TheoCamSensorPixelSizeY"])
		if pixelW and pixelH: break
	cal = Calibration()
	cal.pixelWidth = pixelW*1E5 # XML units are 10*microns
	cal.pixelHeight = pixelH*1E5
	cal.setUnit(u"\u00b5"+"m")
	return cal


def watershed(ip):
	TOL = 0.5
	floatEdm = EDM().makeFloatEDM(ip, 0, False)
	maxIp = MaximumFinder().findMaxima(floatEdm, TOL, ImageProcessor.NO_THRESHOLD, MaximumFinder.SEGMENTED, False, True)
	if (maxIp != None):
		ip.copyBits(maxIp, 0, 0, Blitter.AND)


def getNuclei(ip, s0, s1, minAPx, maxAPx, method):
	proc = ip.duplicate()
	proc.blurGaussian(s0)
	sub = ip.duplicate()
	sub.blurGaussian(s1)
	proc.copyBits(sub, 0,0, Blitter.SUBTRACT)
	stats = proc.getStatistics()
	hist = [int(n) for n in stats.getHistogram()]
	thresh = AutoThresholder().getThreshold(method, hist)
	thresh = (thresh/float(255)) * (stats.max-stats.min) + stats.min
	proc.threshold( int(thresh) )
	mask = proc.convertToByte(False)
	watershed(mask)
	mask.setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE)
	composite = ThresholdToSelection().convert(mask)
	rois = ShapeRoi(composite).getRois()
	nuclei = []
	for roi in rois:
		if minAPx <= roi.getStatistics().area <= maxAPx:
			nuclei.append(roi)
	return nuclei


def getSignalArea(ip, sigma, method):
	proc = ip.duplicate()
	proc.blurGaussian(sigma)
	stats = proc.getStatistics()
	hist = [int(n) for n in stats.getHistogram()]
	thresh = AutoThresholder().getThreshold(method, hist)
	thresh = (thresh/float(255)) * (stats.max-stats.min) + stats.min
	proc.threshold( int(thresh) )
	mask = proc.convertToByte(False)
	mask.setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE)
	roi = ThresholdToSelection().convert(mask)
	return ShapeRoi(roi)


def scorePositive(col, callZ):
	values = []
	for row in range(rt.getCounter()):
		values.append(rt.getValue(col, row))
	mean = sum(values)/float(len(values))
	ssd = sum( (v-mean)**2 for v in values)
	sd = maths.sqrt( ssd/float(len(values)-1) )
	for row in range(rt.getCounter()):
		z = (values[row]-mean)/sd
		rt.setValue(col+" Z", row, z)
		rt.setValue(col+" +ve", row, "Yes" if z>=callZ else "No")


def analyse(imp, pid, show):
	cal = imp.getCalibration()
	width = imp.getWidth()
	height = imp.getHeight()
	stack = imp.getStack()
	ol = Overlay()
	sigma = 1.0 # µm
	k = 5.
	nucsigma0 = sigma/cal.pixelWidth
	nucsigma1 = k*sigma/cal.pixelWidth

	areaSigma = 3.0 # µm
	areaSigmaPx = areaSigma/cal.pixelWidth
	
	minA = 20. # µm²
	maxA = 99999.
	minAPx = minA/(cal.pixelWidth*cal.pixelHeight)
	maxAPx = maxA/(cal.pixelWidth*cal.pixelHeight)
	nucip = stack.getProcessor(1)
	gfpip = stack.getProcessor(2)
	haip = stack.getProcessor(3)
	eduip = stack.getProcessor(4)
	
	nuclei = getNuclei(nucip, nucsigma0, nucsigma1, minAPx, maxAPx, AutoThresholder.Method.Huang)

	gfpStats = gfpip.getStatistics()
	hasGFP = gfpStats.mean >= minGFPMean and gfpStats.skewness <= maxGFPskew
	
	gfpRoi = getSignalArea(gfpip, areaSigmaPx, AutoThresholder.Method.Otsu) if hasGFP else None
	haRoi = getSignalArea(haip, areaSigmaPx, AutoThresholder.Method.Otsu)
	if show:
		for nuc in nuclei:
			nuc.setStrokeColor(Color.CYAN)
			ol.add(nuc)
		if gfpRoi is not None:
			gfpRoi.setStrokeColor(Color.MAGENTA)
			ol.add(gfpRoi)
		if haRoi is not None:
			haRoi.setStrokeColor(Color.YELLOW)
			ol.add(haRoi)
		imp.setOverlay(ol)

	for i,nuc in enumerate(nuclei):
		coord = nuc.getContourCentroid()
		npx = nuc.getStatistics().area
		gfpPos = gfpRoi is not None and gfpRoi.contains(int(coord[0]), int(coord[1]))
		gfpip.setRoi(nuc)
		gfp = gfpip.getStatistics().mean
		haip.setRoi(nuc)
		ha = haip.getStatistics().mean
		eduip.setRoi(nuc)
		edu = eduip.getStatistics().mean
		row = rt.getCounter()
		rt.setValue("Label", row, pid[0])
		rt.setValue("Treatment", row, pid[1])
		rt.setValue("P", row, pid[2])
		rt.setValue("Nucleus", row, i)
		rt.setValue("X", row, coord[0]*cal.pixelWidth)
		rt.setValue("Y", row, coord[1]*cal.pixelHeight)
		rt.setValue("GFP", row, gfp)
		rt.setValue("GFP +ve", row, "Yes" if gfpPos else "No")
		rt.setValue("HA", row, ha)
		rt.setValue("EdU", row, edu)
	
	srow = summaryrt.getCounter()
	summaryrt.setValue("Label", srow, pid[0])
	summaryrt.setValue("Treatment", srow, pid[1])
	summaryrt.setValue("P", srow, pid[2])
	summaryrt.setValue("GFP Area", srow, "0" if gfpRoi is None else gfpRoi.getStatistics().area*cal.pixelWidth*cal.pixelHeight)
	summaryrt.setValue("HA Area", srow, "0" if haRoi is None else haRoi.getStatistics().area*cal.pixelWidth*cal.pixelHeight)
	if gfpRoi and haRoi:
		andRoi = ShapeRoi(gfpRoi).and(ShapeRoi(haRoi))
		summaryrt.setValue("GFP+HA Area", srow, andRoi.getStatistics().area*cal.pixelWidth*cal.pixelHeight)


janicefile = None
path = Prefs.get("janice.path", os.path.expanduser("~"))
chooser = JFileChooser(path)
chooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES)
chooser.setSelectedFile(File(path))
op = chooser.showOpenDialog(None)
if op == JFileChooser.APPROVE_OPTION:
	sel = chooser.getSelectedFile()
	Prefs.set("janice.path", sel.getAbsolutePath())
	janicefile = sel
else:
	exit("cancelled")

positions = []
show = False
if janicefile.isDirectory():
	getPositions(janicefile, positions)
else:
	if janicefile.getName().endswith(".tif"):
		selpos = parsePosition(janicefile)
		parent = janicefile.getParentFile()
		getPositions(parent, positions)
		positions = [pos for pos in positions if pos.pid == selpos.pid] # keep only associated channels if a File was selected
		show = True
	else:
		exit(0)

idset = set()
for pos in positions:
	if pos.pid not in idset:
		idset.add(pos.pid)
idset = list(idset)
idset.sort()

rt = ResultsTable()
summaryrt = ResultsTable()
for pid in idset:
	pidstr = "{} {} P{}".format(*pid)
	IJ.showStatus("Analysing "+pidstr)
	channels = getChannels(positions, pid)
	
	stack = ImageStack()
	cal = None
	for i,cpos in enumerate(channels):
		if cal is None:
			cal = getCalibration(cpos.path)
		chimp = ImagePlus(cpos.path)
		chip = chimp.getProcessor()
		stack.addSlice("C{}".format(i), chip)

	imp = ImagePlus(pidstr, stack)
	imp.setDimensions(stack.getSize(),1,1)
	imp.setCalibration(cal)
	analyse(imp, pid, show)
	if show:
		imp.show()
	else:
		imp.close()

if rt.getCounter()>0:
	scorePositive("HA", haZcall)
	scorePositive("EdU", eduZcall)
	
	for srow,pid in enumerate(idset):
		indices = []
		ncells = 0
		ngfp = 0
		nha = 0
		ngfpha = 0
		nedu = 0
		ngfpedu = 0
		for row in range(rt.getCounter()):
			label = rt.getStringValue("Label", row)
			treatment = rt.getStringValue("Treatment", row)
			p = int(rt.getStringValue("P", row))
			if label == pid[0] and treatment == pid[1] and p == pid[2]:
				gfpPos = rt.getStringValue("GFP +ve", row) == "Yes"
				haPos = rt.getStringValue("HA +ve", row) == "Yes"
				eduPos = rt.getStringValue("EdU +ve", row) == "Yes"
				ncells += 1
				if gfpPos:
					ngfp += 1
					if haPos: ngfpha += 1
					if eduPos: ngfpedu += 1
				if haPos: nha += 1
				if eduPos: nedu += 1

		summaryrt.setValue("n Cells", srow, ncells)
		summaryrt.setValue("n GFP +ve", srow, ngfp)
		summaryrt.setValue("n HA +ve", srow, nha)
		summaryrt.setValue("n EdU +ve", srow, nedu)
		summaryrt.setValue("n GFP HA +/+", srow, ngfpha)
		summaryrt.setValue("n GFP EdU +/+", srow, ngfpedu)

	rt.show("Nuclei "+janicefile.getAbsolutePath())
	summaryrt.show("Results "+janicefile.getAbsolutePath())
	