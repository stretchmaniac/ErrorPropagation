import math
from random import *
from scipy.stats import norm
import matplotlib.pyplot as plt
from matplotlib import collections as mc
import ast
import astor
import re

class Distribution:
	def __init__(self, intervals):
		self.intervals = intervals
		self.numPoints = len(intervals) + 1

	@staticmethod
	def normalDistribution(center, standardDev, numPoints):
		# there are alternatives to constant density in the x direction, but
		# 	for simplicity we'll stick to constant for now
		varience = standardDev**2
		minZ = -6
		maxZ = 6

		intervals = []

		def zTrans(z):
			return center + z*standardDev

		xStep = (zTrans(maxZ) - zTrans(minZ)) / (numPoints - 1)
		prevPoint = [zTrans(-6)-xStep, 0]

		for i in range(numPoints):
			x = zTrans(minZ) + (i+random()/10-.05) * xStep

			# formula comes from wikipedia article on gaussian distribution
			y = (1/math.sqrt(2*math.pi*varience))*math.e**(-(x-center)**2/(2*varience))

			intervals.append([prevPoint, [x,y]])
			prevPoint = [x, y]

		return Distribution(intervals)

	@staticmethod
	def vDistribution(leftRad, center, rightRad, numPoints):
		xDelta = (rightRad - leftRad) / numPoints

		def plot(x):
			if x <= center:
				return 1-(center - x)/(center-leftRad)
			return 1-(x - center)/(rightRad-center)

		intervals = []

		# as in an inverted v, like a caret
		for i in range(numPoints):
			x = leftRad + (i+random()*.1-.05) * xDelta
			nextX = leftRad + (i+random()*.1-.05+1) * xDelta
			y = plot(x)
			nextY = plot(nextX)
			intervals.append([[x,y],[nextX, nextY]])
		newDist = Distribution(intervals)
		newDist.normalize()
		return newDist

	@staticmethod
	def fromString(tableString, precision):
		# same as Distribution.fromTable, but takes a string, like what would be copy-pasted
		# out of excel. Specifically, it is tab separated columns and \n separated rows
		table = [x.strip() for x in tableString.split('\n') if x.strip() != '']

		table = list(map(lambda row: row.split('\t'), table))

		# convert rows 2+ to floats
		for i in range(1, len(table)):
			table[i] = [float(x) for x in table[i]]

		return Distribution.fromTable(table, precision)

	@staticmethod
	def fromTable(table, precision):
		# takes a 2d table like [['m', 'm error', 'v', 'v error'], [1.2,.3,3.2,.1]...]
		# returns a dictionary, with header values as keys and lists of distributions as key (the column)
		# note that an error column must be "var_name error". Example:
		#	m		m error		v		v error
		#	1.2		.3			3.2		.1
		#	1.3		.3			2.1		.2
		# turns into {'m': [dist_1.2, dist_1.3], 'v': [dist_3.2, dist_2.1]}
		headers = list(zip(table[0], range(len(table[0]))))
		values = list(filter(lambda x: not (' error' in x[0]), headers))
		errors = list(filter(lambda x: ' error' in x[0], headers))
		# match values to their error
		valueErrorPairs = []
		for val in values:
			valName = val[0]
			possibleErrors = list(filter(lambda x: x[0] == valName + ' error', errors))
			if len(possibleErrors) != 1:
				raise Exception('you appear to have duplicate columns named "'+valName+' error"')

			err = possibleErrors[0]

			# value, measure column, error column
			valueErrorPairs.append((valName, val[1], err[1]))

		# now create the dictionary
		rows = []

		# now iterate through the rows, making the distributions as we go
		for i in range(1, len(table)):
			rowDict = {}
			row = table[i]
			for pair in valueErrorPairs:
				val = row[pair[1]]
				err = row[pair[2]]
				rowDict[pair[0]] = Distribution.normalDistribution(val, err, precision)

			rows.append(rowDict)

		return rows

	@staticmethod
	def evaluateDataSet(expr, distributionTableOrig, resultFunc):
		# copy distributionTableOrig so that it won't be mutated
		distributionTable = []
		for d in distributionTableOrig:
			newDict = {}
			for key in d:
				newDict[key] = d[key].copy()
			distributionTable.append(newDict)
		trials = len(distributionTable)
		results = []
		for row in distributionTable:
			result = Distribution.evaluateExpression(expr, row)
			results.append(resultFunc(result))
		return results

	@staticmethod
	def evaluateExpression(expr, distDictionary):
		# evaluates an expression, like 'a**2 + b**2 - math.sin(c*a*b)', substituting
		# keys of distDictionary (so 'a', 'b', and 'c') with the corresponding distribution values
		# there cannot be a variable named reserved_key_name0 or 2 or 3...

		rootNode = ast.parse(expr, mode='exec')
		# starting at the root node, break the expression up recursively until every node (subset of expression)
		# contains at most 2 variables

		def getVariables(node):
			variables = []
			children = ast.walk(node)
			for child in children:
				gen = ast.iter_fields(child)
				for val in gen:
					if val[0] == 'id' and val[1] in distDictionary and val[1] not in variables:
						variables.append(val[1])
			return variables

		distributionVariableCount = 0

		def astor_source(node):
			res = astor.to_source(node, indent_with='').strip().replace('\n','')
			while ' ' in res:
				res = res.replace(' ','')
			return res

		def evaluate(node, variableName):
			variables = getVariables(node)
			source = astor_source(node)
			if len(variables) == 1:
				# save some time on simple evaluations:
				if variables[0] == source:
					return distDictionary[variables[0]]
				# do a singleArgCompute on the single distribution
				distribution = distDictionary[variables[0]]
				funcString = 'f = lambda '+variables[0]+':'+source
				exec(funcString)
				return distribution.singleArgCompute(locals()['f'])
			elif len(variables) == 0:
				raise Error('You appear to have no variable in this expression. How am I supposed to make a distribution out of it?')
			elif len(variables) == 2:
				funcString = 'f = lambda '+variables[0]+','+variables[1]+':'+source
				exec(funcString)
				return Distribution.computeDistribution(distDictionary[variables[0]], distDictionary[variables[1]], locals()['f'])
			else:
				# return the evaluation of the daughter nodes
				# replace the nodes with
				daughterNodeGenerator = ast.iter_child_nodes(node)
				varNames = []
				newSource = source
				daughters = []
				for daughter in daughterNodeGenerator:
					daughters.append(daughter)
				if len(daughters) == 1:
					return evaluate(daughters[0], variableName)
				if len(daughters) > 2:
					daughters = [daughters[0], daughters[-1]]

				# there are several possibilities:
				#  1. if any daughter has no variables, then do a single arg compute with the other daughter
				#  2. if no daughter has no variables, then do a full computeDistribution on the evaluation of the daughter

				# sort daughters by variables
				d0Vars = getVariables(daughters[0])
				d1Vars = getVariables(daughters[1])

				switched = False
				if len(d0Vars) > len(d1Vars):
					daughters = [daughters[1], daughters[0]]
					temp = d1Vars
					d1Vars = d0Vars
					d0Vars = temp
					switched = True

				if len(d0Vars) == 0:
					# replace the 2nd daughter with a new variable
					replacementDist = evaluate(daughters[1], variableName+'0')
					daughterSource = astor_source(daughters[1])
					newVarName = variableName+'_0'
					newSource = source.replace(daughterSource, newVarName)
					distDictionary[newVarName] = replacementDist

					funcString = 'f = lambda '+newVarName+':'+newSource
					exec(funcString)
					return replacementDist.singleArgCompute(locals()['f'])

				# replace both daughters with new variables
				rDist0 = evaluate(daughters[0], variableName+'1')
				rDist1 = evaluate(daughters[1], variableName+'2')
				dSource0 = astor_source(daughters[0])
				dSource1 = astor_source(daughters[1])
				vName0 = variableName+'_1'
				vName1 = variableName+'_2'
				# we need to replace in the correct order, from left to right
				newSource = source
				toReplace = ((dSource0, vName0), (dSource1, vName1))
				if switched:
					toReplace = ((dSource1, vName1), (dSource0, vName0))
				for r in toReplace:
					# consider r[0] = 'c' in the string c + tc. Make sure there are no letters around the match
					newSource = re.sub(r'(^|[^a-zA-Z_0-9])'+re.escape(str(r[0]))+r'($|[^a-zA-Z_0-9])', r'\1'+str(r[1])+r'\2', newSource, 1)
				distDictionary[vName0] = rDist0
				distDictionary[vName1] = rDist1

				funcString = 'f = lambda '+vName0+','+vName1+':'+newSource
				exec(funcString)
				return Distribution.computeDistribution(rDist0, rDist1, locals()['f'])

		return evaluate(rootNode, 'reserved_key_name')

	@staticmethod
	def computeDistribution(distA, distB, func):
		globalIntervals = []
		globalIntervalLengths = []
		minX = None
		maxX = None
		for i1 in distA.intervals:
			for i2 in distB.intervals:
				# we can't assume that the first coordinate cooresponds to the first coordinate of the second point 7292752669
				# we will find the widest interval (of all combos)
				xs = None
				try:
					xs = (func(i1[0][0], i2[0][0]), func(i1[0][0], i2[1][0]), func(i1[1][0], i2[0][0]), func(i1[1][0], i2[1][0]))
				except:
					continue

				ys = (i1[0][1]*i2[0][1], i1[0][1]*i2[1][1], i1[1][1]*i2[0][1], i1[1][1]*i2[1][1])

				lowestIndex = None
				lowest = None
				highestIndex = None
				highest = None
				for i in range(len(xs)):
					if lowestIndex == None or xs[i] < lowest:
						lowestIndex = i
						lowest = xs[i]
					if highestIndex == None or xs[i] > highest:
						highestIndex = i
						highest = xs[i]

				interval = [[xs[lowestIndex], ys[lowestIndex]], [xs[highestIndex], ys[highestIndex]]]

				for val in [interval[0][0], interval[1][0]]:
					if minX == None or val < minX:
						minX = val
					if maxX == None or val > maxX:
						maxX = val

				globalIntervals.append(interval)
				globalIntervalLengths.append((i1[1][0]-i1[0][0]) * (i2[1][0]-i2[0][0]))

		numPoints = max(distA.numPoints, distB.numPoints)

		condensed = Distribution.condenseIntervals(globalIntervals, globalIntervalLengths, minX, maxX, numPoints)

		intermediateDist = Distribution(condensed)
		intermediateDist.normalize()

		newXLow = intermediateDist.integrateTo(.0001)
		newXHigh = intermediateDist.integrateTo(.9999)

		# if a distribution is too thin (e.g. 10e-12 kinda thin), then it's possible that integrateTo(.0001) > integrateTo(.9999),
		# in which case I'm going to give up and just make a normal distribution centered around where we are looking at
		newDist = None
		if newXLow + 1e-10> newXHigh:
			newDist = Distribution.normalDistribution(globalIntervals[len(globalIntervals) // 2][0][0], abs(globalIntervals[0][0][0] - globalIntervals[-1][0][0])/5, numPoints)
		else:
			condensed = Distribution.condenseIntervals(globalIntervals, globalIntervalLengths, newXLow, newXHigh, numPoints)

			newDist = Distribution(condensed)
			newDist.normalize()

			newDist.smooth(3)

		return newDist

	# for functions like e^x
	def singleArgCompute(self, func):
		# map line segments, condense as before
		newSegs = []
		newSegLengths = []
		xMin = None
		xMax = None
		for seg in self.intervals:
			try:
				newSeg = [[func(seg[0][0]), seg[0][1]], [func(seg[1][0]), seg[1][1]]]
			except:
				continue
			if newSeg[0][0] > newSeg[1][0]:
				temp = newSeg[0]
				newSeg[0] = newSeg[1]
				newSeg[1] = temp

			for x in [newSeg[0][0], newSeg[1][0]]:
				if xMin == None or x < xMin:
					xMin = x
				if xMax == None or x > xMax:
					xMax = x

			newSegs.append(newSeg)
			newSegLengths.append(seg[1][0] - seg[0][0])

		condensed = Distribution.condenseIntervals(newSegs, newSegLengths, xMin, xMax, self.numPoints)

		intermediateDist = Distribution(condensed)
		intermediateDist.normalize()

		newXLow = intermediateDist.integrateTo(.0001)
		newXHigh = intermediateDist.integrateTo(.9999)

		condensed = Distribution.condenseIntervals(newSegs, newSegLengths, newXLow, newXHigh, self.numPoints)

		newDist = Distribution(condensed)
		newDist.normalize()

		newDist.smooth(3)

		return newDist

	@staticmethod
	def condenseIntervals(globalIntervals, intervalLengths, minX, maxX, numPoints):
		# condense into smaller set of intervals
		finalIntervals = []
		xDelta = (maxX - minX) / (numPoints)
		x = minX
		while x < maxX-xDelta/100:
			finalIntervals.append([[x,0],[x + xDelta, 0]])
			x += xDelta
		def intervalIndex(xVal):
			return math.floor((xVal - minX) / xDelta)

		def intervalOverlap(i1, i2):
			left = None
			right = None
			if i1[0][0] < i2[0][0]:
				left = i2[0][0]
			else:
				left = i1[0][0]
			if i1[1][0] < i2[1][0]:
				right = i1[1][0]
			else:
				right = i2[1][0]
			return right - left if right > left else 0

		# due to the fact that delta x is constant for distA and distB, the weight of each line segment is the same
		# we will arbitrarily designate this weight 1
		count = -1
		for i in globalIntervals:
			count += 1
			# find range of intervals to which this interval extends
			weightPerUnit = 1 / (i[1][0] - i[0][0])
			intList = finalIntervals[intervalIndex(i[0][0]) : intervalIndex(i[1][0]) + 1]

			for iL in intList:
				yDelta = intervalOverlap(iL, i) * .5 * (i[0][1] + i[1][1]) * weightPerUnit * intervalLengths[count]
				iL[0][1] += yDelta
				iL[1][1] += yDelta

		# zip up the interval endpoints (e.g. make them the same as their neighbors)
		i = 1
		while i < len(finalIntervals):
			aveY = .5*(finalIntervals[i-1][1][1] + finalIntervals[i][0][1])
			finalIntervals[i-1][1][1] = aveY
			finalIntervals[i][0][1] = aveY
			i += 1

		return finalIntervals

	def normalize(self):
		# makes the area under the curve 1 again
		cumulativeArea = self.integrate(self.intervals[0][0][0], self.intervals[-1][1][0])
		for i in self.intervals:
			i[0][1] /= cumulativeArea
			i[1][1] /= cumulativeArea

	def draw(self):
		lc = mc.LineCollection(self.intervals, linewidths=1)
		fig, ax = plt.subplots()
		ax.add_collection(lc)
		ax.autoscale()
		ax.margins(.1)

	def integrate(self, min, max):
		return self.integrateLeft(max) - self.integrateLeft(min)

	def smooth(self, kernelRadius):
		newIntervals = []
		kernel = [1/(2*kernelRadius-1)] * (2*kernelRadius-1)
		for i in range(len(self.intervals)):
			count = 0
			sum = 0
			for j in range(-kernelRadius, kernelRadius + 1):
				if 0 <= i + j < len(self.intervals):
					sum += self.intervals[i + j][1][1]
					count += 1
			if count > 0:
				newIntervals.append([[self.intervals[i][0][0], self.intervals[i][0][1]],[self.intervals[i][1][0], sum/count]])

		# reset first coordinate
		for i in range(1, len(newIntervals)):
			newIntervals[i][0][1] = newIntervals[i-1][1][1]

		self.intervals = newIntervals

	@staticmethod
	def intervalArea(interval, prevI, nextI):
		if prevI == None or nextI == None:
			# default to a trapezoid approximation
			return .5 * (interval[0][1] + interval[1][1]) * (interval[1][0] - interval[0][0])

		prevSlope = (prevI[1][1] - prevI[0][1]) / (prevI[1][0] - prevI[0][0])
		nextSlope = (nextI[1][1] - nextI[0][1]) / (nextI[1][0] - nextI[0][0])

		prevX = (prevI[0][0] + prevI[1][0]) / 2
		nextX = (nextI[0][0] + nextI[1][0]) / 2

		(px, py) = (prevX, (prevI[0][1] + prevI[1][1]) / 2)
		a = (nextSlope - prevSlope) / (2*nextX - 2*prevX)
		b = prevSlope - 2*a*prevX
		c = py - a*px**2 - b*px

		# now integrate from interval[0][0] to interval[1][0]
		# int(ax^2 + bx+c) = (1/3 ax^3 + 1/2 bx^2 + cx)
		return (1/3) * a*interval[1][0]**3 + .5*b*interval[1][0]**2 + c*interval[1][0] - (1/3 * a*interval[0][0]**3+.5*b*interval[0][0]**2+c*interval[0][0])

	def integrateLeft(self, max):
		area = 0
		cInterval = None
		prevInterval = None
		c = 0
		for interval in self.intervals:
			currentInterval = interval
			area += Distribution.intervalArea(currentInterval, prevInterval, self.intervals[c+1] if c + 1 < len(self.intervals) else None)
			if currentInterval[1][0] > max:
				break
			prevInterval = interval
			c += 1

		if currentInterval[1][0] > max:
			area -= (currentInterval[1][0] - max) * (currentInterval[0][1] + currentInterval[1][1]) * .5

		return area

	def integrateTo(self, areaFromLeft):
		# integrate from the side till we hit areaFromLeft
		total = 0
		currentInterval = None
		c = 0
		prevInterval = None
		for interval in self.intervals:
			currentInterval = interval
			total += Distribution.intervalArea(interval, prevInterval, self.intervals[c+1] if c+1 < len(self.intervals) else None)
			if total > areaFromLeft:
				break
			prevInterval = interval
			c += 1

		if total > areaFromLeft:
			return currentInterval[1][0] - (total - areaFromLeft) / (.5 * (currentInterval[0][1] + currentInterval[1][1]))

		return currentInterval[1][0]

	def getErrorBounds(self, standardDev):
		area = norm.cdf(standardDev) - .5

		maxs = self.getLocalMaxima()
		if len(maxs) > 0:
			# get the highest maximum
			maxs.sort(key = lambda x: x[1])
			max = maxs[-1]

			middle = max[0]
			middleArea = self.integrate(self.intervals[0][0][0], middle)
			lower = self.integrateTo(middleArea - area)
			upper = self.integrateTo(middleArea + area)
			return (lower, middle, upper)
		else:
			raise RuntimeError("The resulting distribution appears to have no local maxima")

	def getLocalMaxima(self):
		maxs = []
		# use rolle's theorem to find intervals which contain a local maxima, then approximate the location
		# using a quadratic regression of the slopes of the before and after intervals
		for i in range(1, len(self.intervals)-1):
			prevI = self.intervals[i-1]
			int = self.intervals[i]
			nextI = self.intervals[i+1]

			prevSlope = (prevI[1][1] - prevI[0][1]) / (prevI[1][0] - prevI[0][0])
			nextSlope = (nextI[1][1] - nextI[0][1]) / (nextI[1][0] - nextI[0][0])

			if prevSlope > 0 and nextSlope < 0:
				# we have slope = d/dx ax^2 + bx + c = 2ax + b
				prevX = (prevI[0][0] + prevI[1][0]) / 2
				nextX = (nextI[0][0] + nextI[1][0]) / 2

				centerX = prevX + (nextX - prevX) * prevSlope / (prevSlope-nextSlope)
				# we have s = 2ax + b, prevSlope = 2a(prevX) + b, nextSlope = 2a(nextX) + b =>
				# prevSlope - 2a(prevX) = nextSlope - 2a(nextX), a = (nextSlope - prevSlope) / (2 nextX - 2 prevX)
				# b = prevSlope - 2a(prevX)
				# given a point (px, py)...
				# a(px)^2+b(px)+c = py, c = py - a(px)^2 - b(px)
				(px, py) = (prevX, (prevI[0][1] + prevI[1][1]) / 2)
				a = (nextSlope - prevSlope) / (2*nextX - 2*prevX)
				b = prevSlope - 2*a*prevX
				c = py - a*px**2 - b*px

				centerY = a*centerX**2 + b*centerX + c
				maxs.append([centerX, centerY])

		return maxs

	def copy(self):
		intervalCopy = []
		for i in self.intervals:
			intervalCopy.append([[i[0][0], i[0][1]], [i[1][0], i[1][1]]])
		newSelf = Distribution(intervalCopy)
		newSelf.numPoints = self.numPoints
		return newSelf

	def __mul__(self, other):
		if isinstance(other,self.__class__):
			return Distribution.computeDistribution(self, other, lambda x,y: x*y)
		return self.singleArgCompute(lambda x: x*other)
	def __rmul__(self, other):
		return self*other
	def __add__(self, other):
		if isinstance(other, self.__class__):
			return Distribution.computeDistribution(self, other, lambda x,y: x+y)
		return self.singleArgCompute(lambda x: x+other)
	def __radd__(self, other):
		return self + other
	def __sub__(self, other):
		if isinstance(other,self.__class__):
			return Distribution.computeDistribution(self, other, lambda x,y: x-y)
		return self.singleArgCompute(lambda x: x-other)
	def __rsub__(self, other):
		return -self + other
	def __truediv__(self, other):
		if isinstance(other, self.__class__):
			return Distribution.computeDistribution(self, other, lambda x,y: x/y)
		return self.singleArgCompute(lambda x: x/other)
	def __rtruediv__(self, other):
		return self**(-1) * other
	def __pow__(self, other):
		if isinstance(other, self.__class__):
			return Distribution.computeDistribution(self, other, lambda x,y: x**y)
		return self.singleArgCompute(lambda x: x**other)
	def __rpow__(self, other):
		return self**other
	def __neg__(self):
		return self.singleArgCompute(lambda x:-x)

	# same as distZip(dists, copy(dists), func)
	@staticmethod
	def mapSelf(dists, func):
		for i in range(len(dists)):
			dists[i] = Distribution.computeDistribution(dists[i], dists[i].copy(), func)

	@staticmethod
	def foldIntoFirst(dists, func):
		result = dists[0]
		i = 1
		while i < len(dists):
			result = Distribution.computeDistribution(result, dists[i], func)
			i += 1

		return result

	def saveToFile(self, fileName):
		file = open(fileName, 'w')
		file.write(str(self.intervals))
		file.close()

	@staticmethod
	def fromFile(fileName):
		file = open(fileName, 'r')
		return Distribution(eval(file.read()))

	# pts is in the form [[p1x, p1y], [p2x, p2y], ...]
	@staticmethod
	def linearRegression(pts):
		# calculates a distribution representing the slope and y intercept of the calculated linear regression from the points
		# slope = sum(1-n, (xi-x_bar)(yi-y_bar))/sum(1-n,(xi-x_bar)^2)
		xs = [x[0] for x in pts]
		ys = [x[1] for x in pts]
		xBar = Distribution.foldIntoFirst(xs, lambda a,b: a+b)
		yBar = Distribution.foldIntoFirst(ys, lambda a,b: a+b)
		xBar = xBar.singleArgCompute(lambda x: x/len(pts))
		yBar = yBar.singleArgCompute(lambda x: x/len(pts))

		numerators = [Distribution.evaluateExpression('(xi-xbar)*(yi-ybar)', {'xi':j[0], 'xbar':xBar, 'ybar':yBar, 'yi':j[1]}) for j in pts]
		numeratorTotal = Distribution.foldIntoFirst(numerators, lambda a,b: a+b)

		denominators = [Distribution.evaluateExpression('(xi-xbar)**2',{'xi':j[0],'xbar':xBar}) for j in pts]
		denominatorTotal = Distribution.foldIntoFirst(denominators, lambda a,b: a+b)

		slope = Distribution.computeDistribution(numeratorTotal, denominatorTotal, lambda a,b: a/b)

		# y int = y_bar - slope*x_bar
		yInt = Distribution.evaluateExpression('ybar - slope*xbar',{'ybar':yBar,'slope':slope,'xbar':xBar})

		return (slope, yInt)

	# combines two lists of distributions elementwise according to func
	@staticmethod
	def distZip(dists1, dists2, func):
		if len(dists1) != len(dists2):
			raise Error('dist1 and dist2 must be the same length')

		res = []
		for i in range(len(dists1)):
			res.append(Distribution.computeDistribution(dists1[i], dists2[i], func))

		return res

	def show(self):
		plt.show()

	def toFastDist(self):
		bounds = self.getErrorBounds(1)
		return FastDistribution(bounds[1], (bounds[2] - bounds[0])/2)

	# converts all the ordinary Distributions to FastDistributions
	@staticmethod
	def toFastDists(dists):
		newArr = []
		# dists is an array of dictionaries
		for el in dists:
			newArr.append({})
			# gets keys from el
			for key in [*el]:
				errorInterval = el[key].getErrorBounds(1)
				newArr[-1][key] = FastDistribution(errorInterval[1], (errorInterval[1] - errorInterval[0] + errorInterval[2] - errorInterval[1])/2)
		return newArr

# FastDistribution uses the normal derivative propagation method:
# 	delta f(x,y) = sqrt((d/dx f(x,y) delta x)^2 + (d/dy f(x,y)^2 delta y)^2
class FastDistribution:
	# model each distribution as a normal distribution
	def __init__(self, center, stdv):
		self.center = center
		self.stdev = stdv

	def toDistribution(self, numPoints):
		return Distribution.normalDistribution(self.center, self.stdev, numPoints)

	# we want the FastDistribution and normal Distribution objects to share many of the same
	# functions (where applicable)
	def singleArgCompute(self, func):
		stdDev = self.stdev * FastDistribution.D(func, self.center)
		# the center is shifted by the function
		center = func(self.center)
		return FastDistribution(center, stdDev)

	@staticmethod
	def computeDistribution(d1, d2, func):
		newCenter = func(d1.center, d2.center)
		newError = ((d1.stdev*(FastDistribution.D(lambda x:func(x, d2.center),d1.center)))**2 + (d2.stdev*(FastDistribution.D(lambda x:func(d1.center, x),d2.center)))**2)**(1/2)
		return FastDistribution(newCenter, newError)

	# numeric derivative of func with respect to its (one) variable evaluated at pt
	@staticmethod
	def D(func, pt):
		# a symmetric derivative for kicks
		# this is meant to be fast, so an arbitrary delta will probably be sufficient
		delta = abs(pt / 1e8)
		if pt == 0:
			delta = 10**(-8)
		right = func(pt + delta)
		left = func(pt - delta)
		return (right - left) / (2*delta)

	def copy(self):
		return FastDistribution(self.center, self.stdev)

	def getErrorBounds(self, stdDev):
		return (self.center - self.stdev * stdDev, self.center, self.center + self.stdev * stdDev)

	# copied from Distribution
	@staticmethod
	def foldIntoFirst(dists, func):
		result = dists[0]
		i = 1
		while i < len(dists):
			result = FastDistribution.computeDistribution(result, dists[i], func)
			i += 1

		return result

	def saveToFile(self, fileName):
		file = open(fileName, 'w')
		file.write(str(self.center)+','+str(self.stdev))
		file.close()

	@staticmethod
	def fromFile(name):
		file = open(fileName, 'r')
		res = file.read().split(',')
		return FastDistribution(float(res[0]), float(res[1]))

	# mostly copied, with changes to reflect that I don't want to copy evaluateExpression over here
	def linearRegression(pts):
		# calculates a distribution representing the slope and y intercept of the calculated linear regression from the points
		# slope = sum(1-n, (xi-x_bar)(yi-y_bar))/sum(1-n,(xi-x_bar)^2)
		xs = [x[0] for x in pts]
		ys = [x[1] for x in pts]
		xBar = FastDistribution.foldIntoFirst(xs, lambda a,b: a+b)
		yBar = FastDistribution.foldIntoFirst(ys, lambda a,b: a+b)
		xBar = xBar.singleArgCompute(lambda x: x/len(pts))
		yBar = yBar.singleArgCompute(lambda x: x/len(pts))

		numerators = [(j[0]-xBar)*(j[1]-yBar) for j in pts]
		numeratorTotal = FastDistribution.foldIntoFirst(numerators, lambda a,b: a+b)

		denominators = [(j[0]-xBar)**2 for j in pts]
		denominatorTotal = FastDistribution.foldIntoFirst(denominators, lambda a,b: a+b)

		slope = FastDistribution.computeDistribution(numeratorTotal, denominatorTotal, lambda a,b: a/b)

		# y int = y_bar - slope*x_bar
		yInt = yBar - slope*xBar

		return (slope, yInt)

	# this is copied from distribution
	def __mul__(self, other):
		if isinstance(other,self.__class__):
			return FastDistribution.computeDistribution(self, other, lambda x,y: x*y)
		return self.singleArgCompute(lambda x: x*other)
	def __rmul__(self, other):
		return self*other
	def __add__(self, other):
		if isinstance(other, self.__class__):
			return FastDistribution.computeDistribution(self, other, lambda x,y: x+y)
		return self.singleArgCompute(lambda x: x+other)
	def __radd__(self, other):
		return self + other
	def __sub__(self, other):
		if isinstance(other,self.__class__):
			return FastDistribution.computeDistribution(self, other, lambda x,y: x-y)
		return self.singleArgCompute(lambda x: x-other)
	def __rsub__(self, other):
		return -self + other
	def __truediv__(self, other):
		if isinstance(other, self.__class__):
			return FastDistribution.computeDistribution(self, other, lambda x,y: x/y)
		return self.singleArgCompute(lambda x: x/other)
	def __rtruediv__(self, other):
		return self**(-1) * other
	def __pow__(self, other):
		if isinstance(other, self.__class__):
			return FastDistribution.computeDistribution(self, other, lambda x,y: x**y)
		return self.singleArgCompute(lambda x: x**other)
	def __rpow__(self, other):
		return self**other
	def __neg__(self):
		return self.singleArgCompute(lambda x:-x)
