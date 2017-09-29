from errorProp import Distribution
import os.path
import math

# error for velocity, the average of the sum of the quotient of 6 different measured values

vDist = None

if os.path.exists('Physics 1401V Lab 1.3 data/v.txt'):
	vDist = Distribution.fromFile('Physics 1401V Lab 1.3 data/v.txt')
else:
	d1 = Distribution.normalDistribution(.0245, .004, 500)
	d2 = d1.copy()
	d3 = d1.copy()

	t1 = Distribution.normalDistribution(.008518, .00001, 500)
	t2 = Distribution.normalDistribution(.00829, .00001, 500)
	t3 = Distribution.normalDistribution(.008256, .00001, 500)

	quotients = Distribution.distZip([d1, d2, d3], [t1, t2, t3], lambda x,y:x/y)
	sum = Distribution.foldIntoFirst(quotients, lambda x,y: x + y)

	final = sum.singleArgCompute(lambda x: x/3)

	final.draw()
	bounds = final.getErrorBounds(1)
	print(bounds)
	print(bounds[1] - bounds[0], bounds[2] - bounds[1])
			
	final.show()
	vDist = final
	final.saveToFile('Physics 1401V Lab 1.3 data/v.txt')

# now we get to data error measurements

# copy-pasted from excel doc
rawData = [
	#	h_back	h_front	bar_len	h_lnch	h_cork	h1		h2		h3	
	[	0,		0,		20.79,	1.029,	0.015,	1.711,	1.712,	1.711],
	[	19.52,	22.8,	20.79,	1.031,	0.015,	1.827,	1.849,	1.872],
	[	16.47,	23.3,	20.79,	1.033,	0.015,	1.902,	1.892,	1.918],
	[	13.505,	23.8,	20.79,	1.039,	0.015,	1.963,	1.939,	1.958],
	[	11,		24.24,	20.79,	1.048,	0.015,	1.858,	1.914,	1.836],
	[	9.11,	24.8,	20.79,	1.048,	0.015,	1.752,	1.73,	1.746],
	[	7.94,	25.3,	20.79,	1.036,	0.015,	1.624,	1.63,	1.59],
	[	6.87,	25.9,	20.79,	1.036,	0.015,	1.396,	1.404,	1.378],
	[	6.4,	26.6,	20.79,	1.036,	0.015,	1.114,	1.108,	1.114],
	[	6.4,	27,		20.79,	1.04,	0.015,	1.036,	1.048,	1.04]
]

data = []
for row in rawData:
	data.append({
		# name: [data, error]
		'h_back':[row[0], .1],
		'h_front':[row[1], .1],
		'bar_len':[row[2], .1],
		'h_launch':[row[3], .005],
		'h_cork':[row[4], .005],
		'hs':[[row[5], .005], [row[6], .005], [row[7], .005]],
	})

precision = 100
errors = []
for row in data:
	errorRow = {}
	# find error in calculated length from table, that is, L = sqrt(h_ave^2 - (h_launch - h_cork)^2)
	hs = list(map(lambda x: Distribution.normalDistribution(x[0], x[1], precision), row['hs']))
	# average of the hypotenuse
	h_ave = Distribution.foldIntoFirst(hs, lambda x,y: x+y).singleArgCompute(lambda x: x / len(hs))
	
	height = Distribution.normalDistribution(*row['h_launch'], precision)
	cork_height = Distribution.normalDistribution(*row['h_cork'], precision)
	
	intermediate = Distribution.computeDistribution(height, cork_height, lambda x,y: x - y)
	lengthFromTable = Distribution.computeDistribution(h_ave, intermediate, lambda x,y: math.sqrt(x**2 - y**2))
	
	errorRow['L'] = lengthFromTable
	
	# theta error, where theta = asin((h_front - h_back) / bar_len))
	h_front = Distribution.normalDistribution(*row['h_front'], precision)
	h_back = Distribution.normalDistribution(*row['h_back'], precision)
	bar_len = Distribution.normalDistribution(*row['bar_len'], precision)
	
	intermediate2 = Distribution.computeDistribution(h_front, h_back, lambda x,y: x - y)
	theta = Distribution.computeDistribution(intermediate2, bar_len, lambda i, b: math.asin(i / b))
	
	errorRow['theta'] = theta
	
	# theoretical error, given by range = (v cos(theta) / g) (v sin(theta) + sqrt(v^2 sin(theta)^2 + 2 g h))
	
	h_launch = Distribution.normalDistribution(*row['h_launch'], precision)
	h_cork = Distribution.normalDistribution(*row['h_cork'], precision)
	
	h_height = Distribution.computeDistribution(h_launch, h_cork, lambda x,y: x-y)
	
	g = 9.8
	theor1 = Distribution.computeDistribution(vDist, theta, lambda v,theta: v**2 * math.sin(theta)**2)
	theor2 = Distribution.computeDistribution(theor1, h_height, lambda th,h: math.sqrt(th + 2*g*h))
	theor3 = Distribution.computeDistribution(vDist, theta, lambda v, theta: v*math.sin(theta))
	theor4 = Distribution.computeDistribution(theor2, theor3, lambda x,y: x+y)
	theor5 = Distribution.computeDistribution(vDist, theta, lambda v, t: v*math.cos(t) / g)
	range = Distribution.computeDistribution(theor5, theor4, lambda x,y: x*y)
	
	errorRow['theoretical_range'] = range
	
	errors.append(errorRow)
	
for e in errors:
	measuredError = e['L'].getErrorBounds(1)
	expectedError = e['theoretical_range'].getErrorBounds(1)
	thetaError = e['theta'].getErrorBounds(1)
	def toExcel(error):
		return str(error[2] - error[1])+'\t'+str(error[1] - error[0])
	print(toExcel(measuredError)+'\t'+toExcel(expectedError)+'\t'+toExcel(thetaError))
	
	

	



