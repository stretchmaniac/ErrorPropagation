from errorProp import Distribution
import os.path
import math

# error for velocity, the average of the sum of the quotient of 6 different measured values

vDist = None

if os.path.exists('Physics 1401V Lab 1.3 data/v.txt'):
	vDist = Distribution.fromFile('Physics 1401V Lab 1.3 data/v.txt')
else:
	vEquation = '(d/t1 + d/t2 + d/t3) / 3'

	measuredData ='''
		d	d error	t1	t1 error	t2	t2 error	t3	t3 error
		0.0245	0.0012	0.008518	0.00001	0.008256	0.00001	0.008256	0.00001
		0.0245	0.0012	0.00541	0.00001	0.00548	0.00001	0.005458	0.00001
		0.0245	0.0012	0.00403	0.00001	0.00404	0.00001	0.003969	0.00001
	'''

	rows = Distribution.fromString(measuredData, 100)
	results = Distribution.evaluateDataSet(vEquation, rows, lambda r: r)

	vDist = results[0]

	vDist.saveToFile('Physics 1401V Lab 1.3 data/v.txt')
	results[1].saveToFile('Physics 1401V Lab 1.3 data/v2.txt')
	results[2].saveToFile('Physics 1401V Lab 1.3 data/v3.txt')

theoreticalRangeEquation = '(v*math.cos(theta) / 9.8) * (v*math.sin(theta) + math.sqrt(v**2*math.sin(theta)**2 + 2*9.8*(h_r-h_c)))'
measuredRangeEquation = 'math.sqrt((hyp1+hyp2+hyp3)**2/9 - (h_r-h_c)**2)'
measuredThetaEquation = 'math.asin((bar_front - bar_back)/bar_len)'

measuredData = '''
	hyp1	hyp1 error	hyp2	hyp2 error	hyp3	hyp3 error	h_r	h_r error	h_c	h_c error	bar_front	bar_front error	bar_back	bar_back error	bar_len	bar_len error
	1.711	0.005	1.712	0.005	1.711	0.005	1.029	0.005	0.015	0.001	0	0.1	0	0.1	20.79	0.1
	1.827	0.005	1.849	0.005	1.872	0.005	1.031	0.005	0.015	0.001	22.8	0.1	19.52	0.1	20.79	0.1
	1.902	0.005	1.892	0.005	1.918	0.005	1.033	0.005	0.015	0.001	23.3	0.1	16.47	0.1	20.79	0.1
	1.963	0.005	1.939	0.005	1.958	0.005	1.039	0.005	0.015	0.001	23.8	0.1	13.505	0.1	20.79	0.1
	1.858	0.005	1.914	0.005	1.836	0.005	1.048	0.005	0.015	0.001	24.24	0.1	11	0.1	20.79	0.1
	1.752	0.005	1.73	0.005	1.746	0.005	1.048	0.005	0.015	0.001	24.8	0.1	9.11	0.1	20.79	0.1
	1.624	0.005	1.63	0.005	1.59	0.005	1.036	0.005	0.015	0.001	25.3	0.1	7.94	0.1	20.79	0.1
	1.396	0.005	1.404	0.005	1.378	0.005	1.036	0.005	0.015	0.001	25.9	0.1	6.87	0.1	20.79	0.1
	1.114	0.005	1.108	0.005	1.114	0.005	1.036	0.005	0.015	0.001	26.6	0.1	6.4	0.1	20.79	0.1
	1.036	0.005	1.048	0.005	1.04	0.005	1.04	0.005	0.015	0.001	27	0.1	6.4	0.1	20.79	0.1
'''

measuredDistributions = Distribution.fromString(measuredData, 100)
for row in measuredDistributions:
	row['v'] = vDist.copy()

measuredRangeResults = Distribution.evaluateDataSet(measuredRangeEquation, measuredDistributions, lambda r: r)
print('completed measured range results')
measuredThetaResults = Distribution.evaluateDataSet(measuredThetaEquation, measuredDistributions, lambda r: r)
print('completed measured theta results')

# save all of our work
for i in range(len(measuredRangeResults)):
	print('saving point '+str(i)+' out of '+str(len(measuredRangeResults)-1))
	thetaDist = measuredThetaResults[i]
	rangeDist = measuredRangeResults[i]
	# we have the data point (theta, range)
	thetaDist.saveToFile('Physics 1401V Lab 1.3 data/mT+'+str(i)+'.txt')
	rangeDist.saveToFile('Physics 1401V Lab 1.3 data/mR+'+str(i)+'.txt')

thetaValues = []
for res in measuredThetaResults:
	thetaValues.append(res.getErrorBounds(1))

# also, handily enough, works for things besides theta
def getThetaError(center, thetaValues):
	if center < thetaValues[0][1]:
		return (thetaValues[0][2] - thetaValues[0][0]) / 2
	if center > thetaValues[-1][1]:
		return (thetaValues[-1][2] - thetaValues[-1][0]) / 2

	# find sequential thetaValues such that the center falls inbetween the centers
	for i in range(0, len(thetaValues)):
		thetaBefore = thetaValues[i]
		thetaAfter = thetaValues[i+1]
		if thetaBefore[1] <= center <= thetaAfter[1]:
			# return the linear interpolation of the error
			errorBefore = (thetaBefore[2] - thetaBefore[0]) / 2
			errorAfter = (thetaAfter[2] - thetaAfter[0]) / 2
			return errorBefore + (center - thetaBefore[1]) / (thetaAfter[1] - thetaBefore[1]) * (errorAfter - errorBefore)

# because we want more theoretical error points than there are measured theta errors, we will linearly
# interpolate our previous theta error values to get a reasonable estimate
theoreticalResults = []
for c in range(101):
	theta = (c / 100) * math.pi / 2
	error = getThetaError(theta, thetaValues)
	row = {}
	row['theta'] = Distribution.normalDistribution(theta, error, 100)
	row['v'] = vDist.copy()
	row['h_c'] = Distribution.normalDistribution(.015, .001, 100)
	# h_r is empiracly roughly 1.029 + 0.2*sin(theta)
	row['h_r'] = Distribution.normalDistribution(1.029 + 0.2*math.sin(theta), .005, 100)

	theorDist = Distribution.evaluateExpression(theoreticalRangeEquation, row)
	theoreticalResults.append(theorDist)
	row['theta'].saveToFile('Physics 1401V Lab 1.3 data/tT+'+str(c)+'.txt')
	theorDist.saveToFile('Physics 1401V Lab 1.3 data/tR+'+str(c)+'.txt')
	print('completed theoretical data point '+str(c)+' out of 100: '+str(theoreticalResults[-1].getErrorBounds(1)))
