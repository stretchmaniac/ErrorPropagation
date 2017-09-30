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

theoreticalRangeEquation = '(v*math.cos(theta) / 9.8) * (v*math.sin(theta) + math.sqrt(v**2*math.sin(theta)**2 + 2*9.8*(h_r-h_c)))'
measuredRangeEquation = 'math.sqrt((hyp1+hyp2+hyp3)**2/9 - (h_r-h_c)**2)'
measuredTheta = 'math.asin((bar_front - bar_back)/bar_len)'

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
trials = 10
vDistArray = []
for i in range(10):
	vDistArray.append(vDist.copy())
measuredDistributions['v'] = vDistArray

results = Distribution.evaluateDataSet(measuredRangeEquation, measuredDistributions, lambda result: result.getErrorBounds(1))

for r in results:
	print(r)
