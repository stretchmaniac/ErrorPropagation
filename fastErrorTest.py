from errorProp import Distribution, FastDistribution

raw = '''
a	b	a error	b error
1	8	0.1	0.1
2	5	0.1	0.1
3	8	0.1	0.1
4	5	0.1	0.1
4	3	0.1	0.1
'''

dists = Distribution.fromString(raw, 150)
fastDists = Distribution.toFastDists(dists)

(d1, d2) = (fastDists[2]['a'], fastDists[2]['b'])
(realD1, realD2) = (dists[2]['a'], dists[2]['b'])

print('-------getErrorBounds Initial Check---------')
print(d1.getErrorBounds(1), d2.getErrorBounds(1))
print(realD1.getErrorBounds(1), realD2.getErrorBounds(1))

print('-------Addition Check----------')
a1 = d1 + d2
realA1 = realD1 + realD2
print(a1.getErrorBounds(1))
print(realA1.getErrorBounds(1))
print('-------Multiplication Check---------')
m1 = d1*d2
realM1 = realD1*realD2
print(m1.getErrorBounds(1))
print(realM1.getErrorBounds(1))
print('-------Negation Check----------')
n1 = -d1
realN1 = -realD1
print(n1.getErrorBounds(1))
print(realN1.getErrorBounds(1))

pts = [[k['a'], k['b']] for k in fastDists]
realPts = [[k['a'], k['b']] for k in dists]
print('--------Linear Regression Check----------')
print([x.getErrorBounds(1) for x in [*FastDistribution.linearRegression(pts)]])
print([x.getErrorBounds(1) for x in [*Distribution.linearRegression(realPts)]])
