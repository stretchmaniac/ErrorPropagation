from errorProp import Distribution
import math

# so we don't have to recompute the distributions every time

click1Dist = Distribution.fromFile('Physics 1401V Lab 1.3 data/v.txt')
click2Dist = Distribution.fromFile('Physics 1401V Lab 1.3 data/v2.txt')
click3Dist = Distribution.fromFile('Physics 1401V Lab 1.3 data/v3.txt')

thetaDists =  []
rangeDists = []
# get the theta error
for i in range(10):
    thetaDist = Distribution.fromFile('Physics 1401V Lab 1.3 data/mT+'+str(i)+'.txt')
    rangeDist = Distribution.fromFile('Physics 1401V Lab 1.3 data/mR+'+str(i)+'.txt')
    thetaDists.append(thetaDist)
    rangeDists.append(rangeDist)

tThetaDists = []
tRangeDists = []
for i in range(101):
    tThetaDist = Distribution.fromFile('Physics 1401V Lab 1.3 data/tT+'+str(i)+'.txt')
    tRangeDist = Distribution.fromFile('Physics 1401V Lab 1.3 data/tR+'+str(c)+'.txt')
    tThetaDists.append(tThetaDist)
    tRangeDists.append(tRangeDist)

# we want to output this information to a mathematica graph with a shaded region for the theoretical
# and error bars for the measured points 
