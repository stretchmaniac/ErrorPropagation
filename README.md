# ErrorPropagation

Suppose a physical quantity `A` has some measurement probability distribution `a(x)`, 
a function describing the probability that `A` was actually measured as `x`. 
Given another physical quantity `B` and measurement probability
distribution `b(x)`, there exists some distribution `(a+b)(x)` that gives the probability 
that the measurement of `A` plus the measurement of `B` is equal to `x`. 

This is known as arithmetic on random variables. Closed form solutions for normal distributions 
are possible in easy cases (such as the addition example above). More generally one must 
use approximation or numerical solutions. This library is one such numerical tool. I found it 
quite handy doing error analysis in my introductory physics course. 

Setup
===
ErrorPropagation depends on a couple common libraries. Make sure you have the following installed
- [numpy](http://www.numpy.org/)
- [scipy](https://www.scipy.org/)
- [matplotlib](https://matplotlib.org/)

Clone this repo inside your working directory

```
git clone https://github.com/stretchmaniac/ErrorPropagation.git
```

Now you can use ErrorPropagation as a module

```python
from ErrorPropagation.errorProp import Distribution 

x = Distribution.normalDistribution(1, 1, 300)
y = Distribution.normalDistribution(2, 1, 300)
summed = x + y

print(summed.getErrorBounds(1))

# output: (1.5829763108178432, 2.9997301924271893, 4.416104322244424)
```

Simple Usage
=== 
```python
from ErrorPropagation.errorProp import Distribution 
import math

# create a random variable that is normally distributed around 0 
#   with standard deviation 1. Model this distribution with 300 points
x = Distribution.normalDistribution(0, 1, 300)

# create a upside-down "v" shaped distribution from -5 to 1 with peak
#   at 0, again with 300 points
y = Distribution.vDistribution(-5, 0, 1, 300)

# maybe we already computed the sum of x and y
z = Distribution.fromFile('sumXY.txt')

# do some elementary arithmetic 
prodXY = x * y
twoPlusZ = 2 + z
zMinusX = z - x # not equal to y since they are random variables
quotientYX = y / x # recall x is centered at zero

# custom functions
# ... of one variable 
sinZ = z.singleArgCompute(lambda t: math.sin(t))
# ... of two variablesErrorBounds
minXY = Distribution.computeDistribution(x, y, lambda a,b: min(a,b))

# visualize some of our computations
Distribution.draw(y, zMinusX)
Distribution.draw(sinZ)
Distribution.draw(quotientYX)
Distribution.draw(minXY,x,y)

# get a 68% (+/- 1 std. deviation) confidence interval for the mean value for y / x
print(quotientYX.getErrorBounds(1))

# save what we want for later 
quotientYX.saveToFile('quotient.txt')
```

Function Reference
===
Initialization 
---
`constructor Distribution(intervals)`
- Initializes a distribution by a set of intervals. Each interval is a line, together which form the shape of the distribution. In particular, 
the left point of one interval must match the right point of the preceding interval. In general one does not use this constructor.

`static normalDistribution(center, standardDev, precision)`
- Creates a normally distributed random variable with center `center` and standard deviation `standardDev`. It is modelled with 
`precision` number of intervals.

`static vDistribution(left, center, right, precision)`
- Creates a random variable with triangle distribution with minimum `left`, maximum `right` and mean `center`. It is modelled 
with `precision` number of intervals. 

`static fromTable(table, precision)`
- takes an array of data named with headers and returns a dictionary of normally distributed variables, keyed by the 
name of the variable. Each header 
must be accompanied by a matching header with `error` appended to the name of the variable. These columns give the associated 
standard deviation of the variable. Each random variable is modelled with a distribution with `precision` number of intervals.
```python
Distribution.fromTable([
    ['a', 'a error', 'b', 'b error'],
    [1.0, 0.1, 10.0, 0.05], 
    [2.1, 0.12, 5.0, 0.05]
], p)
# returns the dictionary (nd(...) is shorthand for Distribution.normalDistribution(...))
# {
#   a: [nd(1.0, 0.1, p), nd(2.1, 0.12, p)],
#   b: [nd(10.0, 0.05, p), nd(5.0, 0.05, p)]  
# }
```

`static fromString(tableString, precision)`
- Facilitates integration with a spreadsheet program. Copy and paste a table created in MS Excel or Google Sheets 
in the format as `fromTable(table, precision)` directly
into `tableString`. Returns the same dictionary object as `fromTable(table, precision)`.

Computation
---
`static computeDistribution(var1, var2, func)`
- Calculates the resulting random variable `func(var1, var2)`, where `var1` and `var2` are random 
variables (Distributions) and `func` is a function of two variables. 

`static fromFile(file)`
- Retrieves a Distribution object from the file `file` located in the working directory. See `var1.saveToFile(file)`.

`var1.singleArgCompute(func)`
- Calculates the resulting random variable `func(var1)`, where `var1` is a random variable (Distribution)
and `func` is a function of one variable. 

`static evaluateExpression(expr, varDictionary)`
- Computes an expr (string) in python syntax given a dictionary of variables. For example:
```python
Distribution.evaluateExpression(
    'a + b + c', # any python syntax is allowed here
    {
        a: Distribution.normalDistribution(0, 1.0, 300),
        b: Distribution.normalDistribution(1, 0.5, 250),
        c: Distribution.normalDistribution(2, 1.2, 300)
    }
)
# returns a random variable with mean at 0 + 1 + 2 = 3
```

Visualization
---
`static draw(var1[, var2[,...]])`
- Plots the graphs of the probability distributions of `var1`, `var2`, ... on a single plot. Note that 
when calling `draw` more than once one will have to exit the plot window of one plot before the next 
will appear. 

Analysis
---
`var1.getErrorBounds(standardDev)`
- Gives a tuple `(a, m, b)` so that `m` is the mean of `var1` and `var1.cdf(a,b) = normalcdf(-standardDev, standardDev)`.   

`var1.integrate(min, max)`
- Calculates `var1.cdf(min, max)`.

`var1.getLocalMaxima()`
- Gives a best-guess for the local maxima in the probability distribution function of `var1`. This works better 
with calls to `var1.smooth(...)`.

Miscellaneous
--- 
`var1.smooth(kernelRadius)`
- Gives the random variable formed by the convolution of the square function `f(x) = 1` with width `kernelRadius * w` 
and the probability distribution function of `var1`, where `w` is the width of an individual interval making up 
`var1`. This generally has a smoothing effect. 

`var1.copy()`
- Creates a deep copy of the Distribution object `var1`. 

`var1.saveToFile(file)`
- Saves the contents of the Distribution object `var1` to the file `file`, which will be overwritten or 
created in the working directory. See `static fromFile(file)`. 