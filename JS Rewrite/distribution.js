// Alan Koval, started 11/17/2017

function Distribution(center, standardDev, numPoints, minZ, maxZ){
  center = center | 0;
  standardDev = standardDev | 1;
  numPoints = numPoints | 150;
  minZ = minZ | -7;
  maxZ = maxZ | 7;

  // we want model a distribution (a bell curve) from a set number of points. To be intelligent about this,
  // we want to have the maximum density of points where the curvature (2nd derivative) of the distribution
  // has the greatest magnitude.

  // to do this, we'll construct a normal distribution with center 0 and standard deviation 1, then scale it

  // we want to limit the change in first derivative to a set amount. Since the first derivative changes
  // direction at -1 and 1, we'll use these as piecewise waypoints, in a sense

  let firstDeriv = x => -x * math.E ** (-x**2/2) / math.sqrt(2*math.PI);
  // we have the integral of the absolute value of the 2nd derivative of the normal distribution
  // from negative infinity to infinity is 0.9678828981
  let totalDerivAllowance = 0.9678828981;
  // with a point at each endpoint, the number of intervals is numPoints - 1
  let stepDerivAllowance = totalDerivAllowance / (numPoints - 1);
  let xs = [];
  let currentX = minZ;

  while(xs.length <= numPoints){

  }
}

function newtonsMethod(func, seed){

}
