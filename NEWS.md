# multipointR 0.99.10
* clip window to zero if the window extends into negative coordinates. This
can happen if points are close to zero or coordinates are standardised to (0,1)
and `spatstat.geom::ripras` pads too much around it. 

# multipointR 0.99.9
* new default behaviour for window set in `speToPPP`. Before, this was a
rectanglular window given the min and max of the x,y coordiantes of the points.
From now on, the window is estimated by default as the convex hull surrounding 
all points with the function `spatstat.geom::ripras`. This makes the default
more sensible for e.g. TMAs. 
* default is to not standardise the x,y coordiantes to (0,1) anymore as this can
lead to negative edges with the `ripras` functionality and might also be
counter-intuitive to do this by default. 

# multipointR 0.99.8
* wrote a new function to perform fast approximations of the covariance matrix 
of a lgcp model. The functionality was present already in `vcov.kppm` but could
not be called in `summary.kppm`. Wrote a new generic `summary` for 
`multipointRppm` objects that allows to call the fast `vcov.kppm` function.

# multipointR 0.99.7
* complete rewrite of the defaults to log Gaussian Cox process models
* implementing relaxed fitting of regularised models

# multipointR 0.99.6
* adapted vignette to new functionality with `splines` package
* also reducing time to build vignette
* enable elastic net fits to regularise complex models e.g. with many spline
coefficients

# multipointR 0.99.5
* change from `deparse` to `deparse1` to accommodate long formulae

# multipointR 0.99.4
* more conservative estimation of the window from the min and max coordinates

# multipointR 0.99.3
* addressed R CMD check note on `.data` being a global variabel. Imported now
from `rlang` explicitly

# multipointR 0.99.2
* linting and BioCstyle formatting

# multipointR 0.99.1
* adressing [BioC review](https://github.com/Bioconductor/BiocContributions/
issues/63#issuecomment-4833222722)

# multipointR 0.99.0
* initial BioC submission

