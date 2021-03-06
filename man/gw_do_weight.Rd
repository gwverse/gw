% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gw_core_functions.R
\name{gw_do_weight}
\alias{gw_do_weight}
\title{Application of the geographic weighting function}
\usage{
gw_do_weight(i, bw, nearbyMat, dist_mat, weight_func)
}
\arguments{
\item{i}{The i^{th} observation (as `obs_index` in `gw_get_weight` and `gw_get_nearby`)}

\item{bw}{The bandwidth (fixed or adaptive).}

\item{nearbyMat}{A `matrix` or `list` of indices of nearby locations for each observation.}

\item{dist_mat}{An n * n matrix of distances of each location.}

\item{weight_func}{The weighting function returned by `gw_get_weight`.}
}
\value{
A vector of weights, for use in for example a locally weighted regression model.
}
\description{
This function is for use within the application related function, for example to undertake a GWR. The idea here is that rather than working with a single observation location, the analysis is procededing with a vector of locations and multiple indices of nearby locations in matrix form for adaptive bandwidths (ie the nearest *n* observation) or in list form for fixed bandwidths (ie the observations within a specified distance of the different locations, whose number may vary from location to location).
This assembles the  outputs of `gw_get_weight` and `gw_get_nearby` so that they can be used in an `apply` function across the whole spatial dataset rather than for just a single observation at a time.
}
\examples{
# load some packages and data
library(tmap)
data(georgia)
# define a distance matrix, a location and an adaptive bandwidth
dist_mat = as.matrix(dist(st_coordinates(st_centroid(georgia)), upper = T, diag = T))
# create the nearby function - see the help for `gw_get_nearby`
nearby_func = gw_get_nearby(adaptive = TRUE)
# create the weighting function - see the help for `gw_get_weight`
weight_func = gw_get_weight(kernel = "bisquare", adaptive = TRUE)
# create an index of observations
indexMat  = matrix(1:nrow(georgia), ncol = 1) 	
# determine nearby locations for all observation (list for adaptive, matrix for fixed bandwidth) 
nearbyMat = apply(indexMat, 1, function(x) 
  nearby_func(x, dist_mat, bw))
# apply the weight function to the nearby locations 
weightMat = apply(indexMat, 1, function(x) 
  gw_do_weight(x, bw, nearbyMat, dist_mat, weight_func))
# examine
dim(weightMat)
weightMat[, 70]
# map the result
georgia$weight = weightMat[, 70]
tm_shape(georgia)+tm_polygons("weight", palette = "Reds")+
tm_shape(georgia[70,])+tm_borders(lwd = 2) 
}
