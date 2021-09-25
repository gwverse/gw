# from ff_play_v3_gwCore_functions.R
# gw_get_nearby - returns index of nearby locations for any observation location
# gw_get_weight - returns the weight function (adaptive, kernel) 
# gw_do_weight - applies the weight function and rturns a weight matrix of n obs rows and n bws columns

# see https://ourcodingclub.github.io/tutorials/writing-r-package/
# library(devtools)
# setwd("/Users/geoaco/Desktop/my_docs_mac/leeds_work/research/gw_verse/gw_verse_R/gw_/R)
# load_all(".")
# setwd("/Users/geoaco/Desktop/my_docs_mac/leeds_work/research/gw_verse/gw_verse_R/gw_")
# library(roxygen2); # Read in the roxygen2 R package
# roxygenise();

#' Nearby observations 
#'
#' Returns a function for determining an index of nearby observations for a specific location and for a given bandwidth. 
#' The returned function is expected to be used within a 'function factory' approach. 
#' The returned function requires the following inputs: a distance matrix (`dist_mat`), an adaptive or fixed bandwidth, `bw`, and a location for which the index is calculated (`obs_index`).
#' @param adative A logical value TRUE or FALSE to indicate whether an adaptive or fixed bandwidth distance is being used.
#' @return A function that returns an index of nearby observations, given  an observation, (`obs_index`), a distance matrix (`dist_mat`), within the specified bandwidth, `bw`, for a given location (`obs_index`).
#' @examples
#' # load some packages and data
#' library(tmap)
#' data(georgia)
#' # define a distance matrix, a location and an adaptive bandwidth
#' dist_mat = as.matrix(dist(st_coordinates(st_centroid(georgia)), upper = T, diag = T))
#' obs_index = 30
#' bw = 30
#' # create the function
#' nearby_func = gw_get_nearby(adaptive = TRUE)
#' # apply to get an index of locations
#' index = nearby_func(obs_index, dist_mat, bw)
#' # map the result
#' tm_shape(georgia)+tm_polygons()+
#' tm_shape(georgia[index,])+tm_fill("red", alpha = 0.5)+
#' tm_shape(georgia[obs_index,])+tm_fill("black") 
#' @export
gw_get_nearby = function(adaptive){
	if( adaptive) return(function(obs_index, dist_mat, bw) {(sort(order(dist_mat[obs_index, ])[1:bw]))} )
	if(!adaptive) return(function(obs_index, dist_mat, bw) {(which(dist_mat[obs_index,] < bw))} )
}

#' Geographic weighting function
#'
#' Returns a function for calculating weights for each set of nearby observations for a given bandwidth, as identified by the function returned from `gw_get_nearby`. 
#' The returned function is expected to be used within a 'function factory' approach. 
#' The returned function requires the following to have been defined: a vector of distances (`dists`) such as are returned by indexing a distance matrix to get nearby observations (as determined by the function returned by `gw_get_nearby`) for a specific observation location, and an adaptive or fixed bandwidth, `bw`.
#' @param kernel The type of distance weighting to be used: one of "bisquare", "gaussian", "exponential", "tricube" or "boxcar".
#' @param adaptive A logical value (`TRUE` or `FALSE`) to indicate whether an adaptive or fixed bandwidth distance is being used.
#' @return A function that returns a vector weights given an adaptive or fixed bandwdth, `bw` and a vector of distances weighted by the specified `kernel` parameter.
#' @examples
#' # load some packages and data
#' library(tmap)
#' data(georgia)
#' # define a distance matrix, a location and an adaptive bandwidth
#' dist_mat = as.matrix(dist(st_coordinates(st_centroid(georgia)), upper = T, diag = T))
#' obs_index = 70
#' bw = 30
#' # create the nearby function - see the help for `gw_get_nearby`
#' nearby_func = gw_get_nearby(adaptive = TRUE)
#' # apply to get an index of locations and get a vector of distances
#' index = nearby_func(obs_index, dist_mat, bw)
#' dists = dist_mat[obs_index,index]
#' # create the weighting function and weight the nearby locations
#' weight_func = gw_get_weight(kernel = "bisquare", adaptive = TRUE)
#' # apply the weight function to the distances
#' w = weight_func(bw, dists)
#' # map the result
#' g2 = georgia[index,]
#' g2$weight = w
#' tm_shape(georgia)+tm_borders()+
#' tm_shape(g2)+tm_fill("weight", palette = "Reds")+
#' tm_shape(georgia[obs_index,])+tm_borders(lwd = 2) 
#' @export
gw_get_weight = function(kernel, adaptive) {
	if(kernel == "bisquare") {	
		if(adaptive) {
			return(function(bw, dists) {
				bw = max(dists)
				(1-(dists/bw)^2)^2
				}
			)
		} else {
			return(function(bw, dists) {
				(1-(dists/bw)^2)^2
				}
			)
		}
	}
	if(kernel == "gaussian") {	
		if(adaptive) {
			return(function(bw, dists) {
				bw = max(dists)
				{exp(-.5*(dists/bw)^2)}
				}
			)
		} else {
			return(function(bw, dists) {
				{exp(-.5*(dists/bw)^2)}
				}
			)
		}
	}
	if(kernel == "exponential") {	
		if(adaptive) {
			return(function(bw, dists) {
				bw = max(dists)
				{exp(-dists/bw)}
				}
			)
		} else {
			return(function(bw, dists) {
				{exp(-dists/bw)}
				}
			)
		}
	}
	if(kernel == "tricube") {	
		if(adaptive) {
			return(function(bw, dists) {
				bw = max(dists)
				{(1-(dists/bw)^3)^3}
				}
			)
		} else {
			return(function(bw, dists) {
				{(1-(dists/bw)^3)^3}
				}
			)
		}
	}
	if(kernel == "boxcar") 	
		return(function(bw, dists) rep(1, length(dists)))
}

#' Application of the geographic weighting function
#'
#' This function is for use within the application related function, for example to undertake a GWR. The idea here is that rather than working with a single observation location, the analysis is procededing with a vector of locations and multiple indices of nearby locations in matrix form for adaptive bandwidths (ie the nearest *n* observation) or in list form for fixed bandwidths (ie the observations within a specified distance of the different locations, whose number may vary from location to location).
#' This assembles the  outputs of `gw_get_weight` and `gw_get_nearby` so that they can be used in an `apply` function across the whole spatial dataset rather than for just a single observation at a time. 
#' @param i The i^{th} observation (as `obs_index` in `gw_get_weight` and `gw_get_nearby`)
#' @param bw The bandwidth (fixed or adaptive).
#' @param nearbyMat A `matrix` or `list` of indices of nearby locations for each observation.
#' @param dist_mat An n * n matrix of distances of each location.
#' @param weight_func The weighting function returned by `gw_get_weight`.
#' @return A vector of weights, for use in for example a locally weighted regression model.
#' @examples
#' # load some packages and data
#' library(tmap)
#' data(georgia)
#' # define a distance matrix, a location and an adaptive bandwidth
#' dist_mat = as.matrix(dist(st_coordinates(st_centroid(georgia)), upper = T, diag = T))
#' # create the nearby function - see the help for `gw_get_nearby`
#' nearby_func = gw_get_nearby(adaptive = TRUE)
#' # create the weighting function - see the help for `gw_get_weight`
#' weight_func = gw_get_weight(kernel = "bisquare", adaptive = TRUE)
#' # create an index of observations
#' indexMat  = matrix(1:nrow(georgia), ncol = 1) 	
#' # determine nearby locations for all observation (list for adaptive, matrix for fixed bandwidth) 
#' nearbyMat = apply(indexMat, 1, function(x) 
#'   nearby_func(x, dist_mat, bw))
#' # apply the weight function to the nearby locations 
#' weightMat = apply(indexMat, 1, function(x) 
#'   gw_do_weight(x, bw, nearbyMat, dist_mat, weight_func))
#' # examine
#' dim(weightMat)
#' weightMat[, 70]
#' # map the result
#' georgia$weight = weightMat[, 70]
#' tm_shape(georgia)+tm_polygons("weight", palette = "Reds")+
#' tm_shape(georgia[70,])+tm_borders(lwd = 2) 
#' @export
gw_do_weight = function(i, bw, nearbyMat, dist_mat, weight_func) {
	# list test for fixed and adaptive indices
	if (is.list(nearbyMat)) {
		nearby_locs = as.vector(nearbyMat[[i]]) 
		} else {
			nearby_locs = as.vector(nearbyMat[,i])			
		}
	dists = as.vector(dist_mat[i, nearby_locs])
	w = weight_func(bw, dists)
	w.vec = rep(0, nrow(dist_mat))
	w.vec[nearby_locs] = w
	return(w.vec)
}


