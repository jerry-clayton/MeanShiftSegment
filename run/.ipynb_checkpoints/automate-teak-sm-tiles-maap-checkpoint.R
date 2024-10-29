install.packages('lidR', repos="https://cloud.r-project.org")
install.packages('devtools', repos="https://cloud.r-project.org")
install.packages('data.table', repos="https://cloud.r-project.org")
install.packages('sf', repos="https://cloud.r-project.org")
library(devtools)
#install_version('rgdal','1.6-7', repos="https://cloud.r-project.org")
install_version('rgeos','0.6-4', repos="https://cloud.r-project.org")
install_github('niknap/MeanShiftR')


###
# Set parameters for MS segmentation from command line args
# Leaving ms parameters hardcoded for now intentionally; they are set for this site already. 
###


args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

# Plot width for subplots
P_WIDTH <- as.numeric(args[3])
# Plot buffer for subplots
P_BUFFER <- as.numeric(args[4])

#Fraction of Cores to use
FRAC.CORES <- as.numeric(args[5])

# Height-to-crown-depth ratio for Mean Shift Algorithm
H2CD <- 0.85
# Height-to-crown-width ratio for Mean Shift Algorithm
H2CW <- 0.15
# Minimum Z coordinate. Everything below this is discarded before running segmentation. 
MINZ <- 1

library(lidR)
library(MeanShiftR)
library(data.table)
library(sf)


print(paste("input:",input_file))

# use all available processor cores
set_lidr_threads(0)

#####
##   Begin function definitions for merging
#####

compute_min_centroid_dist <- function(cm){
    
	# Initialize a vector to store the minimum distance for each centroid
	centroids <- cm
    min_distances <- numeric(length(centroids))

  # Loop through each centroid and compute distances to all other centroids
  for (i in 1:length(centroids)) {

    # Get the centroid at index i
    current_centroid <- centroids[i, ]

    # Compute distances from current centroid to all other centroids
    distances <- st_distance(current_centroid, centroids)

    # Set the distance to itself (index i) to a very large number to exclude it
    distances[i] <- Inf

    # Get the minimum distance to another centroid
    min_distances[i] <- min(distances)
	}

    return(min_distances)
}

h2cw_ratio <- function(height){
	return(10^(-0.1*(height^0.61)))
}

RMSE <- 1.1231

# gets the first neighbor inside the regression radius based on the X-Y centroid.
# in the future should probably evaluate candidates based on Z centroid too.

getMergeID <- function(treeID, centroid, regRadius, cm){
	returnID <- NULL
	for(tree in cm$ID){
		if(tree == treeID){
		## do not test itself
		  next
		}

		tree <- cm[ID == tree,,]
		distance <- as.numeric(st_distance(centroid, tree$centroid.geometry))
		if(distance < regRadius){
			returnID <- tree$ID
			break
		}
	}
	return(returnID)
	
}

mergeTreesByID <- function(currentID, neighborID, las_dt){

	## assign all points with currentID to have neighborID
	las_dt[ID == currentID, ID := neighborID]
	
	return(las_dt)
}

update_crown_metrics <- function(head, las_dt, cm, currentID, neighborID){
	## update cm to do the following:
		## remove currentID entry
		## re-compute the centroid, area, etc for neighborID
		newPts <- las_dt[ID == neighborID,,]
		cm_new <- crown_metrics(lidR::LAS(newPts, head), geom = "concave", func = NULL, attribute = "ID")
     		cm_new$area <- st_area(cm_new)
     		cm_new$centroid <- st_centroid(cm_new)
		currTree <- cm[ID == currentID,,]
		newTree <- cm[ID == neighborID,,]
		newMaxZ <-  max(currTree$maxZ, newTree$maxZ)
		cm_new$maxZ <- newMaxZ
    		cm_new$minDiam <- (h2cw_ratio(cm_new$maxZ)*cm_new$maxZ)
     		cm_new$minArea <- ((cm_new$minDiam/2)^2)*pi
		#print("new formatted row: ")
		#print(cm_new)
		#print("old formatted row: ")
		#print(newTree)
		#remove current ID
		cm <- cm[ID != currentID]
		cm_new <- as.data.table(cm_new)
		#add this new version of the merged tree as a row to cm and return it
		cm[ID == neighborID, names(cm_new) := cm_new]
		return(cm)

}

mergeBelowThreshold <- function(head, las, heightThreshold){
    
     # save NAs separately to add back later
     las_nas <- filter_poi(las, ID == 99999)
     
     # remove NAs
     las <- filter_poi(las, ID != 99999)

     las_dt <- lidR::payload(las) %>% as.data.table

     #filter to only include trees below the height threshold
     
     maxHeightByID <- las_dt[, .(max(Z)), by = ID]
     treesBelowThreshold <- maxHeightByID[V1 < heightThreshold,,]
     las_subset <- filter_poi(las, ID %in% treesBelowThreshold$ID)

     #compute X-Y polygon of each tree, get its area and centroid for merging ops, and add the max-height
     cm <- crown_metrics(las_subset, attribute = "ID", geom = "concave", func = NULL)
     cm$area <- st_area(cm)
     cm$centroid <- st_centroid(cm)
     cm <- cm %>% as.data.table	
     cm[treesBelowThreshold, maxZ := V1, on = "ID"]
     cm$minDiam <- (h2cw_ratio(cm$maxZ)*cm$maxZ)
     cm$minArea <- ((cm$minDiam/2)^2)*pi
    # print(cm)
     ##
     ##  Wrap everything above here into a separate function

     ## add the min area as a column

     dynamic_trees <- lidR::payload(las) %>% as.data.table
     originalIDs <- cm$ID

     for (tree in 1:length(originalIDs)){
	tree <- originalIDs[tree]
	treeStats <- cm[ID == tree,,]
	## merge if the segmented area is less than the area of the regression equation	
	if(as.numeric(treeStats$area) < as.numeric(treeStats$minArea)){
	## merging
 	mergeID <-getMergeID(tree, treeStats$centroid.geometry, treeStats$minDiam/2, cm)

	#mergeID <-getMergeIDwithAngle(tree, treeStats$centroid.geometry, las_dt, treeStats$minDiam/2, cm)
	if(is.null(mergeID)){
		print(paste0('no merge candidate found for tree ID: ',tree))
		next
	}
	print(paste0('treeID: ',tree,' merge ID found:',mergeID))
	las_dt <- mergeTreesByID(tree, mergeID, las_dt)
	#print('finished mergeTreesByID')
	cm <- update_crown_metrics(head, las_dt, cm, tree, mergeID)
	#minD <- compute_min_centroid_dist(cm$centroid.geometry)
	#centroids <- compute_3d_centroids(las_dt)
	#print("3d centroids")
	#print(centroids)
	#minD3d <- compute_min_centroid_dist(centroids)
	#print("minimum distances for XYZ centroids")
	#print(minD3d)
	#print("minimum distances from XY centroids")
	#print(minD)
	#df <- data.frame( ID = originalIDs,
	#		 Min_XY_centroid_dist = minD,
	#		 Min_XYZ_centroid_dist =minD3d)
	#print(df)
	#write.csv(df, 'teak_043_min_centroid_dist.csv')
	##print(length(originalIDs))
	next
	}
	else{ 
	  next
	}

     }
	## join back as LAS
     	# join las_dt and las_nas, then use LAS()
     	las_nas <- as.data.table(payload(las_nas))
     	# because the header already has 99999 set as the NA value, it seems to be excluding those points.
        # here we assign NA to the column, which will be replaced with 99999 by the LAS() function.
        las_nas[, ID := NA]
        joined <- rbind(las_nas,las_dt)
	#print(las_nas)
	flas <- LAS(joined, head)
	return(flas)
	#	saveRDS(las_dt, "merge_test_043_full.rds")
}


#####
##   End function definitions for merging
#####


###
##   Begin running AMS3D on file
###

# read header and file separately so that we can save the LAS later
f_head <- readLASheader(input_file)
f_las <- readLAS(input_file)

# convert LAS to data.table for MeanShiftR package
f_dt <- lidR::payload(f_las) %>% as.data.table

# subdivide the point cloud to parallelize
point_clouds <- MeanShiftR::split_BufferedPointCloud(f_dt, plot.width = P_WIDTH, buffer.width = P_BUFFER)
print("Point cloud generation done")

# this just has to be here
lib_path <- .libPaths()[1]

# catch errors in files
tryCatch({
  # important that we only use 90% of the available processor cores to prevent crashing
    
	ms_result <- MeanShiftR::parallel_MeanShift(point_clouds, lib.path = lib_path, frac.cores = FRAC.CORES, version = 'classic', H2CW = H2CW, H2CL = H2CD,minz=MINZ)

     print("Mean Shift segmentation done")
     
     # assert 10 point minimum for segmented trees; drop everything with fewer than 10
     byid <- ms_result[, .(.N), by = ID]
     g10 <- byid[N>10,,]
     tg10 <- ms_result[ID %in% g10$ID]
     ms_result <- tg10 
     
     print("joining ms result to original data")
     
    # make IDs from xyz coords
     f_dt[, concat := paste(X,Y,Z, sep = "_")]
     ms_result[, concat := paste(X,Y,Z, sep = "_")]

     ms_result[, ID := ID]
     # left join (update-by-reference-join) (stackoverflow)
     # adds treeID to original lidR payload

     f_dt[ms_result, on = "concat", ID := ID]
 
     # save meanshift-generated IDs to LAS format
     # specifically: use original header with updated data and then force the header to update with add_lasattribute_manual()
     flas <- LAS(f_dt, f_head)
     flas <- add_lasattribute_manual(flas, f_dt[,ID], name = "ID", desc = "tree ID", type = "int64", NA_value = 99999)


     # write, and then read, the LAS file;
     # this makes it more straightforward to save the IDs after merging, and it enables us to reclaim all the memory used to generate the tree IDs
     writeLAS(flas, output_file)

     print("segmented las written")

##### 
##   End AMS3D segmentation, begin merging
#####
     # clear all objects but retain functions,
     # run garbage collection to ensure memory is reclaimed 
     rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])
     gc()
    
     # reload the command line file args                      
     args <- commandArgs(trailingOnly = TRUE)
     input_file <- args[1]
     output_file <- args[2]

                           
     flas <- readLAS(output_file)
     f_head <- readLASheader(output_file)
    
     # merge trees from segmented las
     mergedLAS <- mergeBelowThreshold(f_head, flas, 100)
    
     print("saving merged file")
     # overwrite the earlier segmentation with our merged result
     writeLAS(mergedLAS, output_file)
     },
     error = function(e){
     	print("error in file:")
	print(e)
     	skipped[input_file]
     })
