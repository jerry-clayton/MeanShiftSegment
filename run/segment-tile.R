###
# Set parameters for MS segmentation from command line args
# Leaving ms parameters hardcoded for now intentionally; they are set for this site already. 
###


args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

#Fraction of Cores to use
FRAC.CORES <- as.numeric(args[3])

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
library(plyr)


print(paste("input:",input_file))

# use all available processor cores
set_lidr_threads(0)

###
##   Begin running AMS3D on file
###

# read file as data.table
f_dt <- readRDS(input_file)
dtname <- names(f_dt)[1]
f_dt <- as.data.table(f_dt[[dtname]])
print('data.table: ')
print(f_dt)
pc.list <- list()
pc.list[[1]] <- f_dt
# this just has to be here
lib_path <- .libPaths()[1]
print('enter MS apply')
# catch errors in files
tryCatch({
  # important that we only use 90% of the available processor cores to prevent crashing
    # centroid accuracy is worth considering, default is 2
	ms_result <- MeanShiftR::apply_MeanShift(pc.list, lib.path = lib_path, run.parallel= FALSE, frac.cores = FRAC.CORES, version = 'classic', H2CW = H2CW, H2CL = H2CD,minz=MINZ)

     print("Mean Shift segmentation done")
     
     # assert 10 point minimum for segmented trees; drop everything with fewer than 10
     byid <- ms_result[, .(.N), by = ID]
     g10 <- byid[N>10,,]
     tg10 <- ms_result[ID %in% g10$ID]
     ms_result <- tg10 


    saveRDS(ms_result, output_file)
    #  print("joining ms result to original data")
     
    # # make IDs from xyz coords
    #  f_dt[, concat := paste(X,Y,Z, sep = "_")]
    #  ms_result[, concat := paste(X,Y,Z, sep = "_")]

    #  ms_result[, ID := ID]
    #  # left join (update-by-reference-join) (stackoverflow)
    #  # adds treeID to original lidR payload

    #  f_dt[ms_result, on = "concat", ID := ID]
 
    #  # save meanshift-generated IDs to LAS format
    #  # specifically: use original header with updated data and then force the header to update with add_lasattribute_manual()
    #  flas <- LAS(f_dt, f_head)
    #  flas <- add_lasattribute_manual(flas, f_dt[,ID], name = "ID", desc = "tree ID", type = "int64", NA_value = 99999)


    #  # write, and then read, the LAS file;
    #  # this makes it more straightforward to save the IDs after merging, and it enables us to reclaim all the memory used to generate the tree IDs
    #  writeLAS(flas, output_file)

    #  print("segmented las written")

##### 
##   End AMS3D segmentation, begin merging
#####
     # clear all objects but retain functions,
     # run garbage collection to ensure memory is reclaimed 
 #     rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])
 #     gc()
    
 #     # reload the command line file args                      
 #     args <- commandArgs(trailingOnly = TRUE)
 #     input_file <- args[1]
 #     output_file <- args[2]

                           
 #     flas <- readLAS(output_file)
 #     f_head <- readLASheader(output_file)
    
 #     # merge trees from segmented las
 #     mergedLAS <- mergeBelowThreshold(f_head, flas, 100)
    
 #     print("saving merged file")
 #     # overwrite the earlier segmentation with our merged result
 #     writeLAS(mergedLAS, output_file)
     },
     error = function(e){
     	print("error in file:")
	print(e)
     	print(input_file)
     })
