# Working with ENMeval

library(ENMeval)
library(tidyverse)
library(raster)
library(terra)

# Loading occ points 

## Preparing occ
occ <- select(spatial.thinned[[8]], longitude, latitude)

## Preparing envs
file_list <- list.files("Bioclims/",
                         pattern = ".asc", 
                         full.names = T) 
envs <- stack(file_list)


ENMevaluate(occs = occ,
            envs = envs,
            bg = NULL,
            tune.args = list(fc = c("L","Q"),rm = 1:3),
            partitions = "checkerboard1",
            algorithm = "maxent.jar",
            partition.settings = NULL ,
            other.settings = NULL,
            categoricals = NULL,
            doClamp = TRUE,
            clamp.directions = ,
            user.enm = NULL,
            user.grp = NULL,
            occs.testing = NULL,
            taxon.name = NULL,
            n.bg = 10000,
            overlap = FALSE,
            overlapStat = c("D", "I"),
            user.val.grps = NULL,
            user.eval = NULL,
            rmm = NULL,
            parallel = FALSE,
            numCores = NULL,
            parallelType = "doSNOW",
            updateProgress = FALSE,
            quiet = FALSE,
            occ = NULL,
            env = NULL,
            bg.coords = NULL,
            RMvalues = NULL,
            fc = NULL,
            occ.grp = NULL,
            bg.grp = NULL,
            method = NULL,
            bin.output = NULL,
            rasterPreds = NULL,
            clamp = NULL,
            progbar = NULL) -> enm.eval
