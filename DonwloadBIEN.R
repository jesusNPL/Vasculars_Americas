setwd("C:/Users/jpintole/Desktop/RangesBIEN")


library(BIEN)
library(sp)
library(maptools)
temp_dir <- file.path(tempdir(), "BIEN_temp")#Set a working directory
x <- BIEN_ranges_species("Tabebuia_roseoalba")
BIEN_ranges_genus("Tabebuia")
sppShapes <- list()
spnames <- c("Tabebuia_roseoalba", "Tabebuia_aurea", "Tabebuia_sp.")
for(i in 1:length(spnames)){
  sppShapes[[i]] <- BIEN_ranges_species()
}


##### Download available species ranges from BIEN #####
## BIEN_ranges_list a data.frame containing listing all range maps currently available.

rangesBIEN <- BIEN_ranges_list()
sppnames <- as.character(rangesBIEN$species)

sppRanges2 <- list()
for(r in 20001:length(sppnames)){
  svMisc::progress(r, max.value = (length(sppnames)-20001))
  sppRanges2[[r]] <- BIEN_ranges_species(sppnames[r])
}

sppRangesAll <- do.call(rbind, sppRanges)
head(sppRangesAll)

save.image("rangesBIEN.RData")
