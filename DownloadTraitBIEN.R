require(BIEN)


rangesBIEN <- BIEN_ranges_list()
sppnames <- as.character(rangesBIEN$species)

sppnames2 <- gsub("_", " ", sppnames)

traits <- BIEN_trait_list()



BIEN_trait_mean(species=c("Poa annua","Juncus trifidus"),trait="leaf dry mass per leaf fresh mass") 


sppSLA <- list()
sppWPH <- list()
sppSeed <- list()

for(r in 1:length(sppnames)){
  svMisc::progress(r, max.value = (length(sppnames)))
  sppSLA[[r]] <- BIEN_trait_mean(sppnames2[r], trait = "leaf area per leaf dry mass")
  sppWPH[[r]] <- BIEN_trait_mean(sppnames2[r], trait = "whole plant height")
  sppSeed[[r]] <- BIEN_trait_mean(sppnames2[r], trait = "seed mass")
}

sppSLAAll <- do.call(rbind, sppSLA)
head(sppSLAAll)

sppWPHAll <- do.call(rbind, sppWPH)
head(sppWPHAll)

sppSeedAll <- do.call(rbind, sppSeed)
head(sppSeedAll)

save.image("TraitsBIEN.RData")

sppRanges <- list()
for(i in 30001:length(sppNames)){
  svMisc::progress(i, max.value = (length(sppNames)-30001))
  sppRanges[[i]] <- BIEN_ranges_species(sppNames[i])
}




