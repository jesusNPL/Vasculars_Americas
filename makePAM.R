
makePAM <- function (path, xmn = -180, xmx = 180, ymn = -90, ymx = 90, resol = 1, 
        remove.cells = TRUE, remove.sp = TRUE, show.matrix = FALSE, 
        crs = CRS("+proj=longlat +datum=WGS84"), cover = 0, presence = NULL, 
        origin = NULL, seasonal = NULL, count = FALSE) 
{
  
  if ( ! ("maptools" %in% installed.packages())) {install.packages("maptools", dependencies = T)}
  if ( ! ("letsR" %in% installed.packages())) {install.packages("letsR", dependencies = T)}
  if ( ! ("rgdal" %in% installed.packages())) {install.packages("rgdal", dependencies = T)}
  
require(maptools)
require(letsR)
require(rgdal)
  
  shapes <- list.files(path, pattern = "shp$", full.names = T, 
                       recursive = T)
  r <- raster(xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
  res(r) <- resol
  values(r) <- 0
  valores <- values(r)
  xy <- xyFromCell(r, 1:length(valores))
  nomes <- numeric(length(shapes))
  matriz <- matrix(0, nrow = nrow(xy), ncol = length(shapes))
  matriz <- cbind(xy, matriz)
  n <- length(shapes)
  k <- 0
  if (cover > 0) {
    if (!(xmn == -180 & xmx == 180 & ymn == -90 & ymx == 
            90)) {
      grid <- rasterToPolygons(r)
      areagrid <- try(areaPolygon(grid), silent = TRUE)
    }
    else {
      areagrid <- values(area(r)) * 1e+06
    }
    if (class(areagrid) == "try-error") {
      areagrid <- values(area(r)) * 1e+06
    }
  }
  if (count) {
    dev.new(width = 2, height = 2, pointsize = 12)
    par(mar = c(0, 0, 0, 0))
    for (j in 1:n) {
      plot.new()
      text(0.5, 0.5, paste(paste("Total:", n, "\n", "Runs to go: ", 
                                 (n - j))))
      valores2 <- valores
      shp <- readShapePoly(shapes[j], delete_null_obj = TRUE, 
                           force_ring = T)
      nomes[j] <- levels(shp$Species)[1]
      shp <- lets.shFilter(shp, presence = presence, origin = origin, 
                           seasonal = seasonal)
      if (!is.null(shp)) {
        k <- k + 1
        cell <- extract(r, shp, cellnumber = T, small = T, 
                        weights = T)
        cell <- cell[!sapply(cell, is.null)]
        if (length(cell) > 0) {
          cell <- lapply(cell, function(x) {
            colnames(x) <- 1:3
            return(x)
          })
        }
        cell2 <- do.call(rbind.data.frame, cell)
        if (cover == 0) {
          cell3 <- cell2[which(cell2[, 3] >= cover), 
                         ]
        }
        else {
          areashape <- areaPolygon(shp)
          prop <- numeric()
          for (k1 in 1:length(cell)) {
            prop <- c(prop, cell[[k1]][, 3] * areashape[k1]/areagrid[cell[[k1]][, 
                                                                                1]])
          }
          if (any(prop > 1)) {
            prop[prop > 1] <- 1
          }
          cell3 <- cell2[which(prop >= cover), ]
        }
        valores2[cell3[, 1]] <- 1
        matriz[, (j + 2)] <- valores2
      }
    }
    dev.off()
  }
  if (!count) {
    for (j in 1:n) {
      valores2 <- valores
      shp <- readShapePoly(shapes[j], delete_null_obj = TRUE, 
                           force_ring = T)
      nomes[j] <- levels(shp$Species)[1]
      shp <- lets.shFilter(shp, presence = presence, origin = origin, 
                           seasonal = seasonal)
      if (!is.null(shp)) {
        k <- k + 1
        cell <- extract(r, shp, cellnumber = T, small = T, 
                        weights = T)
        cell <- cell[!sapply(cell, is.null)]
        if (length(cell) > 0) {
          cell <- lapply(cell, function(x) {
            colnames(x) <- 1:3
            return(x)
          })
        }
        cell2 <- do.call(rbind.data.frame, cell)
        if (cover == 0) {
          cell3 <- cell2[which(cell2[, 3] >= cover), 
                         ]
        }
        else {
          areashape <- areaPolygon(shp)
          prop <- numeric()
          for (k1 in 1:length(cell)) {
            prop <- c(prop, cell[[k1]][, 3] * areashape[k1]/areagrid[cell[[k1]][, 
                                                                                1]])
          }
          if (any(prop > 1)) {
            prop[prop > 1] <- 1
          }
          cell3 <- cell2[which(prop >= cover), ]
        }
        valores2[cell3[, 1]] <- 1
        matriz[, (j + 2)] <- valores2
      }
    }
  }
  if (k == 0) {
    stop("after filtering none species distribution left")
  }
  colnames(matriz) <- c("Longitude(x)", "Latitude(y)", nomes)
  riqueza <- rowSums(as.matrix(matriz[, -c(1, 2)]))
  if (remove.cells) {
    matriz <- removeCells(matriz)
  }
  if (remove.sp) {
    matriz <- removeSp(matriz)
  }
  matriz <- unicas(matriz)
  if (show.matrix) {
    return(matriz)
  }
  else {
    values(r) <- riqueza
    final <- list(Presence_and_Absence_Matrix = matriz, Richness_Raster = r, 
                  Species_name = (colnames(matriz)[-(1:2)]))
    class(final) <- "PresenceAbsence"
    return(final)
  }
}


removeCells <- function(x) {
  rem <- which(rowSums(x[, -c(1, 2), drop = FALSE]) == 0)
  if (length(rem) > 0) {
    x <- x[-rem, , drop = FALSE]
  }
  if(nrow(x) == 0) {
    stop("No cells left after removing cells without occurrences")
  }
  return(x)
}

removeSp <- function(x) {  
  
  rem <- which(colSums(x[, -(1:2), drop = FALSE]) == 0) + 2
  
  if (length(rem) > 0) {
    x <- x[, -rem, drop = FALSE]
  }
  
  if (ncol(x) == 2) {
    stop("No species left after removing species without occurrences")
  }
  
  return(x)
}

unicas <- function(resu) {
  nomes <- colnames(resu)
  
  if(!any(duplicated(nomes))) {
    return(resu)
  } else {
    n <- ncol(resu)
    for(i in 1:(n - 1)){  
      nome1 <- nomes[i]
      for(j in 1:n){
        nome2 <- nomes[j]    
        if (nome1 == nome2) {
          divid <- which((resu[, i] != 0 & resu[, j] != 0))    
          resu[,i] <- resu[, i] + resu[, j]
          resu[divid, i] <- resu[divid, i]/2
        }
      }
    }
    pos <- duplicated(nomes)
    resu <- resu[, !pos]
    if (is.vector(resu)) {
      nomes <- names(resu)
      resu <- matrix(resu, ncol = length(resu))
      colnames(resu) <- nomes                
    }
    return(resu)
  }
}

##### USAGE #####
#PAM_vasculars <- makePAM("C:/Users/jpintole/Documents/Jesus/Ranges_vasculars/", count = T)
