#Identify spots near edges of box and tissue
EdgeDetectionVisium <- function(coord_path, seurat.obj = NULL,
                                search = 'radius', neighbors = 7) {
  #data_dir - directory that has tissue_positions.csv file
  #seurat.obj - matching seurat object
  require(RANN)
  require(stringr)
  require(dplyr)
  require(tidyr)
  
  #finds and reads in the data
  coords <- read.delim(paste(coord_path,
                             list.files(coord_path)[which(str_detect(
                               list.files(coord_path),
                               'tissue_position')==TRUE)], sep = '/'),
                       header = F, sep = ',')
  colnames(coords) <- c('barcode', 'in_tissue', 'array_row', 'array_col',
                        'pxl_row_in_fullres', 'pxl_col_in_fullres')
  
  if(is.null(seurat.obj) == FALSE) {
    #double checks everything is in the right order
    coords <- coords[match(colnames(seurat.obj), coords$barcode),]
  }
  
  rownames(coords) <- coords$barcode
  closest <- as.data.frame(nn2(data=coords[,c(3,4)], k=7, searchtype = 'radius',
                               radius = 2)[[1]])
  closest[closest == 0] <- NA
  
  rownames(closest) <- coords$barcode
  closest$barcode <- rownames(closest)
  closest2 <- closest %>% drop_na()
  closest3 <- closest %>% filter(!barcode %in% closest2$barcode)
  coords$Filter <- coords$barcode
  coords[rownames(closest2),]$Filter <- 'Keep'
  coords[closest3$barcode,]$Filter <- 'Filter'
  
  #Now do it again
  coords$Filter2 <- coords$Filter
  coords_red <- coords %>% filter(Filter == 'Keep')
  
  closest <- as.data.frame(nn2(data=coords_red[,c(3,4)], k=neighbors, 
                               searchtype = search,
                               radius = 2)[[1]])
  closest[closest == 0] <- NA
  
  rownames(closest) <- coords_red$barcode
  closest$barcode <- rownames(closest)
  closest2 <- closest %>% drop_na()
  closest3 <- closest %>% filter(!barcode %in% closest2$barcode)
  coords[closest3$barcode,]$Filter2 <- 'Filter'
  
  # Now Again
  coords$Filter3 <- coords$Filter2
  coords_red <- coords %>% filter(Filter2 == 'Keep')
  
  closest <- as.data.frame(nn2(data=coords_red[,c(3,4)], k=7, searchtype = 'radius',
                               radius = 2)[[1]])
  closest[closest == 0] <- NA
  
  rownames(closest) <- coords_red$barcode
  closest$barcode <- rownames(closest)
  closest2 <- closest %>% drop_na()
  closest3 <- closest %>% filter(!barcode %in% closest2$barcode)
  coords[closest3$barcode,]$Filter3 <- 'Filter'
  
  # Last Time
  coords$Filter4 <- coords$Filter3
  coords_red <- coords %>% filter(Filter3 == 'Keep')
  
  closest <- as.data.frame(nn2(data=coords_red[,c(3,4)], k=7, 
                               searchtype = 'radius',
                               radius = 2)[[1]])
  closest[closest == 0] <- NA
  
  rownames(closest) <- coords_red$barcode
  closest$barcode <- rownames(closest)
  closest2 <- closest %>% drop_na()
  closest3 <- closest %>% filter(!barcode %in% closest2$barcode)
  coords[closest3$barcode,]$Filter4 <- 'Filter'
  
  return(coords)
}
