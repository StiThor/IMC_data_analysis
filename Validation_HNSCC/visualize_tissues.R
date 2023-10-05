library(cytomapper)

# Read in the mask for HNSCC01
myMask <- loadImages('HNSCC01.tiff')
myMask <- scaleImages(myMask, 2^16-1)
# Read in images
myStack1 <- loadImages('HNSCC01_CD163.tiff')
myStack2 <- loadImages('HNSCC01_CD8.tiff')
myStack3 <- loadImages('HNSCC01_ECad.tiff')

mcols(myMask)$ImageNb <- c("HNSCC01")
mcols(myStack1)$ImageNb <- c("HNSCC01")
mcols(myStack1)$ImageName <- c("HNSCC01")
mcols(myStack2)$ImageNb <- c("HNSCC01")
mcols(myStack2)$ImageName <- c("HNSCC01")
mcols(myStack3)$ImageNb <- c("HNSCC01")
mcols(myStack3)$ImageName <- c("HNSCC01")

channelNames(myStack1) <- c("CD163",'a','b')
channelNames(myStack2) <- c("CD8",'a','b')
channelNames(myStack3) <- c("ECad",'a','b')

cur_channel1 <- getChannels(myStack1,'CD163')
cur_channel2 <- getChannels(myStack2,'CD8')
cur_channel3 <- getChannels(myStack3,'ECad')
#merge stacks
all_stacks <- mergeChannels(cur_channel1, cur_channel2)
all_stacks <- mergeChannels(all_stacks,cur_channel3)

plotPixels(image =all_stacks,colour_by = c('CD163','CD8','ECad'),img_id = "ImageNb")

sce <- measureObjects(myMask, all_stacks,img_id = "ImageNb")

A <- sce 
A$CellType <- 'OtherCells'

load('plot_data.RData')

a <- which(plot_data$File=='OPC01')
tmp <- plot_data[a,]
A$CellType <- tmp$CellType

A$CellNb <- 1:nrow(tmp)

#here is the old cell types, need to modify according to the new cell types
A$CellType[A$CellType!='Ecad+ Epithelial cells' &
             A$CellType!='CAF-1' &
             A$CellType!='CAF-2' &
             A$CellType!='CD4+ T-cells' &
             A$CellType!='CD8+ T-cells' &
             A$CellType!='CD4+ FOXP3+ Treg' &
             A$CellType!='FOXP3+ Treg' &
             A$CellType!='CD20+ B-cells' &
             A$CellType!='T/B cells' &
             A$CellType!='CD68+ Macrophages' &
             A$CellType!='PDPN+ Ly Endothelial cells' &
             A$CellType!='CD31+ Endothelial cells' &
             A$CellType!='Pericytes'] <- 'OtherCells'

#'gold','firebrick1','deepskyblue1','darkorchid2','forestgreen',
'gray0'

g1 <- plotCells(myMask, object = A,
                img_id = "ImageNb", cell_id = "CellNb",
                colour_by = c('CD163','CD8','ECad'),
                outline_by = "CellType",
                colour = list(CD163 = c("black", "red"),
                              CD8=c('black','blue'),
                              ECad=c('black','yellow'),
                              CellType = c(OtherCells = "gray88",
                                           `Ecad+ Epithelial cells` = "yellow",
                                           `CAF-1` = "red",
                                           `CAF-2` = "red4",
                                           `CD4+ T-cells` = "green",
                                           `CD8+ T-cells` = "green3",
                                           `CD4+ FOXP3+ Treg` = "forestgreen",
                                           `FOXP3+ Treg` = "olivedrab1",
                                           `CD20+ B-cells` = "blue",
                                           `T/B cells` = "blue4",
                                           `CD68+ Macrophages`="olivedrab3",
                                           `PDPN+ Ly Endothelial cells`="darkorchid2",
                                           `CD31+ Endothelial cells`="darkmagenta",
                                           `Pericytes`="darkseagreen1")),
                scale_bar = list(length = 10,
                                 cex = 1,
                                 lwidth = 10,
                                 colour = "white",
                                 position = "bottomleft",
                                 margin = c(5,5),
                                 frame = 3),
                image_title = list(text = c('OPC9'),
                                   position = "topleft",
                                   colour = "white",
                                   margin = c(2,10),
                                   font = 1,
                                   cex = 1),
                legend = list(colour_by.title.font = 2,
                              colour_by.title.cex = 1,
                              colour_by.labels.cex = 0.5
                ),
                return_plot = TRUE,thick=TRUE)

myfile<-paste('OPC1_tissuePlot_v1.pdf',sep='')
pdf(myfile,onefile = TRUE)
print(g1)
dev.off()

A$mock <- '1'

g2 <- plotCells(mask = myMask, object = A,
                cell_id = "CellNb", img_id = "ImageNb", 
                colour_by = "CellType",
                outline_by = "mock",
                colour = list(CellType = c(OtherCells = "ghostwhite",
                                           `Ecad+ Epithelial cells` = "yellow",
                                           `CAF-1` = "tomato1",
                                           `CAF-2` = "red4",
                                           `CD4+ T-cells` = "olivedrab1",
                                           `CD8+ T-cells` = "olivedrab4",
                                           `CD4+ FOXP3+ Treg` = "seagreen3",
                                           `FOXP3+ Treg` = "seagreen1",
                                           `CD20+ B-cells` = "skyblue2",
                                           `T/B cells` = "royalblue4",
                                           `CD68+ Macrophages`="royalblue3",
                                           `PDPN+ Ly Endothelial cells`="darkmagenta",
                                           `CD31+ Endothelial cells`="sienna",
                                           `Pericytes`="rosybrown"),
                              mock = c(`1`='gray3')),
                scale_bar = list(length = 10,
                                 cex = 1,
                                 lwidth = 10,
                                 colour = "white",
                                 position = "bottomleft",
                                 margin = c(5,5),
                                 frame = 3),
                image_title = list(text = c(''),
                                   position = "topleft",
                                   colour = "white",
                                   margin = c(2,10),
                                   font = 1,
                                   cex = 1),
                legend = list(colour_by.title.font = 2,
                              colour_by.title.cex = 1,
                              colour_by.labels.cex = 0.5
                ),
                return_plot = TRUE,thick=TRUE)


myfile<-paste('HNSCC01.pdf',sep='')
pdf(myfile,onefile = TRUE)
print(g2)
dev.off()