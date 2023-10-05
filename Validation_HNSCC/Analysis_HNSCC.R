# Load necessary packages
pkgs2load <- c('readxl','plyr', 'dplyr', 'magrittr', 'tidyr', 'flowCore', 'FlowSOM', 
               'data.table', 'Rtsne', 'ggplot2', 'gridExtra', 
               'ggpubr', 'readr', 'RColorBrewer','openxlsx',
               'matrixStats','reshape2','limma','ggrepel',
               'ConsensusClusterPlus','pheatmap','ggridges','grid','corrplot','progress','knitr')
sapply(pkgs2load, require, character.only = TRUE)

# Data set has to be downloaded from: 
# https://github.com/StiThor/Development-of-antibody-panel-for-imaging-mass-cytometry-to-investigate-cancer-associated-fibroblast

mySamples <- c('HNSCC01','HNSCC02','HNSCC03','HNSCC04','HNSCC05','HNSCC06','HNSCC07','HNSCC08','HNSCC09','HNSCC10')

# STEP 1 Data integration: read and concatenate the segmented single cells into one data frame 
AllData <- data.frame()
for(myname in mySamples){
  
  myStr<-paste('Processing sample: ',myname,sep='')
  print(myStr)
  # Data and code should be in same folder. If not modify accordingly 
  inputFileName <- paste(myname,'.csv',sep='')
  dataValues <- as.data.frame(read_csv(inputFileName)) 
  # List of antibodies. Name of each antibody must be the same as in .csv file(s)
  
  markers <- c('Cell_10_152Sm_FAP',
               'Cell_14_158Gd_ECad',
               'Cell_11_153Eu_ITGA11',
               'Cell_12_155Gd_FOXP3',
               'Cell_13_156Gd_CD4',
               'Cell_15_159Tb_CD68',
               'Cell_16_161Dy_CD20',
               'Cell_17_162Dy_CD8',
               'Cell_18_163Dy_TenascinC',
               'Cell_19_167Er_GranzymeB',
               'Cell_1_141Pr_aSMA',
               'Cell_20_168Er_Ki67',
               'Cell_21_169Tm_Collagen1',
               'Cell_22_172Yb_FSP1',
               'Cell_23_175Lu_CD146',
               'Cell_24_176Yb_CD90',
               'Cell_2_142Nd_EGFR',
               'Cell_3_143Nd_Podoplanin',
               'Cell_4_144Nd_YAP1',
               'Cell_5_145Nd_Caveolin',
               'Cell_6_146Nd_CD16',
               'Cell_7_147Sm_CD163',
               'Cell_8_149Sm_CD140b',
               'Cell_9_151Eu_CD31',
               'Area',
               'Eccentricity',
               'Extent',
               'Number_Neighbors',
               'X_position',
               'Y_position'
               )
  
  data <- dataValues[,markers]
  
  # Give each antibody their correct name
  colnames(data) <- c('FAP',
                      'ECad',
                      'ITGA11',
                      'FOXP3',
                      'CD4',
                      'CD68',
                      'CD20',
                      'CD8',
                      'TenascinC',
                      'GranzymeB',
                      'aSMA',
                      'Ki67',
                      'Collagen1',
                      'FSP1',
                      'CD146',
                      'CD90',
                      'EGFR',
                      'Podoplanin',
                      'YAP1',
                      'Caveolin',
                      'CD16',
                      'CD163',
                      'CD140b',
                      'CD31',
                      'Area',
                      'Eccentricity',
                      'Extent',
                      'No_of_Neighbors',
                      'X',
                      'Y')
  
  # Store all data for further analysis
  a <- data.frame(File=myname,NumCells=nrow(data),data)
  AllData <- rbind(AllData,a)
}

#STEP2: Data transformation and normalization  
# This is only censoring, and other transformation is not performed
# Remove unnecessary columns from data and convert to matrix
a <- as.matrix(AllData[, 3:(ncol(AllData) - 6)])

# Calculate the 99th percentile for each column in the subset
myQuant <- apply(a, 2, function(x) quantile(x, probs = 0.99))

# Initialize a progress bar with the total number of rows in AllData
pb <- progress_bar$new(total = nrow(AllData))

# Loop through each row in AllData
for (i in 1:nrow(AllData)) {
  k <- 1
  # Loop through the columns in the selected range
  for (j in 3:(ncol(AllData) - 6)) {
    # Check if the value in AllData is greater than the corresponding 99th percentile
    if (AllData[i, j] > myQuant[k]) {
      # Replace the value with the 99th percentile
      AllData[i, j] <- myQuant[k]
    }
    k <- k + 1
  }
  # Increment the progress bar
  pb$tick()
}

# Generate plot to inspect total number of cells in each sample
plot_data <- data.frame(table(AllData$File))
colnames(plot_data)[1] <- 'Sample'
colnames(plot_data)[2] <- 'NumCells'

p1 <- ggplot(plot_data) +
  aes(x = Sample, y = NumCells, fill = Sample) +
  geom_bar(stat = 'identity', width = 0.28) +
  geom_text(aes(label = NumCells), vjust = -0.5, size = 3) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 10),
    axis.title.y = element_text(size = 10),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 12, face = 'bold', lineheight = 1),
    legend.position = 'none'
  )

myfile <- paste('Summary_NumCells.pdf', sep = '')
pdf(myfile)
print(p1)
dev.off()

# STEP3: Clustering the data using FlowSOM and ConsensusClusterPlus
# Similar clustering results can be obtained using Phenograph or any other algorithm
# Initialize a grid of relatively big size compared to the available data dimension
n_clusters_x <- 15
n_clusters_y <- 15
n_meta_clusters <- 200
#  Use the 24 antibodies as input. Shape features could also be included if wanted 
selected_markers<- c('FAP',
                     'ITGA11',
                     'FOXP3',
                     'CD4',
                     'ECad',
                     'CD68',
                     'CD20',
                     'CD8',
                     'TenascinC',
                     'GranzymeB',
                     'aSMA',
                     'Ki67',
                     'Collagen1',
                     'FSP1',
                     'CD146',
                     'CD90',
                     'EGFR',
                     'Podoplanin',
                     'YAP1',
                     'Caveolin',
                     'CD16',
                     'CD163',
                     'CD140b',
                     'CD31')

data1 <- as.matrix(AllData[,selected_markers])
# Run flowSOM and consensusClusterPlus
clustering_data <- flowFrame(data1)
start_time = Sys.time()
fsom1 <- FlowSOM(clustering_data,
                 colsToUse = selected_markers,
                 xdim = n_clusters_x, ydim = n_clusters_y,
                 nClus=n_meta_clusters,
                 maxMeta = n_meta_clusters, 
                 rlen = 10,transform = FALSE)
end_time = Sys.time()
end_time-start_time
codes <- fsom1$map$codes
nmc <- n_meta_clusters
mc <- ConsensusClusterPlus(t(codes), maxK = nmc, reps = 100, 
                           pItem = 0.9, pFeature = 1, 
                           clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average", 
                           distance = "euclidean", seed = 1234,plot='pdf')

#res <- calcICL(mc,plot='pdf')
results <- mc
maxK <- nmc
# STEP 4: Estimate the quality of metaclustering given the PAC metric. Other metrics could also be used.
############## PAC implementation ##############
Kvec = 2:maxK
x1 = 0.1; x2 = 0.9 # Threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="") # From 2 to maxK
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}#end for i

# Plot the proportion of ambiguous clustering (PAC) values
df <- data.frame(PAC)
df$K <- c(2:(length(PAC)+1))
o <- ggplot(df, aes(x=K,y=PAC))+
  geom_point(size = 0.68)+
  geom_line()+
  theme_classic() +
  geom_vline(xintercept = 130,linetype="dotted", 
             color = 'red', size=0.88)+
  geom_hline(yintercept = 0.01,linetype="dotted", 
             color = 'red', size=0.88)+
  ylab('Proportion of ambiguous clustering (PAC) ')+
  xlab('Number of metaclusters ')+
  # Scale_color_gradient(colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))+
  theme(axis.text.y = element_text( size = 12 ),
        axis.text.x = element_text( size = 12 ),
        axis.title.y = element_text( size = 12 ),
        axis.title.x = element_text( size = 12 ),strip.text = element_text(size = 12,face='bold',lineheight=1),
        legend.position = "bottom",
        plot.title = element_text(size = 12, face = "bold",hjust = 0.5))+
  guides(colour = guide_legend(override.aes = list(size=4.8)))

myfile<-paste('Clustering_PAC_A.pdf',sep ='')
pdf(myfile,onefile = TRUE)
print(o)
dev.off()

# After visual inspection its necessary to set PAC cutoffs and decide about merging clusters
# Here we use the value of 130, as it is the first point under 1% PAC, but other solutions are also acceptable depending on the level of resolution
nmc <- 130
code_clustering1 <- mc[[nmc]]$consensusClass
clusters <- code_clustering1[fsom1$map$mapping[,1]]
cluster_num <- max(clusters)
# Produce Z scored values
data <- AllData[,3:(ncol(AllData)-6)] 
data_Z <- scale(as.matrix(data))
data_Z <- data_Z[,selected_markers]

# Identify the z-transformed centroids of clusters
centroids<-matrix(0L, nrow = cluster_num, ncol = (ncol(data_Z)))
for (i in 1:cluster_num) {
  
  b<-which(clusters %in% i)
  d<-data_Z[b,1:ncol(data_Z)]
  #df <- data[b,2:(ncol(data)-6)]
  d<-colMeans(d)
  centroids[i,]<-d
}
# Estimate the cluster abundance
clustering_table <- as.numeric(table(clusters))
clustering_table<-100*clustering_table/sum(clustering_table)
clustering_table<-round(clustering_table,2)
# Fix the centroids phenograph row and column labels
b1 <- as.vector(paste('C',1:cluster_num,' (', clustering_table,'%)',sep=''))
b2 <- as.vector(paste(1:cluster_num,sep=''))
c1 <- colnames(data_Z)[1:(ncol(data_Z))]
rownames(centroids) <- b1
# Save data and write excel file with results for annotation
colnames(centroids) <- c(c1)
centroids <- as.data.frame(centroids)
write.xlsx(centroids, "Clusters.xlsx",asTable = FALSE,overwrite=TRUE)


################################################################
# At this step it is necessary to read and annotate the clusters using cell type defining markers
# An example of the annotated file is provided in file Clusters.xlsx that can be downloaded from our repository
################################################################
# Load the annotated clusters -- specify the unassigned
inputFileName <- "Clusters.xlsx"
manualLabels <- read_excel(inputFileName)
manualLabels <- manualLabels[,c('Cluster','Annotations')]
plot_data <- AllData
myClust <- paste('C',clusters,sep='')
plot_data$Clusters <- myClust
dr <- match(plot_data$Clusters, manualLabels$Cluster)
plot_data$CellType <- manualLabels$Annotations[dr]

# Estimate summary statistics about the phenotype abundance of CAFs
samples_count <- plot_data %>% 
  group_by(CellType,File) %>% 
  summarise(samples = n())%>% 
  mutate(perc = samples/sum(samples))

# Generate exploratory plots
samples_count$CellType <- factor(samples_count$CellType,levels=c('Ecad+ Epithelial cells',
                                                                 'CAF-1',
                                                                 'CAF-2',
                                                                 'Pericytes',
                                                                 'CD4+ T-cells',
                                                                 'CD8+ T-cells',
                                                                 'FOXP3+ Treg',
                                                                 'FOXP3+ CD4+ T-cells',
                                                                 'CD20+ B-cells',
                                                                 'CD68+ Macrophages',
                                                                 'T/B cells',
                                                                 'CD31+ Endothelial cells',
                                                                 'PDPN+ Ly Endothelial cells',
                                                                 'Unassigned'))

# Figure with proportion of different cell phenotype in each sample
o1 <- ggplot(samples_count, aes(x=CellType,y=perc,group=File,fill=File))+
  geom_bar(stat = 'identity',width=0.25)+
  ylab('Relative percentage (%)')+
  theme_bw() +
  theme(axis.text.y = element_text( size = 10 ),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8.8),
        axis.title.y = element_text( size = 10 ),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12,face='bold',lineheight=1),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(size = 18, face = "bold",hjust = 0.5))

myfile<-paste('/Summary_clusters_relative_perc.pdf',sep='')
pdf(myfile,onefile = TRUE)
print(o1)
dev.off()

# Figure with counts of each cell type in each sample
o2 <- ggplot(samples_count, aes(x=CellType,y=(samples),group=File,fill=File))+
  geom_bar(stat = 'identity',width=0.25)+
  ylab('Number of cells')+
  theme_bw() +
  theme(axis.text.y = element_text( size = 10 ),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8.8),
        axis.title.y = element_text( size = 10 ),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12,face='bold',lineheight=1),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(size = 18, face = "bold",hjust = 0.5))

summary_table <- samples_count %>%
  group_by(CellType, File) %>%
  summarize(TotalSamples = sum(samples)) %>%
  arrange(CellType)

# Print table with exact numbers of each cell type for each sample + write information as .csv file
kable(summary_table, caption = "Summary Statistics for Cell Types (Aggregated with Sample Info)", format = "markdown")
write.csv(summary_table, file = "summary_table.csv", row.names = FALSE)

myfile<-paste('Summary_clusters_counts.pdf',sep='')
pdf(myfile,onefile = TRUE)
print(o2)
dev.off()

# Calculate the total count for each cell type
totals <- samples_count %>%
  group_by(CellType) %>%
  summarize(Total = sum(samples))

# Merge the totals back into the original data frame
samples_count <- samples_count %>%
  left_join(totals, by = "CellType")

# Create the ggplot
o3 <- ggplot(samples_count, aes(x = CellType, y = Total, fill = CellType)) +
  geom_bar(stat = 'identity', width = 0.25) +
  geom_text(aes(label = Total), vjust = -0.5, size = 3, position = position_dodge(width = 0.25)) +
  scale_fill_manual(limits = levels(samples_count$CellType), values = c('yellow',
                                                                        'tomato1',
                                                                        'red4',
                                                                        'rosybrown4',
                                                                        'olivedrab1',
                                                                        'olivedrab4',
                                                                        'seagreen1',
                                                                        'seagreen3',
                                                                        'skyblue2',
                                                                        'royalblue2',
                                                                        'royalblue4',
                                                                        'sienna',
                                                                        'darkmagenta',
                                                                        'gray80')) +
  ylab('Total Number of Cells') +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8.8),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12, face = 'bold', lineheight = 1),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5))

myfile <- paste('Summary_cellType_counts.pdf', sep = '')
pdf(myfile, onefile = TRUE)
print(o3)
dev.off()

# Calculate summary statistics for each CellType
summary_table <- samples_count %>%
  group_by(CellType) %>%
  summarize(
    Count = sum(samples),
    Mean = mean(samples),
    Median = median(samples),
    SD = sd(samples),
    Min = min(samples),
    Max = max(samples)
  ) %>%
  arrange(CellType)

# Print the summary table
kable(summary_table, caption = "Summary Table for Cell Type Counts", format = "markdown")

# STEP5: Visualize the clusters in heatmap
data <- AllData[,3:(ncol(AllData)-6)] 
data_Z <- scale(as.matrix(data))
data_Z <- data_Z[,selected_markers]
cluster_names <- unique(plot_data$CellType)

# Identify the z-transformed centroids of clusters
centroids<-matrix(0L, nrow = length(cluster_names), ncol = (ncol(data_Z)))
myC <- 1
myClusterNames <- c('Ecad+ Epithelial cells',
                    'CAF-1',
                    'CAF-2',
                    'Pericytes',
                    'CD4+ T-cells',
                    'CD8+ T-cells',
                    'FOXP3+ Treg',
                    'FOXP3+ CD4+ T-cells',
                    'CD20+ B-cells',
                    'CD68+ Macrophages',
                    'T/B cells',
                    'CD31+ Endothelial cells',
                    'PDPN+ Ly Endothelial cells',
                    'Unassigned')
  
for (i in myClusterNames) {
  
  b<-which(plot_data$CellType==i)
  d<-data_Z[b,1:ncol(data_Z)]
  #df <- data[b,2:(ncol(data)-6)]
  d<-colMeans(d)
  centroids[myC,]<-d
  myC <- myC+1
}
# Estimate the abundance of annotated cell types
clustering_table <- table(plot_data$CellType)
clustering_table<-100*clustering_table/sum(clustering_table)
clustering_table<-round(clustering_table,2)
clustering_table <- as.data.frame(clustering_table)
b1 <- myClusterNames
b2 <- paste('',clustering_table$Freq,'%',sep='')
c1 <- colnames(data_Z)[1:(ncol(data_Z))]
rownames(centroids) <- b2
colnames(centroids) <- c1

# Generate info for the row annotation in heatmap
my_row_annot <- data.frame(CellType=b1)
rownames(my_row_annot)<- b2
colnames(my_row_annot)[1] <- 'CellType'
#specify the colors
my_colour = list(
  CellType = c('Ecad+ Epithelial cells' = 'yellow',
               'CAF-1' = 'tomato1',
               'CAF-2' = 'red4',
               'Pericytes' = 'rosybrown4',
               'CD4+ T-cells' = 'olivedrab1',
               'CD8+ T-cells' = 'olivedrab4',
               'FOXP3+ Treg' = 'seagreen1',
               'FOXP3+ CD4+ T-cells' = 'seagreen3',
               'CD20+ B-cells' = 'skyblue2',
               'CD68+ Macrophages' = 'royalblue3',
               'T/B cells' = 'royalblue4',
               'CD31+ Endothelial cells' = 'sienna',
               'PDPN+ Ly Endothelial cells' = 'darkmagenta',
               'Unassigned' = 'gray80'))

breaksList<-seq(min(centroids),max(centroids),by=0.5)
# Plot the cluster centroids
myfile<-paste('Heatmap_cellTypes.pdf',sep ='')
pheatmap(centroids[,1:ncol(data_Z)],labels_row=b2,labels_col = c1,
         annotation_row = my_row_annot,annotation_colors = my_colour,
         cluster_rows=TRUE,cluster_cols = TRUE ,fontsize = 6,
         display_numbers = TRUE, number_color = "black",fontsize_number = 4,
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)),
         filename=myfile,
         cutree_rows = 4,cutree_cols=5,
         cellheight = 16, cellwidth = 16,angle_col = 45)


# STEP6: Perform dimensional reduction using t-sne

df <- AllData[,selected_markers]
duplicates <- duplicated(df[,1:ncol(df)]) | duplicated(df[,1:ncol(df)], fromLast = TRUE) 
tsne_out_high<-Rtsne(df[!duplicates,1:ncol(df)],max_iter=3000,perplexity=30,seed=1234,num_threads=4)

tsne_data_3 <- as.data.frame(data_Z[!duplicates,])
tsne_data_3$TSNE1 <- tsne_out_high$Y[,1]
tsne_data_3$TSNE2 <- tsne_out_high$Y[,2]
tsne_data_3$Sample <- AllData[!duplicates,'File']
myLabel <- plot_data$CellType
tsne_data_3$CellType <- myLabel[!duplicates]
tsne_data_3$CellType<- factor(tsne_data_3$CellType,levels=c('Ecad+ Epithelial cells',
                                                            'CAF-1',
                                                            'CAF-2',
                                                            'Pericytes',
                                                            'CD4+ T-cells',
                                                            'CD8+ T-cells',
                                                            'FOXP3+ Treg',
                                                            'FOXP3+ CD4+ T-cells',
                                                            'CD20+ B-cells',
                                                            'CD68+ Macrophages',
                                                            'T/B cells',
                                                            'CD31+ Endothelial cells',
                                                            'PDPN+ Ly Endothelial cells',
                                                            'Unassigned'))

# Visualize the cell types in 2D
o1 <- ggplot(tsne_data_3,  aes(x = TSNE1, y = TSNE2, color = CellType )) +
  geom_point(size = 0.28) +
  theme_bw() +
  scale_color_manual(limits = levels(tsne_data_3$CellType),values = c('yellow',
                                                                      'tomato1',
                                                                      'red4',
                                                                      'rosybrown4',
                                                                      'olivedrab1',
                                                                      'olivedrab4',
                                                                      'seagreen1',
                                                                      'seagreen3',
                                                                      'skyblue2',
                                                                      'royalblue3',
                                                                      'royalblue4',
                                                                      'sienna',
                                                                      'darkmagenta',
                                                                      'gray80'))+
  guides(fill = FALSE)+
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 3))+
  theme(axis.text.y = element_text( size = 8 ),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8.8),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12,face='bold',lineheight=1),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(size = 18, face = "bold",hjust = 0.5))

myfile<-paste('Figure_tsne.pdf',sep='')
pdf(myfile,onefile = TRUE)
print(o1)
dev.off()

# Similar t-sne but split based on the sample information
o2 <- ggplot(tsne_data_3,  aes(x = TSNE1, y = TSNE2, color = CellType )) +
  geom_point(size = 0.28) +
  theme_bw() +
  scale_color_manual(limits = levels(tsne_data_3$CellType),values = c('yellow',
                                                                      'tomato1',
                                                                      'red4',
                                                                      'rosybrown4',
                                                                      'olivedrab1',
                                                                      'olivedrab4',
                                                                      'seagreen1',
                                                                      'seagreen3',
                                                                      'skyblue2',
                                                                      'royalblue3',
                                                                      'royalblue4',
                                                                      'sienna',
                                                                      'darkmagenta',
                                                                      'gray80'))+
  facet_wrap(~Sample,ncol=2)+
  guides(fill = FALSE)+
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 6))+
  theme(axis.text.y = element_text( size = 8 ),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8.8),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12,face='bold',lineheight=1),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(size = 18, face = "bold",hjust = 0.5))

myfile<-paste('Figure_tsne_facets',sep='')
pdf(myfile,onefile = TRUE)
print(o2)
dev.off()

o3 <- ggplot(tsne_data_3,  aes(x = TSNE1, y = TSNE2, color = Sample )) +
  geom_point(size = 0.28) +
  theme_bw() +
  guides(fill = FALSE)+
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 2))+
  theme(axis.text.y = element_text( size = 8 ),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8.8),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12,face='bold',lineheight=1),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(size = 18, face = "bold",hjust = 0.5))


myfile<-paste('Figure_tsne_samples.pdf',sep='')
pdf(myfile,onefile = TRUE)
print(o3)
dev.off()

# Make a combined t-sne and color points based on antibody intensity
myC <- 1
n <- list()
for(myMarker in selected_markers){
  
  myD <- tsne_data_3[,c('Sample','TSNE1','TSNE2',myMarker)]
  colnames(myD)[4] <- 'Intensity'
  n[[myC]] <- ggplot(myD,  aes(x = TSNE1, y = TSNE2, color = Intensity)) +
    geom_point(size = 0.108) +
    theme_bw() +
    #facet_wrap(~ Sample, nrow = 2) +
    scale_color_gradientn("Intensity",colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))+
    #ylim(-98, 98)+
    #xlim(-98, 98)+
    guides(fill = FALSE)+
    ggtitle(myMarker)+
    scale_shape_manual(values = c(1,5),name='Cell type')+
    guides(shape = guide_legend(override.aes = list(size = 2), ncol = 5))+
    theme(axis.text.y = element_text( size = 8 ),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8.8),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          strip.text = element_text(size = 12,face='bold',lineheight=1),
          legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          plot.title = element_text(size = 18, face = "bold",hjust = 0.5))
  myC <- myC + 1
  
}
# Combine individual plots to one
o4 <- ggarrange(n[[1]],n[[2]],n[[3]],n[[4]],n[[5]],n[[6]],n[[7]],n[[8]],n[[9]],n[[10]],
                n[[11]],n[[12]],n[[13]],n[[14]],n[[15]],n[[16]],n[[17]],n[[18]],n[[19]],n[[20]],
                n[[21]],n[[22]],n[[23]],n[[24]],
                ncol = 3,nrow = 3)

myfile<-paste('Figure_tsne_markers',sep='')
pdf(myfile,onefile = TRUE)
print(o4)
dev.off()

# Save plot data for downstream visualization of tissues 
save(plot_data, file = "plot_data.RData")
