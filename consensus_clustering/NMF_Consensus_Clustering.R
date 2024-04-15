#Running NMF consensus clustering on significant genomic events such as mutations and copy number changes generated from DNA-Seq to generate prognostic sutbypes.
#Author: Sunandini Sharma

########This module load all libraries #########################################
mirror = 'http://cran.mirrors.hoobly.com/'
install.packages('doMC',repos='http://R-Forge.R-project.org')
install.packages('Cairo',repos=mirror)
install.packages('gtools',repos=mirror)
install.packages('gdata',repos=mirror)
install.packages('caTools',repos=mirror)
install.packages('gplots',repos=mirror)
install.packages("dplyr")
install.packages("pdftools")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

library("doMC")
library("Cairo")
library("gtools")
library("gdata")
library("caTools")
library("gplots")
library("dplyr")
library("pdftools")
library("ComplexHeatmap")
library("stringr")

################################################################################
#Run the NMF consensus clustering adapted from Chapuy, Stewart, & Dunford et al. "Molecular subtypes of diffuse large B cell lymphoma are associated with distinct pathogenic mechanisms and outcomes." Nature Medicine, 30 April 2018. https://doi.org/10.1038/s41591-018-0016-8"

setwd("~/consensus_clustering")
source('src/GDAC_TopgenesforCluster/Topgenes_v1.R')
result <- main("-s./src/GDAC_TopgenesforCluster/","-minput_data/WES_iteration_145.txt","-uALL","-oDLBCL")
source("src/GDAC_NmfConsensusClustering/GDAC_CNMF.R", echo = TRUE)
result <- main("-ssrc/GDAC_NmfConsensusClustering/","-mDLBCL.expclu.gct","-oDLBCL",'-u4','-v10') 
source("src/GDAC_selectBestcluster/select_best_cluster_chc.R")
result <- main("-uDLBCL.expclu.gct","-mPearson","-cDLBCL.cophenetic.coefficient.txt","-wDLBCL.membership.txt","-voutput_dir/DLBCL","-pinput_data/DLBCL_significant_event_matrix.txt")

################################################################################
#This module uses the output from the above module and assigns membership to each patient and perform statistical test to get subtype enrichment based on user selected membership threshold. For this script, I used membership.3 as an example.
# Remove X and - in sample names part of data cleaning
setwd("~/consensus_clustering")
data1<-read.table("DLBCL.membership.txt", header=TRUE, sep="\t", row.names=1)
rownames(data1)<-gsub("^X", "", rownames(data1))
rownames(data1)<-gsub("-", ".", rownames(data1))
memb<- "membership.3"
data1 <- data1[order(data1[, memb]), ]
data2 <- read.table("WES_iteration_145.txt", header=TRUE, sep="\t", check.names=FALSE, row.names=1)
m<-match(rownames(data1), colnames(data2))
data3<-data2[,m]
data5<-data.frame()
class<-unique(data1[, memb])

for(i in 1:length(class)){
  l<-c()
  m<-c()
  n<-c()
  o<-c()
  names<-c()
  ################################################################
  a<-which(data1[,memb] == class[i])
  ################################################################  
  start<-a[1]
  end<-a[length(a)]
  ######This loop performs frequency calculation of each genomic events across clusters and uses fisher's test to get final subtype enrichment/depletion
  for(j in 1:nrow(data3)){
    row1<-data3[j, start:end]
    row2<-data3[j, -c(start:end)]
    altClust1<- row1 %>% select_if(function(col) col > 0)
    notaltClust1<- row1 %>% select_if(function(col) col == 0)
    altClustother<-row2 %>% select_if(function(col) col > 0)
    notaltClustother<- row2 %>% select_if(function(col) col == 0)
    l<-c(l, ncol(altClust1))
    m<-c(m, ncol(notaltClust1))
    n<-c(n, ncol(altClustother))
    o<-c(o, ncol(notaltClustother))
    names<-c(names, rownames(data3)[j])
    
  }
  data4<-data.frame(names, l,m, n, o)
  write.csv(data4, file=paste0(i, "_frequency.csv"), row.names=FALSE)
  rnames<-c()
  pvalue<-c()
  odds.ratio<-c()
  for(k in 1:nrow(data4)){
    cont_table=matrix(c(data4[k,2], data4[k,3], data4[k,4], data4[k,5]), nrow=2)
    test<-fisher.test(cont_table)
    pvalue<-c(pvalue, test$p.value)
    odds.ratio<-c(odds.ratio, test$estimate)
    rnames<-c(rnames, data4[k,1])
  }
  result<-data.frame(rnames, pvalue, odds.ratio)
  result<-subset(result, result$pvalue <= 0.1)
  result<-result[order(result$pvalue),]
  result$cluster=c(rep(i, times = nrow(result)))
  write.csv(result, file=paste0(i, "_significant_events.csv"), row.names=FALSE)
  data5<-rbind(data5, result)
  write.csv(data5, file="All.significant.events.csv", row.names=FALSE)
  
}

#retain only one of the duplicate events with highest odds.ratio across different clusters
data<-read.csv("All.significant.events.csv", header=TRUE, sep=",")
data_ordered<- data[order(data$odds.ratio, decreasing = TRUE), ]
data_highest<-data_ordered[!duplicated(data_ordered$rnames),]
data_highest
data_highest<-data_highest[order(data_highest$cluster),]
write.csv(data_highest, file="unique_across_groups.csv", row.names=FALSE)

#Arrange samples and events for heatmap
match_row<-match(data_highest[,1], rownames(data3))
data6<-data3[match_row,]
write.csv(data6, file="Heatmap.csv", row.names=TRUE)

#Read Heatmap,csv to do graphing################################################
df<-read.csv("Heatmap.csv", header=TRUE, sep=",", row.names=1, check.names=FALSE)
amplification_rows <- grepl("Amplification", rownames(df), ignore.case = TRUE)
df[amplification_rows, ] <- lapply(df[amplification_rows, ], function(x) {
  x[x == 2] <- "Broad_Amp;"
  x[x == 1] <- "Low_Amp;"
  x[x == 0] <- ""
  return(x)
})

#deletion
deletion_rows <- grepl("Deletion", rownames(df), ignore.case = TRUE)
df[deletion_rows, ] <- lapply(df[deletion_rows, ], function(x) {
  x[x == 1] <- "Loss;"
  x[x == 0] <- ""
  return(x)
})

#Broad deletion
del <- grepl("Broad_Deletion", rownames(df), ignore.case = TRUE)
df[del, ] <- lapply(df[del, ], function(x) {
  x[x == 1] <- "Broad_Deletion;"
  x[x == 0] <- ""
  return(x)
})

rows_to_modify <- which(!(grepl("Amplification|Deletion|Broad_Deletion|Broad_Gain", rownames(df), ignore.case = TRUE)))

# Replace values based on conditions
df[rows_to_modify, ] <- lapply(df[rows_to_modify, ], function(x) {
  x[x == 0] <- ""
  x[x > 0] <- "MUT;"
  return(x)
})

##Broad Gain
bg <- grepl("Broad_Gain", rownames(df), ignore.case = TRUE)
df[bg, ] <- lapply(df[bg, ], function(x) {
  x[x == 1] <- "Broad_Gain;"
  x[x == 0] <- ""
  return(x)
})

write.csv(df, file="check.csv")

###################This function assigns color to each genomic event############
alter_fun_list = list( 
  
  background = function(x, y, w, h) {
    
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.9, "mm"), gp = gpar(fill = "white", col = NA))
    
  },
  
  Low_Amp = function(x, y, w, h) {
    
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.9, "mm"), gp = gpar(fill = "salmon1", col = NA))
    
  },
  Broad_Amp = function(x, y, w, h) {
    
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.9, "mm"), gp = gpar(fill = "red", col = NA))
  },
  MUT = function(x, y, w, h) {
    
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.9, "mm"), gp = gpar(fill = "black", col = NA))
  },
  Loss = function(x, y, w, h) {
    
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.9, "mm"), gp = gpar(fill = "blue", col = NA))
  },
  Broad_Gain = function(x, y, w, h) {
    
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.9, "mm"), gp = gpar(fill = "maroon1", col = NA))
  },
  Broad_Deletion = function(x, y, w, h) {
    
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.9, "mm"), gp = gpar(fill = "cyan", col = NA))
  }
  
)

col = c("Low_Amp"="salmon1", "Broad_Amp"="red", "MUT"="black", "Loss"="blue", "Broad_Gain"="maroon1", "Broad_Deletion"= "cyan")


# Read the data from CSV files
data3 <- read.csv("check.csv", header = TRUE, sep = ",", check.names = FALSE, row.names = 1)
data4 <- read.csv("unique_across_groups.csv", header = TRUE, sep = ",", check.names = FALSE)

# Remove leading 'X' and replace '-' with '.' in row names
colnames(data3) <- gsub("^X", "", colnames(data3))
colnames(data3) <- gsub("-", ".", colnames(data3))

# Initialize vectors to store final column and row orders
final_cols <- c()
final_rows <- c()

# Loop through clusters
data1$clust<-data1[,memb]
k<-unique(data1$clust)

for (i in 1:length(k)) {
  # Subset data based on cluster
  subset_data <- data3[, data1$clust == k[i]]
  ids_to_extract <- which(data4$cluster == k[i])
  subset_data <- subset_data[ids_to_extract, ]
  
  # Define alter_fun_list and col as you did in the original code
  
  # Create oncoPrint plot
  h <- oncoPrint(subset_data, get_type = function(x) strsplit(x, ";")[[1]], alter_fun = alter_fun_list, col = col, remove_empty_columns = FALSE, row_names_gp = gpar(fontsize = 5), alter_fun_is_vectorized = FALSE, show_column_names = TRUE, column_names_gp = gpar(fontsize = 4))
  
  # Save oncoPrint plot as a PDF
  pdf_file <- paste0(k[i], ".pdf")
  pdf(pdf_file, height = 6, width = 6)
  print(h)  # Print the plot
  dev.off()
  
  # Extract text content from the PDF
  pdf_text_content <- pdf_text(pdf_file)
  
  cleaned_text <- gsub("\\s+", " ", pdf_text_content)
  
  # Extract columns and rows from the cleaned text
  cols <- colnames(data3)
  rows <- rownames(data3)
  
  # Create regex pattern with word boundaries
  cols_pattern <- paste0("\\b(", paste(cols, collapse = "|"), ")\\b")
  rows_pattern <- paste0("\\b(", paste(rows, collapse = "|"), ")\\b")
  
  extracted_cols <- unlist(str_extract_all(cleaned_text, cols_pattern))
  extracted_rows <- unlist(str_extract_all(cleaned_text, rows_pattern))
  
  # Append to final lists
  final_cols <- c(final_cols, extracted_cols)
  final_rows <- c(final_rows, extracted_rows)
  
}

final_cols<-data.frame(final_cols)
final_rows<-data.frame(final_rows)

data3<-read.csv("check.csv", header=TRUE, sep=",", row.names=1, check.names=FALSE)
p<-match(rownames(data3), final_rows[,1])
data3<-data3[p,]
pdf("new10.pdf", height = 15, width = 20)
oncoPrint(data3, get_type = function(x) strsplit(x, ";")[[1]], alter_fun = alter_fun_list, col = col, remove_empty_columns = FALSE, remove_empty_rows = FALSE, row_names_gp = gpar(fontsize = 15), alter_fun_is_vectorized = FALSE, show_column_names = TRUE, column_names_gp = gpar(fontsize = 8), column_order = 1:nrow(final_cols), row_order = 1:nrow(final_rows))
dev.off()

##################Store Results at a particular membership #####################
name<-memb
dir.create(paste0("~/consensus_clustering/", name, sep="_", "dir"))
z<-paste0("~/consensus_clustering/", name, sep="_", "dir")
source_directory <- "~/consensus_clustering/"
destination_directory <- z
files_to_copy <- list.files(source_directory, pattern = "^DLBCL", full.names = TRUE)
file.copy(files_to_copy, destination_directory, overwrite = FALSE)
file.remove(files_to_copy)


