install.packages("gplots")
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
install.packages("grid")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
install.packages("colorRamp2")
install.packages("magick")
install.packages("gplots")
install.packages("heatmap3")
suppressPackageStartupMessages(library(dplyr))
install.packages("imager")
install.packages("tidyverse")
install.packages("fields")

library(xtable)
library(scales)
library("RColorBrewer")
library("gplots")
library("RColorBrewer")
library(lattice)
library(grid)
library(ComplexHeatmap)
library(colorRamp2)
library("gplots")
library("heatmap3")
library(purrr)
library(imager)
require(ggplot2)
library(fields)

suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))

#----------------------------------------------------------------------

files <- list.files("~/Desktop/Datasets", pattern = ".", all.files = FALSE, recursive = TRUE, full.names = TRUE)
files[1]
length(files)


dataset_total_number <- c()
dataset_id_column <- c()
zero_count_column <- c()
dropout_percentage <- c()
zero_increment <- c()

n_genes <- c()
n_samples <- c()

for (k in files) {
  k
  dataset <- readRDS(k)
  dataset_file <- basename(k) #get full dataset name
  dataset_id <- tools::file_path_sans_ext(dataset_file) #get full dataset name without extension

  dataset <- updateObject(dataset)
  experiments(dataset)
  (dataset <- experiments(dataset)[["gene"]])
  dataset <- assays(dataset)[["TPM"]]
  dim(dataset)
  
  nrow(dataset) 
  ncol(dataset) 
  
  i = 1
  j = 1
  nrow(dataset)
  ncol(dataset)
  zero_count = 0
  entry_count = 0
  percentage = 0
  
  progress <- txtProgressBar(min = 0,  max = nrow(dataset), style = 3,  char = "|") #progress bar
  
  for (i in 1:nrow(dataset)) {
    for (j in 1:ncol(dataset)) {
      entry_count = entry_count +1
     if(dataset[i,j] == "0") {
      zero_count = zero_count + 1
     }
    }
    setTxtProgressBar(progress, i)
  }
  #cat("\n", zero_count) 
  cat("\n") 
  percentage = zero_count / entry_count
  
  dataset_id_column <- append(dataset_id_column, dataset_id) 
  zero_count_column <- append(zero_count_column, zero_count) 
  dataset_total_number <- append(dataset_total_number, entry_count) 
  dropout_percentage <- append(dropout_percentage, percent(percentage, accuracy = 0.01)) 
  
  
  
}

zeros_table <- data.frame(
  Dataset = dataset_id_column,
  Total_entries = dataset_total_number,
  Score = zero_count_column,
  Dropout_percentage = dropout_percentage
  )

zeros_table$

col <- colorRamp2(breaks = c(0, 1), colors = c( "white", "dark red"))
pdf(file = paste0("t/", dataset_id, ".pdf"))
Heatmap(dataset[1:20,1:288], 
        name = "Gene Expression", 
         col = col,
         column_title = dataset_id, row_title = "Gene", 
         show_row_names = FALSE, show_column_names = F,
         #use_raster = FALSE,
         #k = 2, # k-means
         row_names_gp = gpar(fontsize = 7) #   Text size for row names
 )
 dev.off()



#Percentages of zeros
SMARTer <- c(74.81, 69.76, 82.09, 85.68, 34.61, 80.6, 81.16, 85.41, 90.15, 46.07, 44.98, 69.77, 73.46, 75.44, 75.26, 78.45, 80.40, 80.60)
Smart_Seq <- c(70.34, 91.64, 91.65, 66.83, 41.14, 82.21, 87.38, 79.92, 88.35, 85.03, 86.71)
Tang <- c(68.87, 62.50, 62.01, 82.15)
mean(SMARTer)
mean(Smart_Seq)
mean(Tang)

pdf(file = paste0("t/", "legend.pdf"))
boxplot(SMARTer)

all_smarter <- c(13157568, 822348, 62955308, 42716410, 4522914, 5478312, 18782784, 7812306, 1233522, 1224384, 6715842, 30000280, 1956540, 47285010, 34218814, 14299718, 4385856, 42456918)
all_zeros_for_smarter <- c(9760666, 573689, 5168291, 36598103, 1565596, 4415355, 15243952, 6672424, 1112033, 5641141, 3020700, 20930161, 1437367, 35669663, 25751841, 11218716, 352608, 34220939)

all_smart_seq <- c(99718322, 7923790, 34355872, 13294626, 56348352, 9274258, 11804458, 4934088, 13157568, 17543424, 38330554)
all_zeros_for_smart_seq <- c(70141082, 7261491, 31486897, 8884261, 23180225, 7624434, 10314995, 3943204, 11624995, 14916651, 33235315)

all_tang <- c(1891322, 776662, 776662, 21391504)
all_zeros_for_tang <- c(1302647, 485439, 481612, 17573267)

sum(all_smarter)
sum(all_smart_seq)
sum(all_tang)

length(all_zeros_for_tang)





