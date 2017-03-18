args <- commandArgs(trailingOnly = TRUE)
gene_data <- args[1]
snp_data <- args[2]
output_file <- args[3]
options(stringsAsFactors = FALSE)

#make sure output is svg
if(tolower(substr(output, nchar(output) - 2, nchar(output))) != "pdf"){
  warning("not a pdf, prepare for arbitrary file corruption")
}

if(!file.exists(gene_data)){
  stop("missing gene data, terminated")
}

if(!file.exists(snp_data)){
  warning("no snp data, no clustering")
  snp_data <- NA
}


list.of.packages <- c("gplots", "RColorBrewer")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.utstat.utoronto.ca/')

#load packages
library(gplots)
library(RColorBrewer)


if(!is.na(snp_data)){
  input <- read.delim(snp_data)
}else{ input <- NA}

heatmap_data <- read.delim(gene_data)
colnames(heatmap_data)[1] <- "Sample.ID"
heatmap_data$Sample.ID <- as.character(heatmap_data$Sample.ID)

max_char <- 0
for(x in heatmap_data$Sample.ID){
  if(nchar(x) > max_char){
    max_char <- nchar(x)
  }
}

#I = yellow = 1
#P = green = 2
#* = other green = 3
#- = red = 4
new_max_char <- max_char

for(x in 1:nrow(heatmap_data)){
  if(nchar(heatmap_data$Sample.ID[x]) < max_char){
    some_number <- (max_char - nchar(heatmap_data$Sample.ID[x]))
    pos <- max(gregexpr("-", heatmap_data$Sample.ID[x], fixed=TRUE)[[1]])
    
    first_part <- substr(heatmap_data$Sample.ID[x], 0, pos)
    filler <- paste(rep.int(0, some_number), collapse = "")
    last_part <- substr(heatmap_data$Sample.ID[x], pos + 1, nchar(heatmap_data$Sample.ID[x]))
    
    heatmap_data$Sample.ID[x] <- paste(first_part, filler, last_part, sep = "")
  }
  
  
  heatmap_data$Sample.ID[x] <- paste(heatmap_data$MLST[x], heatmap_data$Sample.ID[x] , sep = "     ")
  if(nchar(heatmap_data$Sample.ID[x]) > new_max_char){
    new_max_char <- nchar(heatmap_data$Sample.ID[x])
  }
  
  for(y in 3:ncol(heatmap_data)){
    if(as.character(heatmap_data[x,y]) == "I"){
      heatmap_data[x,y] <- 1.5
    }
    else if(as.character(heatmap_data[x,y]) == "P"){
      heatmap_data[x,y] <- 2.5
    }
    else if(as.character(heatmap_data[x,y]) == "S"){
      heatmap_data[x,y] <- 3.5
    }
    else if(as.character(heatmap_data[x,y]) == "-"){
      heatmap_data[x,y] <- 4.5
    }
  }
}


colors <- colorRampPalette(c("yellow", "green", "chartreuse", "red"))(n = 4)
categorize = c(seq(1,2,length=1), 
               seq(2,3,length=1),  
               seq(3,4,length=1),
               seq(4,5, length=1),
               seq(5,5, length=1))

heatmap_data <- heatmap_data[, -2]
rownames(heatmap_data) <- heatmap_data$Sample.ID
heatmap_data <- heatmap_data[, -1]

for(y in 1:ncol(heatmap_data)){ heatmap_data[, y] <- as.numeric(heatmap_data[,y])}

dendro <- NA

if(!is.na(input)){
  
  formatted_input <- input[, c(1:3) * -1]
  rownames(formatted_input) <- paste(input$X.CHROM, input$POS, input$REF, sep = "_")
  formatted_input <- t(formatted_input)
  formatted_input <- formatted_input[grep("ref", rownames(formatted_input)) * -1, ]
  
  distance_calc <- as.data.frame(matrix(ncol = nrow(formatted_input), nrow = nrow(formatted_input)))
  rownames(distance_calc) <- rownames(formatted_input)
  colnames(distance_calc) <- rownames(formatted_input)
  
  for(x in rownames(distance_calc)){
    for(y in colnames(distance_calc)){
      if(!is.na(distance_calc[x,y])){ next }
      value <- sum(formatted_input[x,] != formatted_input[y,])
      distance_calc[x,y] <- distance_calc[y,x] <- value
    }
  }
  
  distance_calc <- readRDS("distance_matrix.RDS")
  for(x in 1:nrow(distance_calc)){
    st <- rownames(distance_calc)[x]
    st2 <- gsub(".cat.fasta", "", st)
    st3 <- gsub(".", "-", st2, fixed = TRUE)
    st3 <- gsub("X", "", st3, fixed = TRUE)
    
    if(nchar(st3) < max_char){
      some_number <- (max_char - nchar(st3))
      pos <- max(gregexpr("-", st3, fixed=TRUE)[[1]])
      
      first_part <- substr(st3, 0, pos)
      filler <- paste(rep.int(0, some_number), collapse = "")
      last_part <- substr(st3, pos + 1, nchar(st3))
      
      st3 <- paste(first_part, filler, last_part, sep = "")
    }
    if(st3 ==  "95-00000000012"){ st3 <- "95-0012_GATCAG"}
    if(st3 == "95-00000000151"){ st3 <- "95-0151_ACTTGA"}
    new_str <- rownames(heatmap_data)[grep(st3, rownames(heatmap_data))]

    rownames(distance_calc)[x] <- new_str
    colnames(distance_calc)[x] <- new_str
  }
  
  if(length(intersect(rownames(heatmap_data), rownames(distance_calc))) != length(rownames(heatmap_data))){
    warning("snp and heatmap data rownames do not match, prep for arbitrary garbo")
  }
  distance <- as.dist(distance_calc)
  clustered <- hclust(distance, method = "complete")
  heatmap_data <- t(heatmap_data)
  heatmap_data <- heatmap_data[, rownames(distance_calc)]
  dendro <- as.dendrogram(clustered)
}

svg(output_file, width = ncol(heatmap_data) * 0.5, height = nrow(heatmap_data) * 0.5)
heatmap.2(as.matrix(heatmap_data), 
          col = colors, 
          breaks = categorize, 
          Colv = dendro, 
          Rowv = NULL,
          dendrogram = "column", 
          key = FALSE, tracecol = NA, margins=c(8,14), lwid = c(1,6))
dev.off()





