options(stringsAsFactors = FALSE)

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

heatmap_data <- read.delim("Genes.tsv")
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
for(x in 1:nrow(heatmap_data)){
  if(nchar(heatmap_data$Sample.ID[x]) < max_char){
    heatmap_data$Sample.ID[x] <- paste(heatmap_data$Sample.ID[x], rep.int(" ", (max_char - nchar(heatmap_data$Sample.ID[x]))), sep = "")
  }
  
  heatmap_data$Sample.ID[x] <- paste(heatmap_data$MLST[x], heatmap_data$Sample.ID[x] , sep = "     ")
  
  for(y in 3:ncol(heatmap_data)){
    if(heatmap_data[x,y] == "I"){heatmap_data[x,y] <- 1.5}
    if(heatmap_data[x,y] == "P"){heatmap_data[x,y] <- 2.5}
    if(heatmap_data[x,y] == "*"){heatmap_data[x,y] <- 3.5}
    if(heatmap_data[x,y] == "-"){heatmap_data[x,y] <- 4.5}
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

input <- read.delim("DC_LMC.tab")

formatted_input <- input[, c(1:3) * -1]
rownames(formatted_input) <- paste(input$X.CHROM, input$POS, input$REF, sep = "_")
formatted_input <- t(formatted_input)
formatted_input <- formatted_input[grep("ref", rownames(formatted_input)) * -1, ]

distance_calc <- as.data.frame(matrix(ncol = nrow(formatted_input), nrow = nrow(formatted_input)))
rownames(distance_calc) <- rownames(formatted_input)
colnames(distance_calc) <- rownames(formatted_input)

for(x in rownames(distance_calc)){
  for(y in colnames(distance_calc)){
    distance_calc[x,y] <- sum(formatted_input[x,] != formatted_input[y,])
  }
}

for(x in 1:nrow(distance_calc)){
  st <- rownames(distance_calc)[x]
  st2 <- gsub(".cat.fasta", "", st)
  st3 <- gsub(".", "-", st2, fixed = TRUE)
  
  new_str <- rownames(heatmap_data)[grep(st3, rownames(heatmap_data))]
  rownames(distance_calc)[x] <- new_str
  colnames(distance_calc)[x] <- new_str
}

distance <- as.dist(distance_calc)
clustered <- hclust(distance, method = "complete")

pdf("shitty_file.pdf", width = nrow(heatmap_data) * 0.5, height = ncol(heatmap_data) * 0.5 * 2)
heatmap.2(t(as.matrix(heatmap_data)), 
          col = colors, 
          distfun = NA,
          hclustfun = NA,
          breaks = categorize, 
          Colv = as.dendrogram(clustered), 
          Rowv = FALSE, 
          dendrogram = "column", 
          key = FALSE, tracecol = NA, margins=c(8,14), lwid = c(1,12))
dev.off()

#pdf(file = "test_out.pdf")
#plot(clustered)
#dev.off()





