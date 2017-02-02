options(stringsAsFactors = FALSE)

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

heatmap_data <- read.delim("../Genes.tsv")
heatmap_data <- as.data.frame(heatmap_data)
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
  
  heatmap_data$Sample.ID[x] <- paste(heatmap_data$Sample.ID[x], heatmap_data$MLST[x], sep = "  ")
  
  for(y in 1:ncol(heatmap_data)){
    if(heatmap_data[x,y] == "I"){heatmap_data[x,y] <- 1}
    if(heatmap_data[x,y] == "P"){heatmap_data[x,y] <- 2}
    if(heatmap_data[x,y] == "*"){heatmap_data[x,y] <- 3}
    if(heatmap_data[x,y] == "-"){heatmap_data[x,y] <- 4}
  }
  heatmap_data[x,] <- as.numeric(heatmap_data[x,])
}

colors <- colorRampPalette(c("yellow", "green", "chartreuse", "red"))(n = 4)
categorize = c(seq(1,1,length=1), 
               seq(2,2,length=1),  
               seq(3,3,length=1),
               seq(4,4, length=1))

heatmap_data <- heatmap_data[, -2]
rownames(heatmap_data) <- heatmap_data$Sample.ID
heatmap_data <- heatmap_data[, -1]

heatmap.2(heatmap_data, col = colors, breaks = categorize, Colv = "NA", Rowv = "NA", dendrogram = "none")

input <- read.delim("../DC_LMC.tab")

formatted_input <- input[, c(1:3) * -1]
rownames(formatted_input) <- paste(input$X.CHROM, input$POS, input$REF, sep = "_")
formatted_input <- t(formatted_input)

distance_calc <- as.data.frame(matrix(ncol = nrow(formatted_input), nrow = nrow(formatted_input)))
rownames(distance_calc) <- rownames(formatted_input)
colnames(distance_calc) <- rownames(formatted_input)

for(x in rownames(distance_calc)){
  for(y in colnames(distance_calc)){
    distance_calc[x,y] <- sum(formatted_input[x,] != formatted_input[y,])
  }
}

distance <- as.dist(distance_calc)
clustered <- hclust(distance)

#pdf(file = "test_out.pdf")
#plot(clustered)
#dev.off()





