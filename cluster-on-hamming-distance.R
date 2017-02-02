input <- read.delim("xd.tab")

formatted_input <- input[, c(1:3) * -1]
rownames(formatted_input) <- paste(input$X.CHROM, input$POS, input$REF, sep = "_")
formatted_input <- t(formatted_input)

distance_calc <- as.data.frame(matrix(ncol = nrow(formatted_input), nrow = nrow(formatted_input)))
rownames(distance_calc) <- rownames(formatted_input)
colnames(distance_calc) <- rownames(formatted_input)

for(x in rownames(distance_calc)){
  for(y in colnames(distance_calc)){
    distance_calc[x,y] <- sum[formatted_input[x,] != formatted_input[y,]]
  }
}

distance <- as.dist(distance_calc)
clustered <- hclust(distance)

pdf(file = "test_out.pdf")
plot(clustered)
dev.off()