# Load libraries
library("gplots")
library(RColorBrewer)

# Prepare the data
T1<-read.table("Table1.txt",header=T,sep="\t")
T2<-read.table("Table2.txt",header=T,sep="\t")
T3<-merge(T1,T2, by.x="ID",by.y="ID",all.x=TRUE,all.y=TRUE)
write.table(T3,"Table3.txt",sep="\t",na="",quote=F)
T3<-read.table("Table3.txt",header=T,sep="\t")

# Define color palette
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

# Replace NA with 0
T3[is.na(T3)] <- 0
View(`T3`)

# Assign labels in column 1 to "rnames"
rnames <- T3[,1]
# Transform column 2-5 into a matrix
mat_data <- data.matrix(T3[,2:ncol(T3)])
# Assign row names
rownames(mat_data) <- rnames

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(-1,0,length=100),  # for red
  seq(0.01,0.8,length=100),           # for yellow
  seq(0.81,1,length=100))             # for green

# creates a 5 x 5 inch image
png("heatmaps_in_r.png",    # create PNG for the heat map
  width = 5*300,        # 5 x 300 pixels
  height = 5*300,
  res = 300,            # 300 pixels per inch
  pointsize = 8)        # smaller font size

heatmap.2(mat_data,
  cellnote = mat_data,  # same data set for cell labels
  main = "Correlation", # heat map title
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier
  breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="row",     # only draw a row dendrogram
  Colv="NA")            # turn off column clustering

dev.off()               # close the PNG device

# Examples

# Creates image
png("heatmaps_in_r_all.png",    # create PNG for the heat map
width = 10000,        # 5 x 300 pixels
height = 10000,
res = 300,            # 300 pixels per inch
pointsize = 8)        # smaller font size

#Plot heatmap
heatmap.2(mat_data, trace="none", key=TRUE, Colv=F, dendrogram="row", cexRow=2, cexCol=1, lwid=c(2, 6), lhei = c(2,12), density.info = "density", margins = c(1, 10), col=r)
dev.off()

heatmap.2(mat_data, trace="none", key=TRUE, Colv=F, dendrogram="row", cexRow=3, cexCol=1, lwid=c(2, 8), lhei = c(2,12), density.info = "density", margins = c(1, 40), col=my_palette)
dev.off()

heatmap.2(mat_data, trace="none", key=TRUE, keysize=0.05, Colv=F, dendrogram="row", cexRow=0.6, cexCol=1, lmat=rbind( c(0, 3, 4), c(2,1,1 ) ), lwid=c(1.5, 4, 2 ), lhei = c(0.5,4),  density.info = "density", margins = c(6, 6), col=r)
dev.off()

heatmap.2(mat_data, cellnote = mat_data, notecol="black", density.info="none",  # turns off density plot inside color legend
  trace="none", margins =c(12,9), col=my_palette, dendrogram="row", Colv="NA")
dev.off()

png("S7.png", width=1000, height=1000)
pheatmap(mat_data, trace="none", key=T, keysize=5, Colv=F, dendrogram="row", cexRow=0.6, cexCol=1, lhei = c(0.05,5), density.info = "none", margins = c(7, 7), col=r)
dev.off()

colsep=c(2:4) #Separate columns
rowsep=(1:62) #Separate rows
sepcolor="white"
sepwidth=c(0.8,0.8)
Rowv=F,Colv=F
scale="none"         

#Alternative palette
my_palette <- colorRampPalette(c("green", "green3", "black", "red3","red"))(n = 100)
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)


