# Load libraries
library(ggplot2)
library(RColorBrewer)

#Input your data
d=read.csv("Table1.csv")
#Select a color palette
coul <- brewer.pal(4, "Paired")
#Plot data
p<-ggplot(aes(x = Rootstock, y = Bark), data = d) + stat_summary(fun.y = "mean", geom = "bar",fill=coul)
#Add error bars
p<-p+stat_summary(data=bb, fun.y=mean, fun.ymin=min, fun.ymax=max, geom="errorbar", width=0.5)
#Change theme and label size
p<-p+theme_classic(base_size=20)
#Edit x asis
p<-p+ylab("Bark percentage")

# Examples
B<-read.table("Book1.txt",header=T,sep="\t")
pie<-ggplot(B, aes(x=Isolate, y=Amount, fill=Category)) + geom_bar(width = 1, size = 1, color = "white", stat = "identity")
pie<-ggplot(B, aes(x=DEGs, y=Amount, fill=Category)) + geom_bar(width = 1, size = 1, color = "white", stat = "identity") + facet_grid(~Order+Host+Sample+CAZY)
pie = pie + scale_fill_manual(values=c("#55DDE0", "#33658A", "#2F4858", "#F6AE2D", "#F26419", "#999999"))
pie = pie + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
pie = pie +  coord_flip()
