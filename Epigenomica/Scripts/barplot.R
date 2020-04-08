library(ggplot2)
df <- read.table(file= "./bed_E2_filtrado_0.7/by_chr/02_summary_segments.txt")

pdf("./bed_E2_filtrado_0.7/barplot_average_segments.pdf") 
ggplot(data = df, mapping = aes(x = reorder(V2, V1), V1)) + 
  geom_bar(stat = "identity", fill= "salmon")  + coord_flip()+ 
  xlab("")+ ggtitle("Average number of segments by chromosome; State 2 f = 0.7")+
  ylab("Number of 200 bp segments")+
  theme(
    plot.title = element_text(size = 14, face = "bold"), 
    axis.text = element_text(size = 12)
  )
dev.off() 
