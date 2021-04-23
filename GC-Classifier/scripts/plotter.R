# library ggplot2
library(ggplot2)

# importing .tsv file as an argument
tsv <- commandArgs(trailingOnly = TRUE)

# reading the .tsv file and assigning its data into a variable (d)
d <- read.table(tsv, sep ="\t", header= TRUE)

# creating a png
png(file = "bin/result.png", width=800, height=400)

# defining the axis (y-axis and x-axis) and plotting the data
plot <- ggplot(data = d) +geom_line(aes(y = GC1, x = Pos, group=1), color="darkred", size = 0.65) +
			geom_line(aes(y = GC2, x = Pos, group=1), color="steelblue", size = 0.65)

# labelling x-axis, y-axis and plot title
plot + labs(title="GC Content Distribution",
        x ="Position (in Mb)", y = "GC Percent (%)")

# saving the plot
dev.off()
