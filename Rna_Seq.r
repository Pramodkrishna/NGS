source("http://bioconductor.org/biocLite.R") # downloads the latest bioconductor script/repos
biocLite("limma")
biocLite("edgeR")
biocLite("org.Mm.eg.db")
biocLite("Glimma")

library(edgeR)
library(limma)
library(ggplot2)
library(Glimma)
library(RColorBrewer)
library(org.Mm.eg.db)
sampleinfo <- read.delim("SampleInfo.txt")

View(sampleinfo)

seqdata <- read.delim("GSE60450_Lactation-GenewiseCounts.txt", stringsAsFactors = FALSE)
head(seqdata)
View(seqdata)
dim(seqdata)

countdata <- seqdata[,-(1:2)]

# Store EntrezGeneID as rownames
rownames(countdata) <- seqdata[,1]

View(countdata)
head(countdata)
colnames(countdata)


# using substr, you extract the characters starting at position 1 and stopping at position 7 of the colnames
colnames(countdata) <- substr(colnames(countdata), 1, 7)
View(countdata)

table(colnames(countdata)==sampleinfo$SampleName)
#checking whether the names match from the seq data file and from the sample file

myCPM <- cpm(countdata)
head(countdata)
thresh <- myCPM > 0.5
View(thresh)
table(rowSums(thresh))
keep <- rowSums(thresh) >= 2
counts.keep <- countdata[keep,]
summary(keep)


plot(myCPM[,1],countdata[,1])
# Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
plot(myCPM[,1],countdata[,1],ylim=c(0,50),xlim=c(0,3))
# Add a vertical line at 0.5 CPM
abline(v=0.5)

dgeObj <- DGEList(counts.keep)
# have a look at dgeObj
dgeObj
# See what slots are stored in dgeObj
names(dgeObj)
# Library size information is stored in the samples slot


#First, we can check how many reads we have for each sample in the dgeObj
dgeObj$samples$lib.size

barplot(dgeObj$samples$lib.size, names=colnames(dgeObj), las=2)
# Add a title to the plot
title("Barplot of library sizes")

logcounts <- cpm(dgeObj,log=TRUE)

# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

plotMDS(dgeObj)

# We specify the option to let us plot two plots side-by-sde
par(mfrow=c(1,2))
# Let's set up colour schemes for CellType
# How many cell types and in what order are they stored?
levels(sampleinfo$CellType)
## Let's choose purple for basal and orange for luminal
col.cell <- c("purple","orange")[sampleinfo$CellType]
data.frame(sampleinfo$CellType,col.cell)

# Redo the MDS with cell type colouring
plotMDS(dgeObj,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend("topleft",fill=c("purple","orange"),legend=levels(sampleinfo$CellType))
# Add a title
title("Cell type")

# Similarly for status
levels(sampleinfo$Status)
col.status <- c("blue","red","dark green")[sampleinfo$Status]
col.status
plotMDS(dgeObj,col=col.status)
legend("topleft",fill=c("blue","red","dark green"),legend=levels(sampleinfo$Status),cex=0.8)
title("Status")