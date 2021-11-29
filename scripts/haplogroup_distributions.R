#library(ggplot2)
#library(ggpubr)
#theme_set(theme_pubr())

in_file <- snakemake@input[[1]]
hap_file <- snakemake@input[[2]]
out_file_dist <- snakemake@output[[1]]
out_file_freq <- snakemake@output[[2]]


data <- read.delim(in_file,sep="\t", header=TRUE)

haplogroups_to_plot <- read.delim(hap_file,sep="\t", header=FALSE)$V1
haplogroups_to_plot

data_pop_haplogroup <- data[,c("population","Major.Haplogroup")]

data_pop_haplogroup$to_plot <- NA

for(i in 1:length(haplogroups_to_plot)){
    hap_i <- as.character(haplogroups_to_plot[i])
    for(j in 1:length(data_pop_haplogroup[,1])){
        hap_sample_j <- data_pop_haplogroup[j,2]
        if(substr(hap_sample_j,1,nchar(hap_i)) == hap_i){
            data_pop_haplogroup[j,]$to_plot <- hap_i
        }
    }
}

#head(data_pop_haplogroup)

data_pop_haplogroup$to_plot <- as.factor(data_pop_haplogroup$to_plot)

data_plot <- data_pop_haplogroup[,c("population","to_plot")]
summary(data_plot)

pdf(out_file_dist)
 
# Make a stacked barplot
cols = c("#A4A4A4","#222222","#A945FF","#1D72F5","#DF0101","#77CE61","#000000","#FF9326","#0089B2","#FDF060","#FFA6B2","#BFF217","#60D5FD","#CC1577","#F2B950","#7FB21D","#EC496F","#326397","#B26314","#027368","#610B5E","#F3C300","#875692","#F38400","#1D72F5","khaki2","darkturquoise","#bfbfbf","#E31A1C","#F2F3F4","#FF9326","#000000")
barplot(t(table(data_plot)), col=cols, border="white", xlab="",las=2)
legend("topright",fill=cols,rownames(t(table(data_plot))))

dev.off()


pdf(out_file_freq)

# Make a frequency plot 
#data_for_freq <- as.data.frame(t(table(data_plot)))

#ggballoonplot(test, fill = "value")+
#  scale_fill_viridis_c(option = "C")

dev.off()