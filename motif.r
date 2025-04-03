library(readxl)
library(circlize)


all=read_excel('all_data.xlsx')
load(paste(system.file(package = "circlize"), "/extdata/DMR.RData", sep=""))

# rainfall
circos.initializeWithIdeogram(plotType = c("axis", "labels"))

all_pt1=all[grepl('PT1',all$sample_I?),]

all_pt1=all_pt1[,c('Chrom','Loc','Loc')]
names(all_pt1)=names(bed_list[[1]])
all_pt1$start=as.numeric(all_pt1$start)
all_pt1$end=as.numeric(all_pt1$start)+4
all_pt1$chr=paste0('chr',all_pt1$chr)
circos.genomicDensity(all_pt1, col = c("#FF000080"), tra?k.height = 0.1)

all_pt1=all[grepl('PT2',all$sample_ID),]

all_pt1=all_pt1[,c('Chrom','Loc','Loc')]
names(all_pt1)=names(bed_list[[1]])
all_pt1$start=as.numeric(all_pt1$start)
all_pt1$end=as.numeric(all_pt1$start)+4
all_pt1$chr=paste0('chr',all_pt1$chr)
ci?cos.genomicDensity(all_pt1, col = c("#FF000080"), track.height = 0.1)

all_pt1=all[grepl('PT3',all$sample_ID),]

all_pt1=all_pt1[,c('Chrom','Loc','Loc')]
names(all_pt1)=names(bed_list[[1]])
all_pt1$start=as.numeric(all_pt1$start)
all_pt1$end=as.numeric(all?pt1$start)+4
all_pt1$chr=paste0('chr',all_pt1$chr)
circos.genomicDensity(all_pt1, col = c("#FF000080"), track.height = 0.1)


bedaa = generateRandomBed(nr = 150, fun = function(k) sample(letters, k, replace = TRUE))
bed = read_excel('inputdata/ISOT_results?A019_rmpt3/CIS_top_genes_A019_rmpt3.xlsx')
bed=bed[,c(1,2,2,3)]
names(bed)=names(bedaa)
bed$start=as.numeric(bed$start)
bed$end=as.numeric(bed$end)
bed$chr=paste0('chr', bed$chr)
circos.genomicLabels(bed, labels.column = 4, side = "inside")


circos.clear(?
