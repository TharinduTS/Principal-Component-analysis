# Principal-Component-analysis

I am using pre filtered vcfs for admix plots for this analysis.

working in

/for_PCA/

copied all vcf files to a sub directory here

then bgzip them

```bash
module load StdEnv/2020  intel/2020.1.217 bcftools/1.11
for i in *.vcf; do bgzip -c ${i}> ${i}.gz ; done
```
index them
```bash 
for i in *.vcf.gz; do tabix -p vcf ${i} ; done
```
get a list file names seperated by a space using
```bash
ls *.vcf.gz | tr "\n" " "
```

create vcf files for different subgenomes using those lists
run this inside the folder with the vcf files
```
#!/bin/sh
#SBATCH --job-name=abba
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=32gb
#SBATCH --output=abba.%J.out
#SBATCH --error=abba.%J.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

module load StdEnv/2020
module load vcftools/0.1.16
module load intel/2020.1.217 bcftools/1.11


bcftools concat DB_chr1L_positions_excluded.vcf.recode.vcf.gz DB_chr1S_positions_excluded.vcf.recode.vcf.gz DB_chr2L_positions_excluded.vcf.recode.vcf.gz DB_chr2S_positions_excluded.vcf.recode.vcf.gz DB_chr3L_positions_excluded.vcf.recode.vcf.gz DB_chr3S_positions_excluded.vcf.recode.vcf.gz DB_chr4L_positions_excluded.vcf.recode.vcf.gz DB_chr4S_positions_excluded.vcf.recode.vcf.gz DB_chr5L_positions_excluded.vcf.recode.vcf.gz DB_chr5S_positions_excluded.vcf.recode.vcf.gz DB_chr6L_positions_excluded.vcf.recode.vcf.gz DB_chr6S_positions_excluded.vcf.recode.vcf.gz DB_chr7L_positions_excluded.vcf.recode.vcf.gz DB_chr7S_positions_excluded.vcf.recode.vcf.gz DB_chr8L_positions_excluded.vcf.recode.vcf.gz DB_chr8S_positions_excluded.vcf.recode.vcf.gz DB_chr9L_positions_excluded.vcf.recode.vcf.gz DB_chr9S_positions_excluded.vcf.recode.vcf.gz -o ../test_combine/l_and_s.vcf.gz

bcftools concat DB_chr1L_positions_excluded.vcf.recode.vcf.gz DB_chr2L_positions_excluded.vcf.recode.vcf.gz DB_chr3L_positions_excluded.vcf.recode.vcf.gz DB_chr4L_positions_excluded.vcf.recode.vcf.gz DB_chr5L_positions_excluded.vcf.recode.vcf.gz DB_chr6L_positions_excluded.vcf.recode.vcf.gz DB_chr7L_positions_excluded.vcf.recode.vcf.gz DB_chr8L_positions_excluded.vcf.recode.vcf.gz DB_chr9L_positions_excluded.vcf.recode.vcf.gz -o ../test_combine/l_only.vcf.gz

bcftools concat DB_chr1S_positions_excluded.vcf.recode.vcf.gz DB_chr2S_positions_excluded.vcf.recode.vcf.gz DB_chr3S_positions_excluded.vcf.recode.vcf.gz DB_chr4S_positions_excluded.vcf.recode.vcf.gz DB_chr5S_positions_excluded.vcf.recode.vcf.gz DB_chr6S_positions_excluded.vcf.recode.vcf.gz DB_chr7S_positions_excluded.vcf.recode.vcf.gz DB_chr8S_positions_excluded.vcf.recode.vcf.gz DB_chr9S_positions_excluded.vcf.recode.vcf.gz -o ../test_combine/s_only.vcf.gz
```
check sample list if needed
```
bcftools query -l l_and_s.vcf.gz
```
remove extra samples
```bash
vcftools --remove-indv BJE2897_Lendu_TTCCTGGA_cuttrim_sorted.bam --remove-indv BJE3252_Cameroon_TATCGGGA_cuttrim_sorted.bam --remove-indv RT5_Botsw_GGATTGGT_cuttrim_sorted.bam --remove-indv amnh17260_Nigeria_GTGAGGGT_cuttrim_sorted.bam --gzvcf l_and_s.vcf.gz --out ../extra_samples_removed_VCFs/samples_removed_l_and_s.vcf.gz --recode

vcftools --remove-indv BJE2897_Lendu_TTCCTGGA_cuttrim_sorted.bam --remove-indv BJE3252_Cameroon_TATCGGGA_cuttrim_sorted.bam --remove-indv RT5_Botsw_GGATTGGT_cuttrim_sorted.bam --remove-indv amnh17260_Nigeria_GTGAGGGT_cuttrim_sorted.bam --gzvcf l_only.vcf.gz --out ../extra_samples_removed_VCFs/samples_removed_l_only.vcf.gz --recode

vcftools --remove-indv BJE2897_Lendu_TTCCTGGA_cuttrim_sorted.bam --remove-indv BJE3252_Cameroon_TATCGGGA_cuttrim_sorted.bam --remove-indv RT5_Botsw_GGATTGGT_cuttrim_sorted.bam --remove-indv amnh17260_Nigeria_GTGAGGGT_cuttrim_sorted.bam --gzvcf s_only.vcf.gz --out ../extra_samples_removed_VCFs/samples_removed_s_only.vcf.gz --recode
```

get a list of samples and populations in the same order as in vcf file

do_PCA.R
```R

#https://github.com/zhengxwen/SNPRelate/issues/13
#http://corearray.sourceforge.net/tutorials/SNPRelate/
#https://github.com/zhengxwen/SNPRelate/wiki/Preparing-Data

# ***************   These are the things you should edit  **********************

input_file_name<-"samples_removed_l_only.vcf.gz.recode.vcf"
plot_title_to_use<-"L subgenome only"
saving_file_name<-"l_only_selected_samples.pdf"
sample_data_file_name<-"selected_samples_list.txt"

# if you added or changed colors in excel sheet with new colors, you have to edit line for colors below  as well
# ******************************************************************************
# ******install all the packages needed if you are running it on computecanada*********

# install.packages("devtools")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("gdsfmt")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("SNPRelate")
# install.packages("plyr")
# install.packages("ggplot2")
# install.packages("ggrepel")
# ************

library("devtools")
library(gdsfmt)
library(SNPRelate)

# set working directory to current script directory 
# **uncomment this if you use this in local computer. KEEP COMMENTED OUT IF YOU ARE ON COMPUTECANADA) ***
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#enter input vcf filename here
vcf.fn <- input_file_name


# ************* if you retry this, you may have to restart r session before trying ************

snpgdsVCF2GDS(vcf.fn, "test.gds", method="copy.num.of.ref",ignore.chr.prefix = "chr")

# *********************************************************************************************

#snpgdsSummary("test.gds")

#genofile = openfn.gds("test.gds", readonly=FALSE)
genofile = snpgdsOpen("test.gds", readonly=FALSE)

#samp.annot<-data.frame(pop.group = c("brunnescens","hecki","hecki","hecki","hecki","hecki","hecki","maura","maura","maura","maura","maura","maura","Borneo","Sumatra","Malay","Sumatra","Borneo","Borneo","Borneo","siberu","nigra","nigra","nigrescens","nigra","nigrescens","ochreata","ochreata","ochreata","togeanus","togeanus","tonkeana","tonkeana","tonkeana","tonkeana","tonkeana","tonkeana","tonkeana","tonkeana","tonkeana"))

#trying to read sample names and populations from a tsv file instead typing all here
#name the tsv file as samples.list with sample name in the first column , population in second column and colour in third column
sample_data<-read.table(sample_data_file_name)

#create a sample name array in the same order
sample_name_order<-as.array(sample_data$V1)
population_order<-as.array(sample_data$V2)
color_order<-as.array(sample_data$V3)

# get simpler names  for those complex sample names
library(plyr)
#first renaming the sample names which does not follow the specific format manually**DO NOT USE"_" HERE AS SPLIT WILL USE THIS IN NEXT LINE

#copy sample name list to edit
full_sample_names<-sample_name_order

#rename the samples that does not follow that format manually
simple_sample_list_changing<-mapvalues(full_sample_names, from = c("Vred_8_Vred_GCTGTGGA_cuttrim_sorted.bam","JM_no_label1_Draken_CCACGT_cuttrim_sorted.bam","JM_no_label2_Draken_TTCAGA_cuttrim_sorted.bam","946_Draken_TCGTT_cuttrim_sorted.bam","993_Draken_GGTTGT_cuttrim_sorted.bam","2014_Inhaca_10_Inhaca_ATATGT_cuttrim_sorted.bam","2014_Inhaca_150_Inhaca_ATCGTA_cuttrim_sorted.bam","2014_Inhaca_152_Inhaca_CATCGT_cuttrim_sorted.bam","2014_Inhaca_24_Inhaca_CGCGGT_cuttrim_sorted.bam","2014_Inhaca_38_Inhaca_CTATTA_cuttrim_sorted.bam","2014_Inhaca_52_Inhaca_GCCAGT_cuttrim_sorted.bam","2014_Inhaca_65_Inhaca_GGAAGA_cuttrim_sorted.bam"), 
                                       to = c("Vred8","JM1","JM2","Draken946","Draken993","Inhaca10","Inhaca150","Inhaca152","Inhaca24","Inhaca38","Inhaca52","Inhaca65"))

#convert facter list into chrs to rename
s_list_chr<-as.character(simple_sample_list_changing)

#extract sample IDs from names
shortened_sample_list<-sapply(strsplit(s_list_chr,split = "_"),`[`, 1)

#samp.annot<-data.frame(pop.group = c("brunnescens","hecki","hecki","hecki","hecki","hecki","hecki","maura","maura","maura","maura","maura","maura","Borneo","Sumatra","Malay","Sumatra","Borneo","Borneo","Borneo","siberu","nigra","nigra","nigrescens","nigra","nigrescens","ochreata","ochreata","ochreata","togeanus","togeanus","tonkeana","tonkeana","tonkeana","tonkeana","tonkeana","tonkeana","tonkeana","tonkeana","tonkeana"))
samp.annot<-data.frame(pop.group = population_order)



#and then load
add.gdsn(genofile, "sample.annot", samp.annot)

snpgdsSummary("test.gds")

# LD prinnung
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,  method = c("composite"),missing.rate=0, verbose = TRUE)
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))


# or without LD prunning
#pca <- snpgdsPCA(genofile, num.thread=2)
#pc.percent <- pca$varprop*100
#head(round(pc.percent, 2))


tab <- data.frame(sample.id = pca$sample.id,
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
head(tab)

#plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")
#text(tab$EV2, tab$EV1,labels=tab$sample.id, cex= 0.4)

library(ggplot2)
#ggplot(...)+...+ theme(axis.text.x = element_text(angle=60, hjust=1))
#devtools::install_github("slowkow/ggrepel")
library(ggrepel)

pdf(saving_file_name,w=8, h=8, version="1.4", bg="transparent")
tab$Species <- population_order
tab$samp.color <- color_order
tab$samp.fieldid <- shortened_sample_list

# get unique colours and pop names as labels to use in next line
pop_colors<-unique(color_order)

# you will have to use a manual array if you are going to apply same color to more than one pop
pop_labels<-unique(population_order)

# to apply correct colors for the values
#getting same colors for the values
#colors_to_apply<-pop_colors

#creating the array for values
#col_vals<-as.array(paste("'",pop_colors,"'","=","'",colors_to_apply,"'",sep=""))

d<-ggplot(data=tab, aes(x=-EV1,y=-EV2, label = samp.fieldid, color = samp.color)) +
  # label axis 
 	labs(x=expression("-Eigenvector 1"), y=expression("-Eigenvector 2"),title =plot_title_to_use,cex=1) +
  # legend details
  
  # ************
  
  # if you add more colors in excel sheet, you have to add them here as well
 	scale_colour_manual(name="Population", values = c("pink"="pink","forestgreen"="forestgreen","orange"="orange","skyblue"="skyblue","green"="green","blue"="blue","lightgreen"="lightgreen","yellow"="yellow","purple"="purple","red"="red","lightgray"="lightgray","darkgray"="darkgray","gray"="gray","brown"="brown"),breaks=pop_colors,labels=pop_labels)+
 	
  # ************
  
  # add points and fieldID labels
	geom_text_repel(aes(-EV1,-EV2, label=(samp.fieldid)), size=6, point.padding = unit(0.5, "lines"),max.overlaps = 0) + geom_point(size=6) + 
  #geom_label(aes(fill = factor(samp.color)), colour = "white", fontface = "bold") +
  # change to cleaner theme
	theme_classic(base_size = 16) +
	# make it clean
	theme_bw()+ theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  
  #this is for theme inside plot
  # theme(legend.key.size = unit(0.15, 'cm'), #change legend key size
  #      legend.key.height = unit(0.15, 'cm'), #change legend key height
  #      legend.key.width = unit(0.15, 'cm'), #change legend key width
  #      legend.title = element_text(size=11), #change legend title font size
  #      legend.text = element_text(size=11,face = "italic"), #change legend text font size
  #      legend.position = c(0.2, .15))+
  
  
  
  
  
  
	# italicize species names
	#theme(legend.text = element_text(face="italic"))+hel
  # make the text bigger
  #theme(text = element_text(size=20)) +
	# move the legend
	#theme(legend.position = c(.18, .15)) +
  # add space around axis
    theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt"),size = 11),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt"),size = 11)) +
  # remove boxes around legend symbols
  theme(legend.key = element_blank())+
  #these are for annotations. not used here
	annotate(geom = "text", x = 0.03, y = .23, label = "", color = "black", angle = 75, size=7)+
	annotate(geom = "text", x = -.12, y = -.10, label = "", color = "black", angle = 0, size=7)+
	annotate(geom = "text", x = -.25, y = .07, label = "", color = "black", angle = -15, size=7)
d
dev.off()




closefn.gds(genofile)
```
