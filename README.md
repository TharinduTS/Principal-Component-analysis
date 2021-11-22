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

then you can make them seperated by commas and add invited commas using word find and replace

find
```text
^p
```
replace
```txt
","
```

****** edit the very begining and end before using


