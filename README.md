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
#SBATCH --time=12:00:00
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

vcf-merge DB_chr1L_positions_excluded.vcf.recode.vcf.gz DB_chr1S_positions_excluded.vcf.recode.vcf.gz DB_chr2L_positions_excluded.vcf.recode.vcf.gz DB_chr2S_positions_excluded.vcf.recode.vcf.gz DB_chr3L_positions_excluded.vcf.recode.vcf.gz DB_chr3S_positions_excluded.vcf.recode.vcf.gz DB_chr4L_positions_excluded.vcf.recode.vcf.gz DB_chr4S_positions_excluded.vcf.recode.vcf.gz DB_chr5L_positions_excluded.vcf.recode.vcf.gz DB_chr5S_positions_excluded.vcf.recode.vcf.gz DB_chr6L_positions_excluded.vcf.recode.vcf.gz DB_chr6S_positions_excluded.vcf.recode.vcf.gz DB_chr7L_positions_excluded.vcf.recode.vcf.gz DB_chr7S_positions_excluded.vcf.recode.vcf.gz DB_chr8L_positions_excluded.vcf.recode.vcf.gz DB_chr8S_positions_excluded.vcf.recode.vcf.gz DB_chr9L_positions_excluded.vcf.recode.vcf.gz DB_chr9S_positions_excluded.vcf.recode.vcf.gz | bgzip -c > l_and_s.vcf.gz

vcf-merge DB_chr1L_positions_excluded.vcf.recode.vcf.gz DB_chr2L_positions_excluded.vcf.recode.vcf.gz DB_chr3L_positions_excluded.vcf.recode.vcf.gz DB_chr4L_positions_excluded.vcf.recode.vcf.gz DB_chr5L_positions_excluded.vcf.recode.vcf.gz DB_chr6L_positions_excluded.vcf.recode.vcf.gz DB_chr7L_positions_excluded.vcf.recode.vcf.gz DB_chr8L_positions_excluded.vcf.recode.vcf.gz DB_chr9L_positions_excluded.vcf.recode.vcf.gz | bgzip -c > l_only.vcf.gz

vcf-merge DB_chr1S_positions_excluded.vcf.recode.vcf.gz DB_chr2S_positions_excluded.vcf.recode.vcf.gz DB_chr3S_positions_excluded.vcf.recode.vcf.gz DB_chr4S_positions_excluded.vcf.recode.vcf.gz DB_chr5S_positions_excluded.vcf.recode.vcf.gz DB_chr6S_positions_excluded.vcf.recode.vcf.gz DB_chr7S_positions_excluded.vcf.recode.vcf.gz DB_chr8S_positions_excluded.vcf.recode.vcf.gz DB_chr9S_positions_excluded.vcf.recode.vcf.gz | bgzip -c > s_only.vcf.gz

```
