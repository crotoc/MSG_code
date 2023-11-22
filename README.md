# MSG

This code is the implemention of MSG methods published as following:

Y. Ji, Q. Wei, R. Chen, Q. Wang, R. Tao, B. Li, Integration of multidimensional splicing data and GWAS summary statistics for risk gene discovery. PLoS Genet. 18, e1009814 (2022).


# Hassle-free way to run the code: singularity

singularity > 3.8

```{bash}
conda create -p /path_to_env/msg

conda activate /path_to_env/msg

conda install -c conda-forge singularity

```


## step to run it

1. Download the data, MSG.sif from the [link](https://www.dropbox.com/sh/cysiek3bbo8eo4r/AAAjoAhwHRkl0wnGIjie9Kv8a?dl=0):

	- MSG.sif is a basic ubuntu environments build including all dependencies needed by MSG

	- gene.500k.id: The file include the coordinates of upstreaming and downstreaming 500k of a gene

	- hg38.vcf: The rsid for each snp in 1k genome.

	- GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz: The genotype data for all the individuals in GTEx. Used for generating x matrix.

	- samples.filtered.used: The samples used for running the MSG methods. Need to be generated for each tissue by users. Users can use all samples in the genotype file. It will be automatically filtered during the excution.

	- LDREF: The linkage disequilibrim file for each chromosome.

	- rawSplice: The splicing expression in GTEx. Used for generating y matrix.

	- sumstats: The summary statistics of a specific disease. Users need to format their own file by this way. 
	

3. Git the code to local direcoty:

	- scripts: include all the scripts for the MSG method.

	- msg.pl: The wrapper for all steps including generate x matrix, y matrix and cov file and run MSG.

	- generate_db_and_cov.R: The R script to generate cov file for each gene.

	- gtex_comp_MSG_ACAT_GBJ_120522.R: The R script run the MSG method.


4. Put the data and the code in the same dir. If not, please specify absolute path. Run for gene ENSG00000000457:


```{bash}
## clone scripts
git clone https://github.com/crotoc/MSG_code.git

cd MSG_code

## !!!! Download the MSG.sif and data dir to MSG_code dir

## activate the singularity env
conda activate /fs0/chenr6/chenr6/opt/condaEnv/msg

## run one gene to test
singularity run -e --bind /nobackup,/home,/fs0  MSG.sif perl scripts/msg.pl -cmd generateXmatrix -cmd generateYmatrix -cmd generateCov -cmd runMSG -cmd cleanIntermediate -gene ENSG00000000457 -file_genePos data/gene.500k.id -file_snp data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz -file_rsid data/hg38.vcf  -file_sample data/samples.filtered.used -file_splicing data/rawSplice/Whole_Blood.v8.leafcutter_phenotypes.bed.gz --file_ldRef data/LDREF/1000G.EUR. -sumstats data/sumstats/clozukscz.sumstats -dir_script scripts/ -e

## check the results in out
cat out/results/all/ENSG00000000457.select.rsid.vcf.ac.results.MSG_GBJ_ACAT.txt

## generate cmds for all genes
singularity run -e --bind /nobackup,/home,/fs0  MSG.sif parallel -j 1 -q echo "singularity run -e --bind /nobackup,/home,/fs0  MSG.sif perl scripts/msg.pl -cmd generateXmatrix -cmd generateYmatrix -cmd generateCov -cmd runMSG -cmd cleanIntermediate -gene {} -file_genePos data/gene.500k.id -file_snp data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz -file_rsid data/hg38.vcf  -file_sample data/samples.filtered.used -file_splicing data/rawSplice/Whole_Blood.v8.leafcutter_phenotypes.bed.gz --file_ldRef data/LDREF/1000G.EUR. -sumstats data/sumstats/clozukscz.sumstats -dir_script scripts/ -e" ::: $(less data/gene.500k.id | grep protein_coding | cut -f 1 | sort -k1,1V ) > all.cmd

```

5. Get helps of msg.pl by this way:

```
cd MSG_code

singularity run -e --bind /nobackup,/home,/fs0 MSG.sif perldoc scripts/msg.pl 
```


### Be aware

#### Make sure -cmd has the correct order: 

- generateXmatrix
- generateYmatrix
- generateCov
- runMSG

generateXmatrix should be generated before generateCov

generateXmatrix, generateYmatrix and generateCov should be generated before runMSG.


## The definition file of how to consturct MSG.sif

scripts/msg.singularity

Be aware of this needs root priviledge, so you have to do it on a local machine.

```{bash, label = "", linewidth = 85, eval=opt$eval}

sudo /home/chenr6/miniforge3/bin/singularity build MSG.sif scripts/msg.singularity

```








## If you need to install the package by your self

Please referring scripts/msg.singularity


### Requirements

1. perl modules

	- Switch
	- IPC::Run
	- DBI
	- List::Uniq
	- Array::Utils
	- String/Random
	- MCE/Hobo.pm
	- Env/Modify.pm

2. Two customized perl modules
	
- Bundle-Wrapper: https://github.com/crotoc/Bundle-Wrapper

3. tabix: https://github.com/samtools/htslib

4. parallel: https://www.gnu.org/software/parallel/

5. R packages:

	- data.table
	- optparse
	- magrittr
	- reshape2
	- pma
	- mvtnorm
	- rcpp
	- devtools
	- rcppeigen
	- GBJ
	- glmnet
	- this.path
	- dplyr
	- plink2R
	- TiscoMM


### step-by-step instruction

```{bash, label = "", linewidth = 85, eval=opt$eval}
## 1. install conda
## install mambaforge - a variant of conda
ENV=/Path_to_Conda/miniconda3

mkdir download/
cd download/
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
mkdir -p $ENV
bash Mambaforge-Linux-x86_64.sh -b -p $ENV

export PATH=$PATH:$ENV/bin

## 2. install R, perl, cpanm, parallel and other r packages.
### activate the env
source $ENV/bin/activate
### install
mamba install -y -c conda-forge r-base perl perl-app-cpanminus parallel r-data.table r-optparse r-magrittr r-reshape2 r-pma r-mvtnorm  r-rcpp r-devtools r-stringi r-rcppeigen r-gbj r-glmnet r-this.path r-dplyr r-devtools

## 3. install perl modules needed by msg.pl
mamba install -y -c bioconda perl-dbi perl-ipc-run perl-mce perl-module-build perl-string-random perl-mce-shared
cpanm Switch
cpanm List::Uniq
cpanm Array/Utils.pm
cpanm MCE/Hobo.pm
cpanm --force Env/Modify.pm
wget https://github.com/crotoc/Bundle-Wrapper/raw/master/Bundle-Wrapper-0.03.tar.gz
cpanm Bundle-Wrapper-0.03.tar.gz

## 4. install tabix
mamba install -y -c bioconda tabix

## 5. install r packages need by MSG
R -e 'library(devtools);install_github("gabraham/plink2R", subdir="plink2R")'

mamba install -y r-bigmemory r-bigmemory.sri  r-bh r-uuid
R -e 'library(devtools);install_github("XingjieShi/TisCoMM")'

## 6. link compilers
parallel echo ln -s $ENV/bin/x86_64-conda_cos6-linux-gnu-cc $ENV/bin/{=s/.*-//=} ::: $ENV/bin/x86_64-conda_cos6-linux-gnu-cc $ENV/bin/x86_64-conda_cos6-linux-gnu-gcc $ENV/bin/x86_64-conda_cos6-linux-gnu-g++ | bash

## activate the ENV
conda activate $ENV


## run MSG
perl scripts/msg.pl -cmd generateXmatrix -cmd generateYmatrix -cmd generateCov -cmd runMSG -cmd cleanIntermediate -gene ENSG00000000457 -file_genePos data/gene.500k.id -file_snp data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz -file_rsid data/hg38.vcf  -file_sample data/samples.filtered.used -file_splicing data/rawSplice/Whole_Blood.v8.leafcutter_phenotypes.bed.gz --file_ldRef data/LDREF/1000G.EUR. -sumstats data/sumstats/clozukscz.sumstats -dir_script scripts/ -e
```





