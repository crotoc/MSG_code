# MSG

This code is the implemention of MSG methods published as following:

Y. Ji, Q. Wei, R. Chen, Q. Wang, R. Tao, B. Li, Integration of multidimensional splicing data and GWAS summary statistics for risk gene discovery. PLoS Genet. 18, e1009814 (2022).


# Hassle-free way to run the code: singularity

singularity > 3.8

```{bash}
conda create -p /fs0/chenr6/chenr6/opt/condaEnv/msg

conda install -c conda-forge singularity

```


## step to run it

1. Download the MSG.sif and overlay.img from the link:

	- MSG.sif is a basic ubuntu environments build based on the config file: BASE.singularity

	- overlay.img is a image include all the required library for running the script

2. Download the data files in the methods:

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
singularity run -e --bind /nobackup,/home,/fs0 --overlay overlay.img MSG.sif perl /fs0/chenr6/chenr6/myscript/script/msg.pl -cmd generateXmatrix -cmd generateYmatrix -cmd generateCov -cmd runMSG  -gene ENSG00000000457 file_genePos data/gene.500k.id -file_snp data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz  -file_sample data/samples.filtered.used -file_splicing data/rawSplice/Whole_Blood.v8.leafcutter_phenotypes.bed.gz --file_ldRef data/LDREF/1000G.EUR. -sumstats data/sumstats/clozukscz.sumstats -dir_script scripts/ -e
```

5. Get helps of msg.pl by this way:

```
singularity run -e --bind /nobackup,/home,/fs0 --overlay overlay.img MSG.sif perldoc /fs0/chenr6/chenr6/myscript/script/msg.pl 
```


### Be aware

#### Make sure -cmd has the correct order: 

- generateXmatrix
- generateYmatrix
- generateCov
- runMSG

generateXmatrix should be generated before generateCov

generateXmatrix, generateYmatrix and generateCov should be generated before runMSG.



## If you need to install the package by your self


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
	
	- MyBase-Mysub: https://github.com/crotoc/MyBase-Mysub
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
	


