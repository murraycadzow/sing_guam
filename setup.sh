# $1 is the account list file

BASE=/projects/teaching/sing_guam

## Make conda available
mkdir -p ${BASE}/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ${BASE}/miniconda3/miniconda.sh
bash ${BASE}/miniconda3/miniconda.sh -b -u -p ${BASE}/miniconda3
rm -rf ${BASE}/miniconda3/miniconda.sh

## Create conda environment
source ${BASE}/miniconda3/etc/profile.d/conda.sh
conda create -p ${BASE}/workshop_env
conda activate ${BASE}/workshop_env
conda install --solver=libmamba -c conda-forge -c bioconda samtools bcftools vcftools tabix eigensoft plink bedtools bwa

# create directory structures
mkdir -p workshop/{fastq,aligned,}


# grab 1000 genomes data
samtools -bh 

## Create a working area for each account
for account in $(cat $1)
do
	echo ${account}
done
