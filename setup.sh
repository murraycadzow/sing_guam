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

conda activate ${BASE}/workshop_env
# create directory structures
mkdir -p staging fastq aligned ref vcf

cd ref
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
gzip human_g1k_v37.fasta.gz
bwa index human_g1k_v37.fasta

cd ../staging

samtools -bh https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00102/alignment/HG00102.mapped.ILLUMINA.bwa.GBR.low_coverage.20130415.bam MT > example.bam
samtools sort -n -o example.read_ordered.bam
bedtools bamtofastq -i example_mt.read_ordered.bam -fq example.R1.fastq -fq2 example.R2.fastq
mv example.R*.fastq ../fastq/

cd ../aligned/
# grab 1000 genomes data
samtools -bh https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam 4:89051323-89053323 > rs2231142_hom_ref_example.bam
samtools -bh https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00102/alignment/HG00102.mapped.ILLUMINA.bwa.GBR.low_coverage.20130415.bam 4:89051323-89053323 > rs2231142_het_example.bam
samtools -bh https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00111/alignment/HG00111.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam 4:89051323-89053323 > rs2231142_hom_alt_example.bam
samtools index -m *.bam


cd ../vcf
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz
mv ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz 1kgp.chrMT.vcf.gz
bcftools index -t 1kgp.chrMT.vcf.gz
bcftools view -r 4:89052323-89052324 -s HG00096,HG00102,HG00111 -o example.3_sample.vcf -v snps -O v http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
rm ALL.chr4.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi

wget -o 1kgp-phase3-samples.txt http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

## Create a working area for each account
for account in $(cat $1)
do
	echo ${account}
done
