# Change this location for outside of the the SING Micronesia
BASE=/projects/teaching/sing_guam

## Make conda available
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O ${BASE}/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -u -p ${BASE}/miniforge3
source ${BASE}/miniforge3/etc/profile.d/conda.sh

## create the conda environment for the workshop
# (mamba is a drop in faster replacement for conda)
mamba create -p ${BASE}/sing_env
conda activate ${BASE}/sing_env
mamba install -c conda-forge -c bioconda samtools bcftools vcftools tabix eigensoft plink bedtools bwa


conda activate ${BASE}/workshop_env
# create directory structures
mkdir -p staging fastq aligned ref vcf

cd ref
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
gunzip human_g1k_v37.fasta.gz
# takes ~40 min on single thread
bwa index human_g1k_v37.fasta

cd ../staging

samtools view -bh https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00102/alignment/HG00102.mapped.ILLUMINA.bwa.GBR.low_coverage.20130415.bam MT > example_mt.bam
samtools sort -n -o example_mt.read_ordered.bam example_mt.bam
# this will generate a lot of warnings about the pair not being next to it in the bam file
bedtools bamtofastq -i example_mt.read_ordered.bam -fq example.R1.fastq -fq2 example.R2.fastq
mv example.R*.fastq ../fastq/

# grab 1000 genomes data
samtools view -bh https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam 4:89051323-89053323 > rs2231142_hom_ref_example.bam
samtools view -bh https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00102/alignment/HG00102.mapped.ILLUMINA.bwa.GBR.low_coverage.20130415.bam 4:89051323-89053323 > rs2231142_het_example.bam
samtools view -bh https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00111/alignment/HG00111.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam 4:89051323-89053323 > rs2231142_hom_alt_example.bam
samtools index -M *.bam
mv rs2231142.* ../aligned


cd ../vcf
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz -O 1kgp.chrMT.vcf.gz
bcftools index -t 1kgp.chrMT.vcf.gz
bcftools view -r 4:89052323-89052324 -s HG00096,HG00102,HG00111 -o example.3_sample.vcf -v snps -O v http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
rm ALL.chr4.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi

wget -O 1kgp-phase3-samples.txt http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

