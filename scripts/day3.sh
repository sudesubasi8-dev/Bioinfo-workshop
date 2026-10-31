conda env create -f envs/variants.yml
conda activate variants
bwa 2>&1 | grep Version
samtools --version | head -1
bcftools --version | head -1
which bgzip tabix
