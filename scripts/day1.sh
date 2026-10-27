conda config --show channels
conda create -y -n deneme seqkit

conda activate deneme
which seqkit

seqkit version

conda deactivate
conda list -n variants-ready | grep -E 

conda install -y -n deneme seqtk
conda activate deneme
which
conda env export -n deneme --from-history > envs/deneme.yml
cat envs/deneme.yml
conda env remove -y -n deneme

conda env list
conda deactivate 
conda env remove -y -n deneme 
conda env list
conda list -n variants-ready | grep -E "^(# Name|bcftools|samtools|bwa|htslib) "
cat envs/rnaseq.yml
conda env list
cat envs/rnaseq.yml
conda env create -f envs/rnaseq.yml
conda activate rnaseq
salmon --version
fastp --version
fastqc --version
multiqc --version


