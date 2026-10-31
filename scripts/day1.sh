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
bash check/day1-envs.sh
git add envs/ scripts/
git commit -m "Day 1: environment files"
git push
git remote -v
echo $GITHUB_TOKEN | head -c 5
git push
git config --global credential.helper store
git push
rm -f ~/.git-credentials
git push
git config --show-origin --get-all credential.helper
env | grep -i token
git -c credential.helper= push
curl -L -o data/airway_samples.csv https://github.com/onedimkurt/bioinfoworkshop/releases/download/data-v1/airway_samples.csv
head -3 data/airway_samples.csv
curl -s https://api.github.com/repos/onedimkurt/bioinfoworkshop/releases | grep -E '"tag_name"|"browser_download_url"'
curl -s https://api.github.com/repos/onedimkurt/bioinfoworkshop | head -5
head -3 data/airway_samples.csv
ls -lh data/airway_samples.csv
bash setup/day1.sh 
ls -lh data/ | head -30
ls data/quant
zcat data/SRR1039508_1.sub.fastq.gz | head -8
for f in data/SRR10395*.sub.fastq.gz; do
echo "$f $(( $(zcat "$f" | wc -l) / 4 )) reads"done
zcat data/SRR1039508_1.sub.fastq.gz | head -4000 | awk 'NR % 4 == 2 { print length($0) }' |sort | uniq -c
mkdir -p results/fastqc
fastqc -t 2 --extract -o results/fastqc data/SRR1039508_1.sub.fastq.gz
fastqc -t 2 --extract -o results/fastqc data/SRR1039508_2.sub.fastq.gz data/SRR1039509_1.sub.fastq.gz data/SRR1039509_2.sub.fastq.gz
cat results/fastqc/*_fastqc/summary.txt | sort -k2,2 -t$'\t' | cut -f1,2,3
mkdir -p results/fastp
for s in SRR1039508 SRR1039509; do
fastp -i data/${s}_1.sub.fastq.gz -I data/${s}_2.sub.fastq.gz \
-o results/fastp/${s}_1.trim.fastq.gz -O results/fastp/${s}_2.trim.fastq.gz \
--detect_adapter_for_pe -w 2 \
-j results/fastp/${s}.fastp.json -h results/fastp/${s}.fastp.html
done
mkdir -p results/fastqc_trimmed
fastqc -t 2 --extract -o results/fastqc_trimmed results/fastp/*.trim.fastq.gz
cat results/fastqc_trimmed/*_fastqc/summary.txt | sort -k2,2 -t$'\t' | cut -f1,2,3
salmon index -t data/gencode.v50.transcripts.chr.fa.gz -i results/salmon_index --gencode -p 2
for s in SRR1039508 SRR1039509; do
salmon quant -i results/salmon_index -l A \
-1 results/fastp/${s}_1.trim.fastq.gz -2 results/fastp/${s}_2.trim.fastq.gz \
-p 2 -o results/salmon/${s}
done
grep -o '"percent_mapped": [0-9.]*' results/salmon/*/aux_info/meta_info.json
grep -h '"expected_format"\|"strand_mapping_bias"' results/salmon/*/lib_format_counts.json
head -5 results/salmon/SRR1039508/quant.sf
for d in data/quant/SRR*; do
echo "$(basename $d) $(grep -o '"percent_mapped": [0-9.]*' $d/aux_info/meta_info.json)"
done
quant_508 <- read.delim("results/salmon/SRR1039508/quant.sf")
long_tx <- quant_508$Length >= 3000
length_gap <- quant_508$Length[long_tx] - quant_508$EffectiveLength[long_tx]
summary(length_gap)
multiqc -f -q -o results/multiqc results/fastqc results/fastqc_trimmed results/fastp results/salmon
awk -F'\t' 'NR == 1 || $2 != ""' results/multiqc/multiqc_data/multiqc_general_stats.txt | cut -f1-4
git add scripts/
git commit -m "Day 1: reads to counts"
git push
git remote set-url origin https://sudesubasi8-dev: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@github.com/sudesubasi8-dev/Bioinfo-workshop.git
git push
git add scripts/
git commit -m "Day 1: differential expression"
git push
git log --oneline
git status
mkdir -p results/robustness
salmon quant -i results/salmon_index -l A \
-1 data/SRR1039508_1.sub.fastq.gz -2 data/SRR1039508_2.sub.fastq.gz \
-p 2 -o results/robustness/SRR1039508_raw > /dev/null 2>&1
echo "untrimmed: $(grep -o '"percent_mapped": [0-9.]*' results/robustness/SRR1039508_raw/aux_info/meta_info.json)"
echo "trimmed: $(grep -o '"percent_mapped": [0-9.]*' results/salmon/SRR1039508/aux_info/meta_info.json)"
raw_vs_trimmed <- tximport(c(trimmed = "results/salmon/SRR1039508/quant.sf", raw = "results/robustness/SRR1039508_raw/quant.sf"), type = "salmon", tx2gene = tx2gene)$counts
