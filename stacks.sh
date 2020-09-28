

###Copy raw reads files from CGRB server to personal directory for analysis
SGE_Batch -c 'cp /nfs2/hts/illumina/191003_J00107_0216_AHF3LGBBXY_1476/L12356/lane5-s044-index--GBS0239_S44_L005_R1_001.fastq.gz /raid1/home/bpp/gannondu/HJA_Columbine/raw_reads/raw_round2' -r reads_plate1 -m 16G -P 1

SGE_Batch -c 'cp /nfs2/hts/illumina/191003_J00107_0216_AHF3LGBBXY_1476/L12356/lane6-s045-index--GBS0240_S45_L006_R1_001.fastq.gz /raid1/home/bpp/gannondu/HJA_Columbine/raw_reads/raw_round2' -r reads_plate2 -m 16G -P 1

##Process Radtags through Stacks using default filtering parameters

SGE_Batch -c 'process_radtags -f /raid1/home/bpp/gannondu/HJA_Columbine/raw_reads/raw_round2/lane5-s044-index--GBS0239_S44_L005_R1_001.fastq.gz -o /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/ -b /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/barcodes1_round2.tsv -i gzfastq -c -q -r -e apeKI' -r process_radtags1 
SGE_Batch -c 'process_radtags -f /raid1/home/bpp/gannondu/HJA_Columbine/raw_reads/raw_round2/lane6-s045-index--GBS0240_S45_L006_R1_001.fastq.gz -o /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/ -b /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/barcodes2_round2.tsv -i gzfastq -c -q -r -e apeKI' -r process_radtags2 







######################################################################
 ## Trial runs to align sequences to A. coerulea reference genome 
######################################################################

#Bowtie2
  
  # Arbitrarily selected 10 samples to check alignment rates with default parameters
  
   SGE_Batch -c "bowtie2 -x A_coerulea -U /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s193.fq.gz -S /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s193.sam" -r bowtie2_s193_def -P 1
   SGE_Batch -c "bowtie2 -x A_coerulea -U /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s201.fq.gz -S /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s201.sam" -r bowtie2_s201_def -P 1
   SGE_Batch -c "bowtie2 -x A_coerulea -U /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s210.fq.gz -S /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s210.sam" -r bowtie2_s210_def -P 1
   SGE_Batch -c "bowtie2 -x A_coerulea -U /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s428.fq.gz -S /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s428.sam" -r bowtie2_s428_def -P 1
   SGE_Batch -c "bowtie2 -x A_coerulea -U /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s231.fq.gz -S /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s231.sam" -r bowtie2_s231_def -P 1
   SGE_Batch -c "bowtie2 -x A_coerulea -U /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s235.fq.gz -S /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s235.sam" -r bowtie2_s235_def -P 1
   SGE_Batch -c "bowtie2 -x A_coerulea -U /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s250.fq.gz -S /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s250.sam" -r bowtie2_s250_def -P 1
   SGE_Batch -c "bowtie2 -x A_coerulea -U /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s267.fq.gz -S /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s267.sam" -r bowtie2_s267_def -P 1
   SGE_Batch -c "bowtie2 -x A_coerulea -U /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s281.fq.gz -S /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s281.sam" -r bowtie2_s281_def -P 1
   SGE_Batch -c "bowtie2 -x A_coerulea -U /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s472.fq.gz -S /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s472.sam" -r bowtie2_s472_def -P 1

   SGE_Batch -c "samtools flagstat /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s193.sam > /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/bowtie2/default_parms/s193_def_stats.txt" -r s193_stats -P 1
   SGE_Batch -c "samtools flagstat /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s201.sam > /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/bowtie2/default_parms/s201_def_stats.txt" -r s201_stats -P 1
   SGE_Batch -c "samtools flagstat /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s210.sam > /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/bowtie2/default_parms/s210_def_stats.txt" -r s210_stats -P 1
   SGE_Batch -c "samtools flagstat /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s428.sam > /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/bowtie2/default_parms/s428_def_stats.txt" -r s428_stats -P 1
   SGE_Batch -c "samtools flagstat /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s231.sam > /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/bowtie2/default_parms/s231_def_stats.txt" -r s231_stats -P 1
   SGE_Batch -c "samtools flagstat /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s235.sam > /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/bowtie2/default_parms/s235_def_stats.txt" -r s235_stats -P 1
   SGE_Batch -c "samtools flagstat /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s250.sam > /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/bowtie2/default_parms/s250_def_stats.txt" -r s250_stats -P 1
   SGE_Batch -c "samtools flagstat /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s267.sam > /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/bowtie2/default_parms/s267_def_stats.txt" -r s267_stats -P 1
   SGE_Batch -c "samtools flagstat /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s281.sam > /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/bowtie2/default_parms/s281_def_stats.txt" -r s281_stats -P 1
   SGE_Batch -c "samtools flagstat /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s472.sam > /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/bowtie2/default_parms/s472_def_stats.txt" -r s472_stats -P 1

# Try local alignment too

   SGE_Batch -c "bowtie2 --local -x A_coerulea -U /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s193.fq.gz -S /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s193.sam" -r bowtie2_s193_loc -P 1
   SGE_Batch -c "bowtie2 --local -x A_coerulea -U /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s201.fq.gz -S /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s201.sam" -r bowtie2_s201_loc -P 1
   SGE_Batch -c "bowtie2 --local -x A_coerulea -U /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s210.fq.gz -S /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s210.sam" -r bowtie2_s210_loc -P 1
   SGE_Batch -c "bowtie2 --local -x A_coerulea -U /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s428.fq.gz -S /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s428.sam" -r bowtie2_s428_loc -P 1
   SGE_Batch -c "bowtie2 --local -x A_coerulea -U /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s231.fq.gz -S /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s231.sam" -r bowtie2_s231_loc -P 1
   SGE_Batch -c "bowtie2 --local -x A_coerulea -U /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s235.fq.gz -S /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s235.sam" -r bowtie2_s235_loc -P 1
   SGE_Batch -c "bowtie2 --local -x A_coerulea -U /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s250.fq.gz -S /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s250.sam" -r bowtie2_s250_loc -P 1
   SGE_Batch -c "bowtie2 --local -x A_coerulea -U /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s267.fq.gz -S /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s267.sam" -r bowtie2_s267_loc -P 1
   SGE_Batch -c "bowtie2 --local -x A_coerulea -U /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s281.fq.gz -S /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s281.sam" -r bowtie2_s281_loc -P 1
   SGE_Batch -c "bowtie2 --local -x A_coerulea -U /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s472.fq.gz -S /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s472.sam" -r bowtie2_s472_loc -P 1 

# sort alignments
   
   SGE_Batch -c "samtools sort -m 16G -@ 2 -o /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s193.sorted.bam /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s193.sam" -r sort_s193_loc -P 2
   SGE_Batch -c "samtools sort -m 16G -@ 2 -o /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s201.sorted.bam /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s201.sam" -r sort_s201_loc -P 2
   SGE_Batch -c "samtools sort -m 16G -@ 2 -o /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s210.sorted.bam /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s210.sam" -r sort_s210_loc -P 2
   SGE_Batch -c "samtools sort -m 16G -@ 2 -o /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s428.sorted.bam /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s428.sam" -r sort_s428_loc -P 2
   SGE_Batch -c "samtools sort -m 16G -@ 2 -o /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s231.sorted.bam /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s231.sam" -r sort_s231_loc -P 2
   SGE_Batch -c "samtools sort -m 16G -@ 2 -o /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s235.sorted.bam /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s235.sam" -r sort_s235_loc -P 2
   SGE_Batch -c "samtools sort -m 16G -@ 2 -o /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s250.sorted.bam /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s250.sam" -r sort_s250_loc -P 2
   SGE_Batch -c "samtools sort -m 16G -@ 2 -o /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s267.sorted.bam /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s267.sam" -r sort_s267_loc -P 2
   SGE_Batch -c "samtools sort -m 16G -@ 2 -o /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s281.sorted.bam /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s281.sam" -r sort_s281_loc -P 2
   SGE_Batch -c "samtools sort -m 16G -@ 2 -o /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s472.sorted.bam /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/s472.sam" -r sort_s472_loc -P 2
  
# stacks runs

   SGE_Batch -c "gstacks -I /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/ -M /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/popmap_tests.tsv -O /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/gstacks2/ -t 16 --min-mapq 20" -r gstacks -P 16
   
   SGE_Batch -c "populations -P /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/ -O /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/poptest1/ -M /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/popmap_tests.tsv -t 16 -R 1 --write-random-snp --vcf" -r populations1 -P 16
  
  
# Try bwa

   SGE_Batch -c "bwa index -p A_coerulea /raid1/home/bpp/gannondu/HJA_Columbine/A_coerulea/Acoerulea_195_v1.fa" -r bwa_index -P 1
  
   SGE_Batch -c "bwa mem A_coerulea /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s193.fq.gz > /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s193.sam" -r bwa_s193 -P 1
   SGE_Batch -c "bwa mem A_coerulea /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s201.fq.gz > /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s201.sam" -r bwa_s201 -P 1
   SGE_Batch -c "bwa mem A_coerulea /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s210.fq.gz > /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s210.sam" -r bwa_s210 -P 1
   SGE_Batch -c "bwa mem A_coerulea /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s428.fq.gz > /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s428.sam" -r bwa_s428 -P 1
   SGE_Batch -c "bwa mem A_coerulea /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s231.fq.gz > /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s231.sam" -r bwa_s231 -P 1
   SGE_Batch -c "bwa mem A_coerulea /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s235.fq.gz > /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s235.sam" -r bwa_s235 -P 1
   SGE_Batch -c "bwa mem A_coerulea /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s250.fq.gz > /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s250.sam" -r bwa_s250 -P 1
   SGE_Batch -c "bwa mem A_coerulea /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s267.fq.gz > /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s267.sam" -r bwa_s267 -P 1
   SGE_Batch -c "bwa mem A_coerulea /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s281.fq.gz > /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s281.sam" -r bwa_s281 -P 1
   SGE_Batch -c "bwa mem A_coerulea /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/s472.fq.gz > /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s472.sam" -r bwa_s472 -P 1 
  
# sort the alignments

   SGE_Batch -c "samtools sort -m 16G -@ 2 -o /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s193.bam /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s193.sam" -r sort_s193 -P 2
   SGE_Batch -c "samtools sort -m 16G -@ 2 -o /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s201.bam /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s201.sam" -r sort_s201 -P 2
   SGE_Batch -c "samtools sort -m 16G -@ 2 -o /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s210.bam /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s210.sam" -r sort_s210 -P 2
   SGE_Batch -c "samtools sort -m 16G -@ 2 -o /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s428.bam /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s428.sam" -r sort_s428 -P 2
   SGE_Batch -c "samtools sort -m 16G -@ 2 -o /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s231.bam /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s231.sam" -r sort_s231 -P 2
   SGE_Batch -c "samtools sort -m 16G -@ 2 -o /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s235.bam /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s235.sam" -r sort_s235 -P 2
   SGE_Batch -c "samtools sort -m 16G -@ 2 -o /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s250.bam /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s250.sam" -r sort_s250 -P 2
   SGE_Batch -c "samtools sort -m 16G -@ 2 -o /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s267.bam /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s267.sam" -r sort_s267 -P 2
   SGE_Batch -c "samtools sort -m 16G -@ 2 -o /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s281.bam /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s281.sam" -r sort_s281 -P 2
   SGE_Batch -c "samtools sort -m 16G -@ 2 -o /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s472.bam /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/s472.sam" -r sort_s472 -P 2
  
 # stacks
 
   SGE_Batch -c "gstacks -I /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/ -M /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/popmap_tests.tsv -O /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/ -t 16 --min-mapq 20" -r gstacks -P 16
   
   SGE_Batch -c "populations -P /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/ -O /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bwa/ -M /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/Tests/bowtie2/local/popmap_tests.tsv -t 16 -R 1 --write-random-snp --vcf" -r populations1 -P 16

######################################################################
 ## BWA seems to work a bit better (more SNPs out on the back end)  
######################################################################






######################################################################
 ## Continue with full set of samples 
######################################################################

#Run shell script to loop through sequence files and align them to reference genome

  SGE_Batch -c "/raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/bwa/align_all.sh" -r align_all -P 8 -m 16G
  

#Run shell script to loop through sam files and sort them left-to-right on the reference genome

  SGE_Batch -c "/raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/bwa/sort_all.sh" -r sort_all -P 2
  
# Run shell script to compute mean coverage across the genome for each sample

  SGE_Batch -c "/raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/bwa/coverage.sh" -r coverage
  
# Run shell script to compute mean depth across sequenced loci

  SGE_Batch -c "/raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/bwa/depth_all.sh" -r depth
  
# Run stacks, filtering for high-quality read alignments

  SGE_Batch -c "gstacks -I /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/bwa/ -M /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/popmap_round2.tsv -O /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/ -t 16 --min-mapq 20" -r gstacks -P 16

  SGE_Batch -c "populations -P /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/ -O /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/populations_r80merge/ -M /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/popmap_round2.tsv -t 16 -r 0.85 -p 4 --write-random-snp --merge-sites -e apeKI --vcf --ordered-export" -r populations1 -P 16

  SGE_Batch -c "populations -P /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/ -O /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/populations_nomiss -M /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/popmap_round2.tsv -t 16 -R 1 --write-random-snp --merge-sites -e apeKI --vcf" -r populations2 -P 16


##########################################################################################

## VCFTools for summary stats

SGE_Batch -c 'vcftools --vcf /raid1/home/bpp/gannondu/HJA_Columbine/stacks_populations/batch_1.vcf --out populations1 --depth' -r vcftools_depth -P 8 -F 8G -m 16G

SGE_Batch -c 'vcftools --vcf /raid1/home/bpp/gannondu/HJA_Columbine/stacks_populations/batch_1.vcf --out populations1 --hardy' -r vcftools_hardy -P 8 -F 8G -m 16G

SGE_Batch -c 'vcftools --vcf /raid1/home/bpp/gannondu/HJA_Columbine/stacks_populations/batch_1.vcf --out populations1 --het' -r vcftools_het -P 8 -F 8G -m 16G

##########################################################################################

##structure

SGE_Batch -c 'structure' -r structure_k1 -P 16 -m 32G

SGE_Batch -c 'structure -K 2 -o struct_k2' -r structure_k2 -P 16 -m 32G

SGE_Batch -c 'structure -K 3 -o struct_k3' -r structure_k3 -P 16 -m 32G

SGE_Batch -c 'structure -K 4 -o struct_k4' -r structure_k4 -P 16 -m 32G

SGE_Batch -c 'structure -K 5 -o struct_k5' -r structure_k5 -P 16 -m 32G





##### Runs for deltaK method

for i in {1..20}; do echo "SGE_Batch -c 'structure -m mainparams2 -o 58struct_k1.$i' -r structure58_k1.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 2 -o 58struct_k2.$i -m mainparams2' -r structure58_k2.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 3 -o 58struct_k3.$i -m mainparams2' -r structure58_k3.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 4 -o 58struct_k4.$i -m mainparams2' -r structure58_k4.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 5 -o 58struct_k5.$i -m mainparams2' -r structure58_k5.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 6 -o 58struct_k6.$i -m mainparams2' -r structure58_k6.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 7 -o 58struct_k7.$i -m mainparams2' -r structure58_k7.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 8 -o 58struct_k8.$i -m mainparams2' -r structure58_k8.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 9 -o 58struct_k9.$i -m mainparams2' -r structure58_k9.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 10 -o 58struct_k10.$i -m mainparams2' -r structure58_k10.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 11 -o 58struct_k11.$i -m mainparams2' -r structure58_k11.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 12 -o 58struct_k12.$i -m mainparams2' -r structure58_k12.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 13 -o 58struct_k13.$i -m mainparams2' -r structure58_k13.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 14 -o 58struct_k14.$i -m mainparams2' -r structure58_k14.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 15 -o 58struct_k15.$i -m mainparams2' -r structure58_k15.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 16 -o 58struct_k16.$i -m mainparams2' -r structure58_k16.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 17 -o 58struct_k17.$i -m mainparams2' -r structure58_k17.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 18 -o 58struct_k18.$i -m mainparams2' -r structure58_k18.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 19 -o 58struct_k19.$i -m mainparams2' -r structure58_k19.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 20 -o 58struct_k20.$i -m mainparams2' -r structure58_k20.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 21 -o 58struct_k21.$i -m mainparams2' -r structure58_k21.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 22 -o 58struct_k22.$i -m mainparams2' -r structure58_k22.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 23 -o 58struct_k23.$i -m mainparams2' -r structure58_k23.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 24 -o 58struct_k24.$i -m mainparams2' -r structure58_k24.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 25 -o 58struct_k25.$i -m mainparams2' -r structure58_k25.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 26 -o 58struct_k26.$i -m mainparams2' -r structure58_k26.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 27 -o 58struct_k27.$i -m mainparams2' -r structure58_k27.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 28 -o 58struct_k28.$i -m mainparams2' -r structure58_k28.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 29 -o 58struct_k29.$i -m mainparams2' -r structure58_k29.$i -P 2"; done

for i in {1..20}; do echo "SGE_Batch -c 'structure -K 30 -o 58struct_k30.$i -m mainparams2' -r structure58_k30.$i -P 2"; done





##########################################################################################


##Bayescan for identification of outlier loci

#Run 1 with population IDs = to meadows
SGE_Batch -c 'bayescan bayescan_680snps -o bayescan1 -out_pilot -out_freq' -r bayescan1 -P 8 -m 16G


#Run 2 with population IDs = to structure clusters
SGE_Batch -c 'bayescan bayescan2_structpop -o bayescan2 -out_pilot -out_freq' -r bayescan2 -P 8 -m 16G


#Run 3 with population IDs = to meadow complexes

SGE_Batch -c 'bayescan bayescan3_complexes -o bayescan3 -out_pilot -out_freq' -r bayescan3 -P 8 -m 16G


#Run 3 with population IDs = to xeric/mesic/fern meadows

SGE_Batch -c 'bayescan bayescan4_xericmesic -o bayescan4 -out_pilot -out_freq' -r bayescan4 -P 8 -m 16G























