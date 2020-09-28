

###Copy raw reads files from CGRB server to personal directory for analysis
  SGE_Batch -c 'cp /nfs2/hts/illumina/191003_J00107_0216_AHF3LGBBXY_1476/L12356/lane5-s044-index--GBS0239_S44_L005_R1_001.fastq.gz /raid1/home/bpp/gannondu/HJA_Columbine/raw_reads/raw_round2' -r reads_plate1 -m 16G -P 1

  SGE_Batch -c 'cp /nfs2/hts/illumina/191003_J00107_0216_AHF3LGBBXY_1476/L12356/lane6-s045-index--GBS0240_S45_L006_R1_001.fastq.gz /raid1/home/bpp/gannondu/HJA_Columbine/raw_reads/raw_round2' -r reads_plate2 -m 16G -P 1



###Process Radtags through Stacks using default filtering parameters
  SGE_Batch -c 'process_radtags -f /raid1/home/bpp/gannondu/HJA_Columbine/raw_reads/raw_round2/lane5-s044-index--GBS0239_S44_L005_R1_001.fastq.gz -o /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/ -b /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/barcodes1_round2.tsv -i gzfastq -c -q -r -e apeKI' -r process_radtags1 
  SGE_Batch -c 'process_radtags -f /raid1/home/bpp/gannondu/HJA_Columbine/raw_reads/raw_round2/lane6-s045-index--GBS0240_S45_L006_R1_001.fastq.gz -o /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/ -b /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/barcodes2_round2.tsv -i gzfastq -c -q -r -e apeKI' -r process_radtags2 




######################################################################
 ## Aligning, sorting, and basic coverage stats 
######################################################################

#Run shell script to loop through sequence files and align them to reference genome
  SGE_Batch -c "/raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/bwa/align_all.sh" -r align_all -P 8 -m 16G
  
#Run shell script to loop through sam files and sort them left-to-right on the reference genome
  SGE_Batch -c "/raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/bwa/sort_all.sh" -r sort_all -P 2
  
# Run shell script to compute mean coverage across the genome for each sample
  SGE_Batch -c "/raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/bwa/coverage.sh" -r coverage
  
# Run shell script to compute mean depth across sequenced loci
  SGE_Batch -c "/raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/bwa/depth_all.sh" -r depth
  
  
  
  
######################################################################
 ## Stacks commands for creating locus catalog and vcf
######################################################################
  
# Run gstacks, filtering for high-quality read alignments
  SGE_Batch -c "gstacks -I /raid1/home/bpp/gannondu/HJA_Columbine/alignments_round2/bwa/ -M /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/popmap_round2.tsv -O /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/ -t 16 --min-mapq 20" -r gstacks -P 16

# Stacks populations
  #merge loci at adjacent cut sites
  #1 random snp per locus
  #snps must be present in 85% of samples in a population, and all four populations
  SGE_Batch -c "populations -P /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/ -O /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/populations_r80merge/ -M /raid1/home/bpp/gannondu/HJA_Columbine/stacks_round2/popmap_round2.tsv -t 16 -r 0.85 -p 4 --write-random-snp --merge-sites -e apeKI --vcf --ordered-export" -r populations1 -P 16


##########################################################################################



