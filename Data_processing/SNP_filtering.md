SNP filtering using vcfR
================
D. G. Gannon
September 2020

Previous quality filters up to this point include:

1)  Demultiplexing phase: Removing reads for which Phred score falls
    below 10 in a 15% sliding window ( program)

2)  Demultiplexing phase: Removing reads with an uncalled base (
    program)

3)  Catalog phase: Discard alignments with Phred scores \< 20 (i.e. only
    keep alignments for which ![\\text{Pr}(\\text{Incorrect alignment}
    \< 0.01)](https://latex.codecogs.com/png.latex?%5Ctext%7BPr%7D%28%5Ctext%7BIncorrect%20alignment%7D%20%3C%200.01%29
    "\\text{Pr}(\\text{Incorrect alignment} \< 0.01)");  program)

4)  Genotyping: Loci must be present in at least 80% of individuals
    sampled and have a minor allele frequency of at least 0.05

5)  Locus filters: vcf file is filtered to include only bi-allelic loci
    using .

<!-- end list -->

``` bash
  vcftools --vcf ./AQFO_snps_R80.vcf --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all \
  --out ./AQFO_snps_R80_bi
```

### Filtering SNPs based on depth

Repetitive or similar regions in the genome may result in reads that
come from different locations in the genome getting erroneously merged
and aligned. This can result in exceptionally deep coverage on certain
loci and could result in spurious SNP calls. We set a depth threshold
for snps as the
![99^\\text{th}](https://latex.codecogs.com/png.latex?99%5E%5Ctext%7Bth%7D
"99^\\text{th}") percentile of the empirical distribution of mean SNP
depth across samples. This is based on the quantile plot in Figure 1
where one can see that a few SNPs have extremely high coverage.

We do not remove SNPs with low depth because the model that calls
genotypes in  avoids calling genotypes when uncertainty in the call is
high (Maruki and Lynch 2017).

``` r
  dp <- extract.gt(col.vcf, element = "DP", as.numeric = T)

  # heatmap.bp(dp, rlabels = F, clabels = F)
  
  mean_dp_snps <- apply(dp, 1, mean, na.rm=T)
  
  plot(log(quantile(mean_dp_snps, probs=c(seq(0,1, by=0.01)))),
       xlab="quantile", ylab="log(depth)")
```

![](SNP_filtering_files/figure-gfm/mean%20depth-1.png)<!-- -->

**Figure 1**: Empirical CDF of the mean depth across samples for each
SNP.

``` r
# Remove snps that have extremely high average read depths
 # what value represents the 98% quantile?
  quant <- quantile(mean_dp_snps, probs=0.99)

# list the snps which have an average depth above the 98% quantile
  snps2_rmv <- which(mean_dp_snps > quant)
  
# we also want to remove snps without any neighbors
#  for the sake of imputation
  snps_per_chrom <- as.data.frame(col.vcf@fix) %>% group_by(., CHROM) %>%
    summarise(n=n())
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

``` r
  scaffs2_rmv <- snps_per_chrom$CHROM[which(snps_per_chrom$n == 1)]
  snps2_rmv <- c(snps2_rmv, 
                 which(col.vcf@fix[,1] %in% scaffs2_rmv))
  
  site2rm <- col.vcf@fix[snps2_rmv, 3]

  # write.table(site2rm, "~/Documents/Columbine/Data/GBS_Data/exclude_sites.txt",
  #             quote = F, sep = "\t",
  #           row.names = F, col.names = F)
```

![\~](https://latex.codecogs.com/png.latex?~ "~")

Using the list of SNPs to exclude, we use vcftools to rewrite the vcf
file with only the loci we wish to keep. We do not use `vcfR` to do
this, as BEAGLE does not recognize files written using `vcfR`.

![\~](https://latex.codecogs.com/png.latex?~ "~")

``` bash
  vcftools --vcf ./AQFO_snps_R80_bi.vcf --exclude ./exclude_sites.txt --recode --recode-INFO-all \
  --out ./AQFO_snps_postfilt
```

![\~](https://latex.codecogs.com/png.latex?~ "~")

With this vcf file containing 7.76 percent missing information, 11878
SNPs and 192 individual plants, we impute the missing genotypes using
BEAGLE (Browning and Browning 2016).

``` bash
# running BEAGLE through the terminal
  java -Xmx8g -jar ./beagle.25Nov19.28d.jar gt=AQFO_snps_postfilt.vcf out=AQFO_snps_impute nthreads=2
```

![\~](https://latex.codecogs.com/png.latex?~ "~")

### Filtering SNPs based on exon annotations

Read in the annotations file from the *A. coerulea* genome and store as
a dataframe.

``` r
# read in the imputed vcf
  col_vcf_gdimp <- read.vcfR(file = "~/Documents/Columbine/Data/GBS_Data/AQFO_snps_impute.vcf")
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 8
    ##   header_line: 9
    ##   variant count: 11878
    ##   column count: 201
    ## Meta line 8 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 11878
    ##   Character matrix gt cols: 201
    ##   skip: 0
    ##   nrows: 11878
    ##   row_num: 0
    ## Processed variant 1000Processed variant 2000Processed variant 3000Processed variant 4000Processed variant 5000Processed variant 6000Processed variant 7000Processed variant 8000Processed variant 9000Processed variant 10000Processed variant 11000Processed variant: 11878
    ## All variants processed

``` r
  ann <- read.table("~/Documents/Columbine/Data/GBS_Data/Acoerulea_322_v3.1.gene_exons.gff3",
                    header = F, sep = "\t", as.is = T)[,c(1,3:5,7)]

# name columns
  names(ann) <- c("Chrom", "annot", "start", "stop", "strand")
  
# structure of the data
  head(ann)
```

    ##    Chrom          annot start stop strand
    ## 1 Chr_01           gene  2657 4987      +
    ## 2 Chr_01           mRNA  2657 4987      +
    ## 3 Chr_01           exon  2657 2841      +
    ## 4 Chr_01 five_prime_UTR  2657 2841      +
    ## 5 Chr_01           exon  4435 4987      +
    ## 6 Chr_01 five_prime_UTR  4435 4439      +

``` r
# annotate snps
  snps <- as.data.frame(col_vcf_gdimp@fix)
  snps$POS <- as.numeric(snps$POS)
  snps$strand <- str_split(snps$ID, ":", simplify = T)[,3]
  snps$SOFA <- NA
  
  for(i in 1:nrow(snps)){
    annots_i <- which((ann$Chrom == snps$CHROM[i]) &
           (ann$start < snps$POS[i]) &
            (ann$stop > snps$POS[i]))
    if(length(annots_i) == 1){
      snps$SOFA[i] <- ann$annot[annots_i]
    } else if(length(annots_i) > 1){
      snps$SOFA[i] <- paste(sort(unique(ann[annots_i,]$annot)), collapse = ":")
    } else {
      snps$SOFA[i] <- "NC"
    }
  }
  
# create "neutral" vcf
  
  col_vcf_nc <- col_vcf_gdimp[which(snps$SOFA == "NC"), ]
  col_vcf_cd <- col_vcf_gdimp[which(snps$SOFA != "NC"), ]
  
  # write.vcf(col_vcf_nc, file = "~/Documents/Columbine/Data/GBS_Data/AQFO_snps_nc.vcf")
  # write.vcf(col_vcf_cd, file = "~/Documents/Columbine/Data/GBS_Data/AQFO_snps_cd.vcf")
```

![\~](https://latex.codecogs.com/png.latex?~ "~")

<div id="refs" class="references">

<div id="ref-browning2016">

Browning, Brian L., and Sharon R. Browning. 2016. “Genotype Imputation
with Millions of Reference Samples.” *The American Journal of Human
Genetics* 98 (1): 116–26. <https://doi.org/10.1016/j.ajhg.2015.11.020>.

</div>

<div id="ref-maruki2017">

Maruki, Takahiro, and Michael Lynch. 2017. “Genotype Calling from
Population-Genomic Sequencing Data.” *G3: Genes, Genomes, Genetics* 7
(5): 1393–1404. <https://doi.org/10.1534/g3.117.039008>.

</div>

</div>
