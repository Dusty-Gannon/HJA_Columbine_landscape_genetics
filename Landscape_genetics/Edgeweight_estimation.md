Estimating edge weights connecting HJA columbine genetic networks
================
D. G. Gannon
October 2020

## Data

We have SNP data on
![n=192](https://latex.codecogs.com/png.latex?n%3D192 "n=192")
*Aquilegia formosa* individuals from XX meadows, which we consider to be
“sub-populations” of the H.J. Andrews *A. formosa* population. We work
with a centered, euclidean distance matrix, ![{\\bf D}\_{(n\\times
n)}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20D%7D_%7B%28n%5Ctimes%20n%29%7D
"{\\bf D}_{(n\\times n)}"), computed from a genetic loading matrix
![{\\bf L}\_{(n\\times
\\ell)}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20L%7D_%7B%28n%5Ctimes%20%5Cell%29%7D
"{\\bf L}_{(n\\times \\ell)}"), where
![\\ell](https://latex.codecogs.com/png.latex?%5Cell "\\ell") is the
number of SNP loci sequenced. Genetic loadings are coded as ![{\\bf
L}\_{ik}
= 0](https://latex.codecogs.com/png.latex?%7B%5Cbf%20L%7D_%7Bik%7D%20%3D%200
"{\\bf L}_{ik} = 0") if individual
![i](https://latex.codecogs.com/png.latex?i "i") is homozygous for the
reference allele at locus ![k](https://latex.codecogs.com/png.latex?k
"k"), ![{\\bf
L}\_{ik}=1](https://latex.codecogs.com/png.latex?%7B%5Cbf%20L%7D_%7Bik%7D%3D1
"{\\bf L}_{ik}=1") if individual
![i](https://latex.codecogs.com/png.latex?i "i") has a heterozygous
genotype at site ![k](https://latex.codecogs.com/png.latex?k "k"), and
![{\\bf
L}\_{ik}=2](https://latex.codecogs.com/png.latex?%7B%5Cbf%20L%7D_%7Bik%7D%3D2
"{\\bf L}_{ik}=2") if individual
![i](https://latex.codecogs.com/png.latex?i "i") is homozygous with two
copies of the alternative allele at site
![k](https://latex.codecogs.com/png.latex?k "k"). The data were
previously filtered to include only bi-allelic loci.

Letting ![\\tilde{\\bf
L}](https://latex.codecogs.com/png.latex?%5Ctilde%7B%5Cbf%20L%7D
"\\tilde{\\bf L}") be the column-centered loadings matrix, we compute
![{\\bf D}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20D%7D
"{\\bf D}") as

  
![&#10;{\\bf D} = \\text{diag}({\\bf G}){\\bf 1}^\\text{T} - 2{\\bf G} +
\\text{diag}({\\bf
G}){\\bf 1}^{\\text{T}},&#10;](https://latex.codecogs.com/png.latex?%0A%7B%5Cbf%20D%7D%20%3D%20%5Ctext%7Bdiag%7D%28%7B%5Cbf%20G%7D%29%7B%5Cbf%201%7D%5E%5Ctext%7BT%7D%20-%202%7B%5Cbf%20G%7D%20%2B%20%5Ctext%7Bdiag%7D%28%7B%5Cbf%20G%7D%29%7B%5Cbf%201%7D%5E%7B%5Ctext%7BT%7D%7D%2C%0A
"
{\\bf D} = \\text{diag}({\\bf G}){\\bf 1}^\\text{T} - 2{\\bf G} + \\text{diag}({\\bf G}){\\bf 1}^{\\text{T}},
")  

where ![{\\bf G} = \\tilde{\\bf L}^\\text{T} \\tilde{\\bf
L}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20G%7D%20%3D%20%5Ctilde%7B%5Cbf%20L%7D%5E%5Ctext%7BT%7D%20%5Ctilde%7B%5Cbf%20L%7D
"{\\bf G} = \\tilde{\\bf L}^\\text{T} \\tilde{\\bf L}") and
![{\\bf 1}](https://latex.codecogs.com/png.latex?%7B%5Cbf%201%7D
"{\\bf 1}") is an ![n](https://latex.codecogs.com/png.latex?n
"n")-dimensional vector of ones.

## Model

We assume that ![{\\bf
D}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20D%7D "{\\bf D}")
