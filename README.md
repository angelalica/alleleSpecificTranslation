###init README.md

THIS SCRIPT determines if two alleles are significantly differentially translated given ddPCR data.

There are two required arguments:
1) ```--filename``` = the name of the input ddPCR file
2) ```--expected_ratios``` = expected mutant/wildtype ratio (usually 1:1 for a diploid heterozygote)

Optional arguments include:
1) ```--fractions```: a vector that specifies whether the groups analyzed will be the individual fractions or groups based on number of ribosomes
    Default value = NULL (groups generated will correspond to each individual polysome fraction)
    If vector is specified (i.e. ```7,3,1,1,1,5```), then the individual fractions will be further grouped based on number of ribosomes

2) ```--diff```: a numeric value that specifies the maximum acceptable difference between inputted and calculated mutant/wildtype expected ratios
    Default value = 0.10 (10%)

3) ```--nsims```: the number of boostrap simulations performed 
    Default value = 1000

4) ```--output```: filename of the output file (where the generated pvalues are)
    Default value = "output.csv", located in the same directory as this script

