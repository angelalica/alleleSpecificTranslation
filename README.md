###init README.md

THIS SCRIPT determines if two alleles are significantly differentially translated given ddPCR data.

There are two required arguments:<br/>
1) ```--filename``` = the name of the input ddPCR file<br/>
2) ```--expected_ratios``` = expected mutant/wildtype ratio (usually 1:1 for a diploid heterozygote)<br/>

Optional arguments include:<br/>
1) ```--fractions```: a vector that specifies whether the groups analyzed will be the individual fractions or groups based on<br/> number of ribosomes<br/>
    Default value = NULL (groups generated will correspond to each individual polysome fraction)<br/>
    If vector is specified (i.e. ```7,3,1,1,1,5```), then the individual fractions will be further grouped based on number of ribosomes<br/>

2) ```--diff```: a numeric value that specifies the maximum acceptable difference between inputted and calculated<br/> mutant/wildtype expected ratios
    Default value = 0.10 (10%)<br/>

3) ```--nsims```: the number of boostrap simulations performed <br/>
    Default value = 1000<br/>

4) ```--output```: filename of the output file (where the generated pvalues are)<br/>
    Default value = "output.csv", located in the same directory as this script<br/>

