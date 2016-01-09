#!/usr/bin/Rscript
install.packages("argparse", repos="http://cran.rstudio.com/")
require("argparse")

### THIS SCRIPT determines if two alleles are significantly differentially translated given ddPCR data

# Required args: 
# filename = name of the input file containing the ddPCR data. The ddPCR data is a matrix where each row is a sample, and the 4 columns contain the Ch1 and Ch2 pos/neg counts (where Ch1 and Ch2 correspond to the mutant and wildtype allele)
# expected_ratios = user's expected mutant/wildtype ratio of counts (usually 1:1 if heterozygous with 2 alleles) = usually generated from totalRNA ddPCR sample

# Optional arguments: fractions, diff, nsims, output (see help for more information)
parser <- ArgumentParser(description='Process cmd line arguments')

parser$add_argument('--filename', default = "~/Documents", help='name of input ddPCR file')
parser$add_argument("--expected_ratios", default = 1.0, type = "double", help="expected mutant/wildtype ratio")
parser$add_argument("--fractions", default = '', help = "vector containing grouped fractions (inferred from polysome profiling data)")
parser$add_argument("--diff", default = 0.10, type = 'double', help = "maximum acceptable difference between inputted and calculated mutant/wildtype expected ratios")
parser$add_argument("--nsims", default = 1000, type = 'integer', help = "number of bootstrap simulations")
parser$add_argument("--output", default = "output.csv", help = "filename of output file with pvalues")

args <- parser$parse_args()

filename = args$filename
expected_ratios = args$expected_ratios
fractions = as.numeric(unlist(strsplit(args$fractions, ",")))
diff = args$diff
nsims = args$nsims
output = args$output

#INPUT: name of file to read in
#OUTPUT: condensed data (mean of each fraction, summarizing the replicates)
summarize_data = function(filename) {
  
  full_data = read.table(filename, header = T, sep = ",")
  
  #add the ones where there's positive in both channels to each respective channel
  for (i in 1:length(rownames(full_data))) {
    full_data[i,3] = full_data[i,3] + full_data[i,2]
    full_data[i,4] = full_data[i,4] + full_data[i,2]
    full_data[i,6] = c(full_data[i,3] / full_data[i,4])
  }
  colnames(full_data) = c("fraction", "both", "mutant", "wildtype", "neither", "MT/WT ratio")
  
  #initialize empty data frame
  condensed_data = NULL
  #condense by group id   
  for (fraction in unique(full_data[, 1])) {
    condensed_row = c(mean(full_data[which(full_data[, 1] == fraction), 3]), 
                      mean(full_data[which(full_data[, 1] == fraction), 4]),
                      mean(full_data[which(full_data[, 1] == fraction), 5]),
                      mean(full_data[which(full_data[, 1] == fraction), 6]))
    rbind(condensed_data, condensed_row) -> condensed_data
  }
  rownames(condensed_data) = unique(full_data[, 1])
  colnames(condensed_data) = c("mutant", "wildtype", "neither", "ratio")
  
  return (as.matrix(condensed_data))
}

#INPUT: condensed data (as shown above), expected MT/WT ratio inputted by the user (i.e. totalRNA generated from ddPCR data)
#OUTPUT: value of expected ratio to be used for the remainder of the program (either inputted or calculated)
#function calculates the expected counts from summing up the data and outputs both expected ratios, then asks user which one to use moving forward
get_expected_counts = function(condensed_data, expected_ratio, diff) {
  
  #sum up all values to obtain totalRNA MT/WT ratio as calculated from ddPCR data
  total_calculated_ratio = sum(condensed_data[, 1]) / sum(condensed_data[, 2])
  cat(sprintf("Inputted Expected Ratio: %f\n", expected_ratio))
  cat(sprintf("Calculated Expected Ratio: %f\n", total_calculated_ratio))
  
  return (expected_ratios)
#   if ((abs(total_calculated_ratio - expected_ratio) > (expected_ratio * diff)) |
#       (abs(total_calculated_ratio - expected_ratio) > (total_calculated_ratio * diff))) {
#     cat(sprintf("Warning: the inputted expected MT/WT ratio is more than %d percent away from the calculated MT/WT ratio.\n", accepted_diff*100))
#     print("Would you like to proceed? Please enter 'y' or 'n'.\n")
#     if (readLines(con = stdin(), n = 1) == "n") {
#       stop("Now exiting program.")
#     }
#   }
#   
#   cat("Press '1' to use the inputted ratio as the expected ratio, press '2' to use the calculated ratio as the expected ratio:\n")
#   whichRatio = readLines(file("stdin"), n = 1)
#   if (whichRatio == "1") { 
#     return (expected_ratio)
#   } else if (whichRatio == "2") { 
#     return (total_calculated_ratio) 
#   }   
}

#INPUT: condensed_data, vector containing how many fractions correspond to the number of ribosomes/mRNA 
#OUTPUT: condensed_data further summarized by fractions corresponding to ribosomal density (free, 1 ribosome/mRNA, etc.)
group_data = function(condensed_data, grouped_fractions) {
  grouped_data = NULL
  for (i in 1:length(grouped_fractions)) {
    if (grouped_fractions[i] > 1) {
      temp = apply(as.matrix(condensed_data)[1:grouped_fractions[i], ], 2, FUN = mean)
    } else { temp = condensed_data[1, ] }
    condensed_data = condensed_data[-(1:grouped_fractions[i]), ]
    grouped_data = rbind(grouped_data, temp)
  }
  rownames(grouped_data) = c(1:(length(grouped_fractions)))
  return(grouped_data)
}

#INPUT: input vector has Ch1+ , Ch2+, and Ch1-Ch2- counts
#OUTPUT: We output bootstrap ratios (MT/WT ratios)
bootstrap_ratio = function (x, nsims) { 
  observation = c(rep(1, x[1]), rep(2, x[2]), rep(0, x[3]) )
  ratios = c()
  for (j in 1:nsims) { 
    boot_tmp = sample(observation, replace =T)
    ratios = c(ratios , length(which(boot_tmp == 1)) / length(which(boot_tmp == 2))  )
  }
  return (ratios)
}

#INPUT: summarized data (grouped either by replicates or by ribosomal density), expected_value (as decided by user in the get_expected_value function), and number of simulations for bootstrap to run)
#OUTPUT: graphically displays mean and distribution of boostrapped ratios for each group (either polysome fraction number or ribosome density)
#writes all pvalues for each group to an output file
is_significant = function (summarized_data, expected_value, n) {
  
  num_iterations = dim(summarized_data)[1]
  
  ratio_matrix = NULL
  all_pvalues = NULL
  
  for (i in 1:num_iterations) { #generate pvalues for each group
    ratios = bootstrap_ratio(summarized_data[i, 1:3], n) 
    ratio_matrix = cbind(ratio_matrix, ratios)
    
    num_extreme = ratios[which(ratios <= expected_value)]
    pvalue = length(num_extreme) / n
    if (pvalue > 0.5) { pvalue = 1 - pvalue }
    pvalue = pvalue * 2 #double sided test
    
    all_pvalues = append(all_pvalues, pvalue)
  }
  
  #plot the values
  x = c(1:num_iterations)
  colnames(ratio_matrix) = x
  y = apply(ratio_matrix, 2, mean)
  sd = apply(ratio_matrix, 2, sd)
  epsilon = 0.06 #arbitrary width of the tops of the error bars on the graph
  
  pdf("group_graphs.pdf")
  
  plot(x, y, 
       xaxt="n",
       ylim = c(0, max(y) + 2*max(sd)), 
       xlab = "Groups",
       ylab = "Bootstrapped MT/WT Ratios", 
       main = "Mean and SDs of Bootstrapped MT/WT Ratios")
  
  axis(1, 
       at = 1:length(rownames(summarized_data)), 
       labels = rownames(summarized_data), 
       cex.axis = 0.5,
       las=2)
  par(mar = c(6, 4.1, 4.1, 2.1), mgp = c(3, 0.5, 0))
  
  segments(x, y-sd, x, y+sd)
  segments(x-epsilon,y-sd,x+epsilon,y-sd)
  segments(x-epsilon,y+sd,x+epsilon,y+sd)
  
  segments(x, y-2*sd, x, y+2*sd)
  segments(x-epsilon,y-2*sd,x+epsilon,y-2*sd)
  segments(x-epsilon,y+2*sd,x+epsilon,y+2*sd)
  
  abline(a = expected_ratios, b = 0)
  names(all_pvalues) = rownames(summarized_data)
  
  dev.off()
  return (all_pvalues)
  
}

#ACTUAL PROGRAM:

condensed_data = summarize_data(args$filename)
expected = get_expected_counts(condensed_data, expected_ratios, diff)
if (!is.null(fractions)) {
  data = group_data(condensed_data, fractions) 
} else { 
  data = condensed_data 
}

pvalues = is_significant(data, expected, nsims) 
write.table(pvalues, output, row.names = TRUE)
