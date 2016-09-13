#
# Load ci_dcdx_c0_from_experiments.txt data table and export it as a LaTeX tabular.
#

# Load libraries
library(data.table)

# Set working directory
this.dir = dirname(parent.frame(2)$ofile)
setwd(this.dir)

# Include export_latex() function
source('include_export_latex.R')

# Load data
dt = fread('../data/ci_dcdx_c0_from_experiments.txt')

# Write LaTeX table (remember to fix header by putting latex equations)
export_table('../latex/table_ci_snr_others_.tex')