#
#  Generate a simulation of receptor occupancy for some 
#  numbers of cell receptors.
#

# Clear all variables
rm(list=ls(all=TRUE))

# Load plotting packages
library(ggplot2)
library(data.table)

# Set working directory
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# Load model parameters and import into variables
params.data = read.table('../data/model_parameters.txt', 
                         header=T, sep='\t', row.names=1)
Km = params.data['Km',]$value
Kd = params.data['Kd',]$value
Rt = params.data['Rt',]$value
rc = params.data['rc',]$value
Is = params.data['Is',]$value
Sb = params.data['Sb',]$value
w  = params.data['w', ]$value
Dc = params.data['Dc',]$value
k2 = params.data['k2',]$value

# Override otal number of receptors per cell
N = 70000

# Generate x and y coordinates
# circle: x^2 + y^2 <= r^2
x  = seq(-1, 1, length.out = sqrt(N)*4/pi)
dt = data.table(expand.grid(x, x))
names(dt) = c('x', 'y')

# Gradient parameters (exponential gradient)
cL = 3
l  = 0.5

# Generate concentrations
dt[, c := cL*exp(-(x+l)/l)]

# Generate receptor occupancies
dt[, rec := rbinom(length(c), 1, c/(c+1))]
dt[x^2 + y^2 > 1, rec := NA] 

# Plot receptors
rec.plot = ggplot(dt, aes(x=x, y=y, fill=factor(rec))) +
    geom_raster() + coord_fixed() +
    scale_fill_manual(values=c('0' = 'gray', '1' = 'black')) +
    theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none'
    )

# Print figure to PDF
pdf('../figures/fig1_receptor_distribution.pdf', width=5, height=5, useDingbats=F)
rec.plot
dev.off()