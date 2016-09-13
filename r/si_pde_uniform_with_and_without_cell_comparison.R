#
#  ---------------------------------------------------------------------------
#  3D Model: Compare SNR with and without cell to analyze effects of cell size
#  ---------------------------------------------------------------------------
#
#  Comparing the gradients and SNR values with and without cells 
#  for a range of exponential gradients, as suggested by a referee.
#  Maybe using linear gradients isn't appropriate.
#
#  --
#  Author: Igor Segota
#  Date: 2015-01-10
#
#  --
#  2015-02-17:  I noticed that some cAMP concentrations obtained
#   from COMSOL are negative (presumably it's because we're getting
#   into the regime with numerical errors). Let's check this. For 
#   example, I have -1.83e-13 c_back. The units are mol/m^3 and I know
#   that 1 nM = 1e-6

# Clear all variables
rm(list=ls(all=T))

#  Load required libraries
library(ggplot2)
library(scales)
library(directlabels)

# Set working directory
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# Load functions
source('include_cAMPPDEfunctions.R')
source('include_multiplot.R')

# Load model parameters and import into variables
params.data <- read.table('../data/model_parameters.txt', header=T, sep='\t', row.names=1)
Km <- params.data['Km',]$value
Kd <- params.data['Kd',]$value
Rt <- params.data['Rt',]$value
rc <- params.data['rc',]$value
Is <- params.data['Is',]$value
Sb <- params.data['Sb',]$value
w  <- params.data['w', ]$value
Dc <- params.data['Dc',]$value
k2 <- params.data['k2',]$value

plot_text_size = 18
snr_max = 20.1

# File names
cell.fileName = '../data/matlab_output_pde_uniform_k2_13300_s-1_Km_10_uM.txt'

# Load file with simulation results with cell
cell.data   <- read.table(cell.fileName, header=T, sep='\t')

# cell.data
# c0_left = cAMP conc on the left boundary
# c0_right = cAMP conc on the right boundary (0 here always)
# c_front  = cAMP conc at the front edge of the cell (facing the gradient)
#            x12 = +-r


#cAMP analytical result for uniform PDE (~'exponential' gradients)

# Calculate the cAMP concentrations with no cells numerically
theta0.front <- 0
theta1.front <- pi/2
thetaN.front <- 5
phi0.front   <- pi/2
phi1.front   <- 3*pi/2
phiN.front   <- 5
# now for each of these combinations, evaluate c(x)
theta.front  <- seq(theta0.front, theta1.front, length.out=thetaN.front)
phi.front    <- seq(phi0.front,   phi1.front,   length.out=phiN.front)

theta0.back <- 0
theta1.back <- pi/2
thetaN.back <- 5
phi0.back   <- -pi/2
phi1.back   <- pi/2
phiN.back   <- 5
# now for each of these combinations, evaluate c(x)
theta.back  <- seq(theta0.back, theta1.back, length.out=thetaN.back)
phi.back    <- seq(phi0.back,   phi1.back,   length.out=phiN.back)
# cAMP.back <- expand.grid(theta.back, phi.back)
#colnames(cAMP.back) <- c('theta', 'phi')

cell.length             <- length(cell.data$c0_left)
cell.pfrontNoCell       <- rep(0, cell.length)
cell.pbackNoCell        <- rep(0, cell.length)
cell.pfrontNoCellCheck  <- rep(0, cell.length)
cell.pbackNoCellCheck   <- rep(0, cell.length)

for (i in seq(1, cell.length)) {
  # for each row in cell.data, calculate c_front and c_back
  # numerically by averaging cAMPconc over theta,phi pairs.

  L <- sqrt((Km*Dc)/(k2*cell.data[i,]$p0))
  
  cell.pfrontNoCell[i] <- probRecOccAvg(
    cL=cell.data[i,]$c0_left, w=w, L=L, r=rc, Kd=Kd, theta.front, phi.front
  )
  cell.pbackNoCell[i]  <- probRecOccAvg(
    cL=cell.data[i,]$c0_left, w=w, L=L, r=rc, Kd=Kd, theta.back,  phi.back
  )
  
  cell.pfrontNoCellCheck[i] <- probRecOcc(
    cL=cell.data[i,]$c0_left, w=w, L=sqrt((Km*Dc)/(k2*cell.data[i,]$p0)), r=rc, Kd=Kd, x=-rc/4
  )
  cell.pbackNoCellCheck[i] <- probRecOcc(
    cL=cell.data[i,]$c0_left, w=w, L=sqrt((Km*Dc)/(k2*cell.data[i,]$p0)), r=rc, Kd=Kd, x=rc/4
  )
  
}

#
# I'm getting some negative values for the cbackNoCell. Something is wrong. Crosscheck this
# With the concentration at x = +- r/4 for example as a rough check to see the numbers are
# fine.
#
# Fixed 2L -> 2*L, now we don't get negative concentrations. Still, double check the
# results.
#


cell.data$pfavg_nocell <- cell.pfrontNoCell
cell.data$pbavg_nocell <- cell.pbackNoCell

cell.data$pfavg_nocell_simple <- cell.pfrontNoCellCheck
cell.data$pbavg_nocell_simple  <- cell.pbackNoCellCheck

# Now find the SNR and SNR_nocell and compare them in a graph for a range of parameters
# cL and P0 (or their fraction) to check what's the prefactor on SNR due to cell size.

# Calculate 3 SNRs (first for our simulation with the cell, second for the analytical
# result with no cell but assuming exponential gradient and integrated over the same
# surface, third for the case where we just find concentration difference between
# x=-r/4 and x=r/4. 

cell.data$SNRwithCell <- SNR(
  pf=cell.data$pfavg, pb=cell.data$pbavg, cMean=cell.data$c_mean, Kd=Kd, I=Is, Sb=Sb, Rt=Rt
)

L = function (p0) sqrt((Km*Dc)/(k2*p0))
snr_exact <- function(cleft, pde0) {
    l = L(pde0)
    cfront <- cleft*sinh(0.5*w/l - 0.5*rc/l)/sinh(w/l)
    cback  <- cleft*sinh(0.5*w/l + 0.5*rc/l)/sinh(w/l)
    Rfront <- (Rt/2)*(cfront/(cfront + Kd))
    Rback  <- (Rt/2)*(cback /(cback  + Kd))
    sigmaRfSq <- (Rt/2)*(cfront*Kd/(cfront+Kd)^2)
    sigmaRbSq <- (Rt/2)*(cback *Kd/(cback +Kd)^2)
    snr <- (Rfront - Rback) / sqrt(sigmaRfSq/Is + sigmaRbSq/Is + Sb^2)
    return(-snr)
}


cell.data$snr_no_cell = snr_exact(cell.data$c0_left, cell.data$p0)

#cell.data$SNRnoCellAnalytic <- SNR(
#  pf=cell.data$pfavg_nocell, pb=cell.data$pbavg_nocell, cMean=cell.data$c_mean, Kd=Kd, I=Is, 
#  Sb=Sb, Rt=Rt
#)

#cell.data$SNRnoCellAnalyticSimpleAverage <- SNR(
#  pf=cell.data$pfavg_nocell_simple, pb=cell.data$pbavg_nocell_simple, cMean=cell.data$c_mean, Kd=Kd, I=Is, 
#  Sb=Sb, Rt=Rt
#)


cell.plotdata <- cell.data[, c('c0_left', 'p0', 'SNRwithCell', 'snr_no_cell')]

#
# Plot SNR with cell
#

x_range = -5:1
y_range = -2:5

snr_with_cell.plot = ggplot(cell.plotdata, aes(x=p0*1e6, y=c0_left/Kd, fill=SNRwithCell, z=SNRwithCell)) +
    geom_raster(show.legend = T) + theme_basic() +
    stat_contour(mapping=aes(color=..level..), color='black', size=0.7, breaks=c(1, 2, 5)) +
    annotation_logticks(color='white') +
    scale_x_continuous(
        name = expression(paste(p[0], ' [nM]')), 
        breaks = 10^x_range, trans = 'log10', 
        labels = sapply(x_range, function (i) as.expression(bquote(10^ .(i)))), 
        expand = c(0,0)
    ) +
    scale_y_continuous(
        name = expression(paste('c(-w/2) [', K[d], ']')), 
        breaks = 10^y_range, trans = 'log10',
        labels = sapply(y_range, function (i) as.expression(bquote(10^ .(i)))), 
        expand = c(0,0)
    ) +
    scale_fill_gradientn('SNR with cell',
        trans = 'sqrt',
        limits = c(0, snr_max),
        breaks = c(0, 1, 4, 9, 16, 25),
        colours = rev(c(
          '#661f19', '#8c2a22', '#e9622e', '#f7e92b', 
          '#86c98e', '#00b2ed', '#2061ae', '#303592')),
        guide = 'colorbar'
    ) +
    scale_color_gradient(guide = 'none') +
    coord_cartesian(xlim=c(1e-5, 1)) +
    theme(
        text = element_text(size=plot_text_size, color='black'),
        axis.text = element_text(size=plot_text_size, color='black'),
        axis.title = element_text(size=plot_text_size, color='black'),
        legend.text = element_text(size=0.7*plot_text_size, color='white'),
        legend.title = element_text(size=0.7*plot_text_size, color='white'),
        plot.margin = margin(10, 10, 2, 2, unit='mm'),
        legend.key.size = unit(10, 'mm'),
        legend.position = c(0.04, 1),
        legend.background = element_rect(fill=alpha('black', 0)),
        legend.justification = c(0, 1)
    )
snr_with_cell.plot = direct.label(snr_with_cell.plot, method=list('bottom.pieces', cex=2, colour='black'))
# snr_with_cell.plot


# Plot SNR with no cell for the same values of p0 and c0_left
snr_no_cell.plot = ggplot(cell.plotdata, aes(x=p0*1e6, y=c0_left/Kd, fill=snr_no_cell, z=snr_no_cell)) +
    geom_raster(show.legend = T) + theme_basic() +
    stat_contour(mapping=aes(color=..level..), color='black', size=0.7, breaks=c(1, 2, 5)) +
    annotation_logticks(color='white') +
    scale_x_continuous(
        name = expression(paste(p[0], ' [nM]')), 
        breaks=10^x_range, trans = 'log10',
        labels=sapply(x_range, function (i) as.expression(bquote(10^ .(i)))), 
        expand=c(0,0)
    ) +
    scale_y_continuous(
        name = expression(paste('c(-w/2) [', K[d], ']')), 
        breaks = 10^y_range, trans = 'log10',
        labels = sapply(y_range, function (i) as.expression(bquote(10^ .(i)))), 
        expand = c(0,0)
    ) +
    scale_fill_gradientn(
        name = 'SNR without cell',
        trans = 'sqrt',
        limits = c(0, snr_max),
        breaks = c(0, 1, 4, 9, 16, 25),
        colours = rev(c(
            '#661f19', '#8c2a22', '#e9622e', '#f7e92b', 
            '#86c98e', '#00b2ed', '#2061ae', '#303592')),
        guide = 'colorbar'
    ) +
    scale_color_gradient(guide = 'none') +
    coord_cartesian(xlim=c(1e-5, 1)) +
    theme(
        text = element_text(size=plot_text_size, color='black'),
        axis.text = element_text(size=plot_text_size, color='black'),
        axis.title = element_text(size=plot_text_size, color='black'),
        legend.text = element_text(size=0.7*plot_text_size, color='white'),
        legend.title = element_text(size=0.7*plot_text_size, color='white'),
        plot.margin = margin(10, 10, 2, 2, unit='mm'),
        legend.key.size = unit(10, 'mm'),
        legend.position = c(0.04, 1),
        legend.background = element_rect(fill=alpha('black', 0)),
        legend.justification = c(0, 1)
    )
snr_no_cell.plot = direct.label(snr_no_cell.plot, method=list('bottom.pieces', cex=2, colour='black'))
# snr_no_cell.plot


# Plot ratio of SNR with cell / SNR without cell
cell.plotdata$snr_with_over_without_cell = cell.plotdata$SNRwithCell/cell.plotdata$snr_no_cell

snr_cell_ratio.plot = ggplot(cell.plotdata[cell.plotdata$snr_no_cell > 1e-3,], aes(x=p0*1e6, y=c0_left/Kd, fill=snr_with_over_without_cell, z=SNRwithCell)) +
    geom_raster(show.legend = T) + theme_basic() +
    # stat_contour(mapping=aes(color=..level..), color='black', size=0.7, breaks=c(1, 2, 5), linetype='dashed') +
    annotation_logticks(color='white') +
    scale_x_continuous(
        name = expression(paste(p[0], ' [nM]')), 
        breaks=10^x_range, trans = 'log10',
        labels=sapply(x_range, function (i) as.expression(bquote(10^ .(i)))), 
        expand=c(0,0)
    ) +
    scale_y_continuous(
        name = expression(paste('c(-w/2) [', K[d], ']')), 
        breaks = 10^y_range, trans = 'log10',
        labels = sapply(y_range, function (i) as.expression(bquote(10^ .(i)))), 
        expand = c(0,0)
    ) +
    scale_fill_gradientn(
        name = expression(over('SNR'['with cell'], 'SNR'['without cell'])),
        trans = 'log10',
        limits = c(0.8, 141),
        colours = rev(c(
             '#661f19', '#8c2a22', '#e9622e', '#f7e92b', 
             '#86c98e', '#00b2ed', '#2061ae', '#303592')),
        guide = 'colorbar',
        na.value = 'gray'
    ) +
    scale_color_gradient(guide = 'none') +
    coord_cartesian(xlim=c(1e-5, 1)) +
    theme(
        text = element_text(size=plot_text_size, color='black'),
        axis.text = element_text(size=plot_text_size, color='black'),
        axis.title = element_text(size=plot_text_size, color='black'),
        legend.text = element_text(size=0.7*plot_text_size, color='black'),
        legend.title = element_text(size=0.7*plot_text_size, color='black'),
        plot.margin = margin(10, 10, 2, 2, unit='mm'),
        legend.key.size = unit(10, 'mm'),
        panel.background = element_rect(fill='gray'),
        legend.position = c(0, 1),
        legend.background = element_rect(fill='white', color='black'),
        legend.justification = c(0, 1),
        legend.margin = unit(2, 'mm')
    )
# snr_cell_ratio.plot = direct.label(snr_cell_ratio.plot, method=list('bottom.pieces', cex=2, colour='black'))
# snr_cell_ratio.plot


# Plot all three
pdf('../figures/si_pde_uniform_compare_with_without_cell.pdf', width=15.1, height=4.5, useDingbats=F)
multiplot(snr_with_cell.plot, snr_no_cell.plot, snr_cell_ratio.plot, cols = 3, annotate = T, hjust=-0.26, vjust=-0.5)
dev.off()

# Export 4.0 x 3.0 in
