#
# Fixed PDE concentration, spherical model (point source of cAMP)
#


# Clear all variables
rm(list=ls(all=T))

# Load necessary packages
library(ggplot2)
library(data.table)
library(directlabels)

# Set working directory
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

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

# Plotting parameters
plot_text_size = 22

# Define functions

L = function (p0) sqrt((Km*Dc)/(k2*p0))
c_spherical = function (C0, r, p0) C0*exp(-r/L(p0))/r
p = function (c) c/(c+Kd)

snr_spherical_exact = function (C0, r, p0) {
    c_front = c_spherical(C0, r-0.5*rc, p0)
    c_back = c_spherical(C0, r+0.5*rc, p0)
    delta_r = 0.5*Rt*(p(c_front) - p(c_back))
    sigma_dr2 = 0.5*Rt*(p(c_front)*(1-p(c_front)) + p(c_back)*(1-p(c_back)))
    sigma = sqrt(Sb^2 + sigma_dr2)
    delta_r/sigma
}

# Generate data frame with this functionj
c0_range = 10^seq(-10, -4, length.out=50)    # units: mol/m^2
p0_range = 10^seq(-11, -6, length.out=50)    # units: mol/(m^2 s)
# r0_range = R*seq(2, 200, length.out=20)
r0_range = c(10, seq(75, 1150, 75))*1e-6

snr.dt = data.table(expand.grid(c0_range, p0_range, r0_range))
names(snr.dt) = c('C0', 'p0', 'r0')

snr.dt[, snr := snr_spherical_exact(C0, r0, p0), by=names(snr.dt)]
snr.dt[, r0_labels := paste0('paste(r == ', round(r0/1e-6), ', " ", mu, m)')]
snr.dt[, r0_labels := factor(r0_labels, levels=unique(r0_labels))]

# snr.df$r0_labels = paste0('paste(r[0] == ', round(snr.df$r0/1e-6), ', \' \', mu, m)')
# snr.df = snr.df[order(snr.df$r0),]
# snr.df$r0_labels = factor(snr.df$r0_labels, levels=unique(snr.df$r0_labels))

# Plot SNR
x_labels = sapply(p0_range, function(p0) as.expression(p0_range))
x_range = -5:-1
y_range = -10:-4
r0_ranges = r0_range[seq(1, length(r0_range), 1)]
# plot_text_size = 30

snr.plot = ggplot(
        snr.dt[r0 %in% r0_ranges,], 
        aes(x=p0*1e6, y=C0, fill=snr, z=snr)
    ) +
    geom_raster(show.legend = T) +
    stat_contour(mapping=aes(color=..level..), color='black', size=0.7, breaks=c(1, 2, 5, 9)) +
    scale_color_gradient(guide = 'none') +
    facet_wrap( ~ r0_labels, ncol=4, labeller = label_parsed) +
    scale_fill_gradientn(
        name = 'SNR', trans = 'sqrt',
        colours = rev(c('#000000', '#661f19', '#8c2a22', '#e9622e', '#f7e92b', 
                        '#86c98e', '#00b2ed', '#2061ae', '#303592')),
        limits = c(0, 52), breaks=(0:13)^2, 
        values = c((0:7)/12, 1),
        guide = 'colourbar'
    ) +
    scale_x_continuous('PDE concentration [nM]', trans='log10', expand=c(0,0),
                       breaks=10^x_range, 
                       labels=sapply(x_range, function (i) as.expression(bquote(10^ .(i))))
                           ) +
    scale_y_continuous(expression(paste(C[0], ' [mol ', m^-2, ']')), trans='log10', expand=c(0,0), 
                       breaks=10^y_range,
                       labels=sapply(y_range, function (i) as.expression(bquote(10^ .(i))))) +
    annotation_logticks(sides='bl', color='white') +
    theme(
        text = element_text(size=plot_text_size, color='black'),
        axis.text = element_text(size=plot_text_size, color='black'),
        axis.ticks = element_blank(),
        axis.title = element_text(size=plot_text_size, color='black'),
        legend.key.size = unit(12, 'mm'),
        legend.key.height = unit(40, 'mm'),
        legend.position = c(1, 0.5),
        legend.justification = c(0, 0.5),
        legend.background = element_rect(fill=alpha('black', 0)),
        legend.text = element_text(size=plot_text_size, color='black'),
        legend.title = element_text(size=plot_text_size, color='black'),
        plot.margin = margin(5, 40, 2, 2, unit='mm'),
        panel.border = element_rect(colour = 'black', fill=alpha('blue',0))
    )  
snr.plot = direct.label(snr.plot, method=list('bottom.pieces', cex=2, colour='black'))

# Generate figure
pdf('../figures/si_pde_uniform_spherical_model_snr.pdf', useDingbats=F, width=15, height=15)
snr.plot
dev.off()
