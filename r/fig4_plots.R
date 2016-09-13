#
# Figure 4 plots. Compare Fixed PDE flux model to Uniform PDE without/with cell 
# and spherical model.
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

# Load multiplot
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

# Plotting parameters
plot_text_size = 25
snr_max = 20.1
legend_text_size = 19

# Define functions
p = function (c) c/(c+Kd)
L = function (p0) sqrt((Km*Dc)/(k2*p0))
c_sph = function (C0, r, p0) C0*exp(-r/L(p0))/r
c_con = function (x, p0, cL) cL*sinh(0.5*w/L(p0) - x/L(p0))/sinh(w/L(p0))


# Generate Fixed PDE concentration in microfluidic geometry datas
# c0 = c(-w/2)/2, not the mean concentration at x=0.
c0s = 10^seq(-2, 5, 0.1)*Kd
p0s = 10^seq(-11, -6, length.out = 100)
snr.dt = data.table(expand.grid(c0s, p0s))
names(snr.dt) = c('c0', 'p0')
snr.dt[, c_front := c_con(-0.5*rc, p0, c0*2), by=names(snr.dt)]
snr.dt[, c_back := c_con(0.5*rc, p0, c0*2), by=names(snr.dt)]
snr.dt[, delta_r := 0.5*Rt*(p(c_front) - p(c_back))]
snr.dt[, sigma_delta_r := sqrt(Sb^2 + 0.5*Rt*p(c_front)*(1-p(c_front))/Is + 0.5*Rt*p(c_back)*(1-p(c_back))/Is)]
snr.dt[, snr := delta_r/sigma_delta_r]
snr.dt[, c0_over_kd := c0/Kd]

# Generate spherical model with uniform PDE concentration data
c0_range = 10^seq(-10, -4, length.out=100)    # units: mol/m^2
p0_range = 10^seq(-11, -6, length.out=100)    # units: mol/(m^2 s)
r0_range = c(10, seq(75, 1150, 75))*1e-6
sph_snr.dt = data.table(expand.grid(c0_range, p0_range, r0_range))
names(sph_snr.dt) = c('C0', 'p0', 'r0')
sph_snr.dt[, c_front := c_sph(C0, r0 - 0.5*rc, p0)]
sph_snr.dt[, c_back  := c_sph(C0, r0 + 0.5*rc, p0)]
sph_snr.dt[, delta_r := 0.5*Rt*(p(c_front) - p(c_back))]
sph_snr.dt[, sigma_delta_r := sqrt(Sb^2 + 0.5*Rt*p(c_front)*(1-p(c_front))/Is + 0.5*Rt*p(c_back)*(1-p(c_back))/Is)]
sph_snr.dt[, snr := delta_r/sigma_delta_r]
sph_snr.dt[, r0_labels := paste0('paste(r[0] == ', round(r0/1e-6), ', " ", mu, m)')]
sph_snr.dt[, r0_labels := factor(r0_labels, levels=unique(r0_labels))]

# Plot Fixed PDE concentration SNR
x_range = -5:0
y_range = -2:5
snr.plot = ggplot(snr.dt, aes(y = 2*c0_over_kd, x = p0*1e6, z = snr, fill = snr)) +
    geom_raster(show.legend = T) +
    stat_contour(mapping=aes(color=..level..), color='black', size=0.7, breaks=c(1, 2, 5)) +
    scale_y_continuous(
        name = expression(paste('c(-w/2) [', K[d], ']', sep='')), trans='log10', expand = c(0, 0),
        breaks = 10^y_range,
        labels = sapply(y_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    scale_x_continuous(
        name = expression(paste(p[0], ' [nM]')), trans='log10', expand = c(0, 0),
        breaks = 10^x_range,
        labels = sapply(x_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    scale_fill_gradientn(
        name = 'SNR', trans = 'sqrt',
        colours = rev(c('#661f19', '#8c2a22', '#e9622e', '#f7e92b', 
                        '#86c98e', '#00b2ed', '#2061ae', '#303592')),
        breaks = c(0, (1:round(sqrt(snr_max)))^2),
        limits = c(0, snr_max),
        guide = 'colourbar'
    ) +
    scale_color_gradient(guide = 'none') +
    annotation_logticks(sides='bl', color='white') +
    theme(
        text = element_text(size=plot_text_size, color='black'),
        axis.text = element_text(size=plot_text_size, color='black'),
        axis.title = element_text(size=plot_text_size, color='black'),
        axis.ticks = element_blank(),
        plot.margin = margin(5, 10, 2, 6, unit='mm'), # t r b l
        legend.position = c(0, 1),
        legend.text = element_text(size=legend_text_size, color='white'),
        legend.title = element_text(size=legend_text_size, color='white'),
        legend.background = element_rect(fill=alpha('black', 0)),
        legend.justification = c(0, 1),
        legend.key.size = unit(13, 'mm'),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = 'black', fill=alpha('blue',0))        
    )
snr.plot = direct.label(snr.plot, method=list('bottom.pieces', cex=2, colour='black'))

# SNR Fixed PDE concentration, horizonstal slices
c0_plot = snr.dt[, sort(unique(c0_over_kd))]
c0_plot = c0_plot[seq(8, 71, by=5)]
# c0_plot = c0_plot[seq(1, length(c0_plot), 5)]
y_range = -3:1
snr_hor.plot = ggplot(snr.dt[c0_over_kd %in% c0_plot], aes(x=p0*1e6, y=snr, group=factor(2*c0_over_kd), color=factor(2*c0_over_kd))) +
    geom_line(size=1.5) +
    scale_x_continuous(
        name = expression(paste(p[0], ' [nM]', sep='')), trans='log10',
        breaks = 10^x_range, expand = c(0, 0),
        labels = sapply(x_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    scale_y_continuous(
        name = 'SNR',
        trans = 'log10', 
        breaks = 10^y_range,
        labels = sapply(y_range, function (i) as.expression(bquote(10^ .(i)))),
        expand = c(0, 0)
    ) +
    coord_cartesian(ylim=c(1e-2, 3e2)) +
    scale_color_manual(
        name = expression(paste('c(-w/2) [', K[d], ']', sep='')),
        values = rainbow(length(c0_plot), end=0.8),
        labels = sapply(round(log10(c0_plot*2), 1), function (i) as.expression(bquote(10^ .(i)))),
        breaks = c0_plot*2
    ) +
    annotation_logticks(sides='bl') +
    guides(color=guide_legend(ncol=5, byrow=F, override.aes = list(size=10))) +
    theme(
        text = element_text(size=plot_text_size, color='black'),
        axis.text = element_text(size=plot_text_size, color='black'),
        axis.title = element_text(size=plot_text_size, color='black'),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color='gray89', linetype='dashed'),
        panel.border = element_rect(colour = 'black', fill=alpha('blue',0)),
        panel.background = element_rect(fill='gray'),
        plot.margin = margin(10, 10, 2, 2, unit='mm'),
        legend.position = c(-0.025, 1.10),
        legend.background = element_rect(fill=alpha('white', 1), color='black'),
        legend.key = element_rect(color=NA, fill=alpha('black', 0)),
        legend.justification = c(0, 1),
        legend.key.size = unit(9, 'mm'),
        legend.direction = 'vertical',
        legend.text.align = 0,
        legend.margin = unit(2, 'mm')
    )

# SNR heatmap spherical model at r0 = 225 um, to compare to 
x_range = -5:0
y_range = -10:-4
r0_plot = 0.000225
r0_ranges = r0_range[seq(1, length(r0_range), 1)]
sph_snr.plot = ggplot(
    sph_snr.dt[r0 == r0_plot,], 
    aes(x=p0*1e6, y=C0, fill=snr, z=snr)
    ) +
    geom_raster(show.legend = T) +
    stat_contour(mapping=aes(color=..level..), color='black', size=0.7, breaks=c(1, 2, 5, 9)) +
    scale_color_gradient(guide = 'none') +
    scale_fill_gradientn(
        name = 'SNR', trans = 'sqrt',
        colours = rev(c('#661f19', '#8c2a22', '#e9622e', '#f7e92b', 
                        '#86c98e', '#00b2ed', '#2061ae', '#303592')),
        breaks = c(0, (1:round(sqrt(snr_max)))^2),
        limits = c(0, snr_max),        
        guide = 'colourbar'
    ) +
    scale_x_continuous(name = expression(paste(p[0], ' [nM]')), trans='log10', expand=c(0,0),
                       breaks=10^x_range, 
                       labels=sapply(x_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    scale_y_continuous(name=expression(paste(C[0], ' [mol ', m^-2, ']')), trans='log10', expand=c(0,0), 
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
        legend.position = 'none',
        # legend.position = c(1, 0.5),
        # legend.justification = c(0, 0.5),
        # legend.background = element_rect(fill=alpha('black', 0)),
        # legend.text = element_text(size=plot_text_size, color='black'),
        # legend.title = element_text(size=plot_text_size, color='black'),
        plot.margin = margin(5, 10, 2, 0, unit='mm'),
        panel.border = element_rect(colour = 'black', fill=alpha('blue',0))
    )  
sph_snr.plot = direct.label(sph_snr.plot, method=list('bottom.pieces', cex=2, colour='black'))

# SNR spherical model, horizonstal slices
C0_plot = sph_snr.dt[, sort(unique(C0))]
C0_plot = C0_plot[seq(1, length(C0_plot), 9)]
y_range = -3:1
sph_snr_hor.plot = ggplot(sph_snr.dt[C0 %in% C0_plot & r0 == r0_plot], 
                          aes(x=p0*1e6, y=snr, group=factor(C0), color=factor(C0))) +
    geom_line(size=1.5) +
    scale_x_continuous(
        name = expression(paste(p[0], ' [nM]', sep='')), trans='log10',
        breaks = 10^x_range, expand = c(0, 0),
        labels = sapply(x_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    scale_y_continuous(
        name = 'SNR',
        trans = 'log10', 
        breaks = 10^y_range,
        labels = sapply(y_range, function (i) as.expression(bquote(10^ .(i)))),
        expand = c(0, 0)
    ) +
    coord_cartesian(ylim=c(1e-2, 3e2)) +
    scale_color_manual(
        name = expression(paste(C[0], ' [mol ', m^-2, ']')),
        values = rainbow(length(C0_plot), end=0.8),
        labels = sapply(round(log10(C0_plot), 1), function (i) as.expression(bquote(10^ .(i)))),
        breaks = C0_plot
    ) +
    annotation_logticks(sides='bl') +
    guides(color=guide_legend(ncol=3, byrow=F, override.aes = list(size=10))) +
    theme(
        text = element_text(size=plot_text_size, color='black'),
        axis.text = element_text(size=plot_text_size, color='black'),
        axis.title = element_text(size=plot_text_size, color='black'),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color='gray89', linetype='dashed'),
        panel.border = element_rect(colour = 'black', fill=alpha('blue',0)),
        panel.background = element_rect(fill='gray'),
        plot.margin = margin(10, 10, 2, 2, unit='mm'),
        legend.position = c(-0.025, 1.10),
        legend.background = element_rect(fill=alpha('white', 1), color='black'),
        legend.key = element_rect(color=NA, fill=alpha('black', 0)),
        legend.justification = c(0, 1),
        legend.key.size = unit(9, 'mm'),
        legend.direction = 'vertical',
        legend.text.align = 0,
        legend.margin = unit(2, 'mm')
    )


# This generates figure 4
pdf('../figures/fig4_analytical_pde_uniform_snr.pdf', useDingbats=F, width=14, height=11.7)
multiplot(snr.plot, snr_hor.plot, sph_snr.plot, sph_snr_hor.plot,
          cols=2, annotate=T, vjust=0, hjust=-0.20, t=-5)
dev.off()