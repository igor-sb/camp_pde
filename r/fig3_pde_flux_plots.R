#
#   Fixed PDE flux model. How does PDE change gradient, mean concentration, signal, noise.
#

# Clear all variables
rm(list=ls(all=T))

# Load libraries
library(ggplot2)
library(data.table)
library(directlabels)

# Set working directory
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# Include multiplot function for generating panels of ggplots
source('include_multiplot.R')

# Plot parameters
plot_text_size = 25
legend_text_size = 19
snr_max = 20

# Load Model parameters
params.data <- read.table('../data/model_parameters.txt', header=T, sep='\t', row.names=1)
Km = params.data['Km',]$value # 10 uM in mol/m^3
Kd = params.data['Kd',]$value # (in units mol/m^3) 30 nM Kd for cAMP
Rt = params.data['Rt',]$value # total number of cAR1 receptors per cell
rc = params.data['rc',]$value # cell radius (m)
Is = params.data['Is',]$value # sampling fold (1.8 from van Haastert, Biophys J., 93, 2007
Sb = params.data['Sb',]$value # non-receptor noise, (25 from van Haastert ^)
w  = params.data['w', ]$value # width of our microfluidic device in meters (1mm))
Dc = params.data['Dc',]$value # [m^2/s] cAMP diffusion
k2 = params.data['k2',]$value # s^-1


# Define functions for SNR, CI calculations
p = function (c) c/(c+Kd)
ci_f = function (snr) snr/(snr + 1)

# Load data from MATLAB/COMSOL
input_files = c('../data/matlab_output_pde_flux_k2_13300_s-1_Km_10_uM_fine_1.txt',
                '../data/matlab_output_pde_flux_k2_13300_s-1_Km_10_uM_fine_2.txt')

snr.dt = rbindlist(lapply(input_files, fread, sep='\t', select=c(1, 2, 5, 6, 7)))
snr.dt = unique(snr.dt)
snr.dt = snr.dt[order(p0, c0_left)]

# Calculate neccessary quantities
snr.dt[, c_front := c_mean + dcdx_mean*rc/2]
snr.dt[, c_back := c_mean - dcdx_mean*rc/2]
snr.dt[, delta_r := 0.5*Rt*(p(c_front) - p(c_back))]
snr.dt[, sigma_delta_r := sqrt(Sb^2 + 0.5*Rt*p(c_front)*(1-p(c_front))/Is + 0.5*Rt*p(c_back)*(1-p(c_back))/Is)]
snr.dt[, snr := delta_r/sigma_delta_r]
snr.dt[, c0_over_kd := 0.5*c0_left/Kd]

p0_min = snr.dt[, min(p0)]
snr.dt[, snr_rel := snr/.SD[p0==p0_min, snr]]
snr.dt[, c0_rel := c_mean/.SD[p0==p0_min, c_mean]]
snr.dt[, dcdx_rel := dcdx_mean/.SD[p0==p0_min, dcdx_mean]]
snr.dt[, delta_r_over_delta_r0 := delta_r/.SD[p0==p0_min, delta_r]]

#
# Generate plots
#

# Plot the gradient , relative gradient and mean concentration as a function of PDE concentration
x_range = -14:-9
y_range = -3:2
c_all.plot = ggplot(snr.dt, aes(x=p0)) +
    geom_line(mapping = aes(y=c0_rel, linetype='dashed'), size = 1.5) +
    geom_line(mapping = aes(y=dcdx_rel, linetype='dotted'), size=1.5) +
    geom_line(mapping = aes(y=dcdx_rel/c0_rel, linetype='solid'), size=1.5) +
    scale_x_continuous(
        name = expression(paste(P[0], ' [mol ', m^-2, s^{-1}, ']', sep='')),
        trans = 'log10', expand = c(0, 0),
        breaks = 10^x_range,
        labels = sapply(x_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    scale_y_continuous(
        name = '', trans = 'log10',
        breaks = 10^y_range,
        labels = sapply(y_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    scale_linetype_manual(
        values = c('dashed', 'solid', 'dotted'),
        labels = c( 
            expression(over( paste(c[0], ' '),
                             c[paste('0, ', p[0], '=0')])),
            expression(over( paste('|', nabla, 'c|'), 
                             paste('|', nabla, 'c|'[paste(p[0], '=0')]))), 
            expression(over( paste('|', nabla, 'c|/', c[0]),
                             paste('(|', nabla, 'c|/', c[0], ')'[paste(p[0], '=0')])))
        )
    ) +
    coord_cartesian(ylim=c(1e-3, 1e2)) +
    guides(linetype=guide_legend(ncol=1, byrow=T, override.aes = list(size=2))) +
    annotation_logticks(sides='bl') +
    theme(
        text = element_text(size=plot_text_size, color='black'),
        axis.text = element_text(size=plot_text_size, color='black'),
        axis.title = element_text(size=plot_text_size, color='black'),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype='dashed', color='gray'),
        panel.border = element_rect(colour = 'black', fill=alpha('blue',0)),
        panel.background = element_rect(fill='white'),
        plot.margin = margin(2, 10, 2, 5, unit='mm'),
        plot.background = element_blank(), 
        legend.position = c(0, 0.5),
        # legend.background = element_rect(fill=alpha('white', 1)),
        # legend.background = element_blank(), 
        legend.key = element_rect(color=NA, fill=alpha('black', 0)),
        legend.justification = c(0, 0.5),
        legend.key.size = unit(11, 'mm'),
        legend.key.width = unit(20, 'mm'),
        legend.key.height = unit(25, 'mm'),
        legend.direction = 'vertical',
        legend.text.align = 0,
        legend.text = element_text(size=1.1*plot_text_size),
        legend.title = element_blank()
    )      

# c_all.plot


# Absolute receptor difference (signal)

cl_plot = snr.dt[, sort(unique(c0_left/Kd))]
cl_plot = cl_plot[seq(1, 36, 5)]
cl_kd_plot = cl_plot*Kd

c0_plot = snr.dt[, sort(unique(c0_over_kd))]
c0_plot = c0_plot[seq(3, length(c0_plot), 5)]
x_range = -14:-9
y_range = -1:3

rec_dif.plot = ggplot(snr.dt[c0_left %in% cl_kd_plot], 
        aes(x=p0, y=delta_r, group=factor(c0_left/Kd), color=factor(c0_left/Kd))) +
    geom_line(size=1.5) +
    scale_x_continuous(name = expression(paste(P[0], ' [mol ', m^-2, s^{-1}, ']', sep='')), trans='log10',
                       breaks = 10^x_range, expand = c(0, 0),
                       labels = sapply(x_range, function (i) as.expression(bquote(10^ .(i))))) +
    scale_y_continuous(name = expression(paste('signal ', Delta, R)), trans = 'log10',
                       breaks = 10^y_range, expand = c(0, 0),
                       labels = sapply(y_range, function (i) as.expression(bquote(10^ .(i))))
                       
    ) +
    coord_cartesian(ylim=c(1e-1, 3e3)) +
    scale_color_manual(
        name = expression(paste('c(-w/2) [', K[d], ']', sep='')),
        values = rainbow(length(cl_plot), end=0.8),
        labels = sapply(log10(cl_plot), function (i) as.expression(bquote(10^ .(i)))),
        breaks = cl_plot
    ) +
    annotation_logticks(sides='bl') +
    guides(color=guide_legend(ncol=2, byrow=T, override.aes = list(size=10))) +
    theme(
        text = element_text(size=plot_text_size, color='black'),
        axis.text = element_text(size=plot_text_size, color='black'),
        axis.title = element_text(size=plot_text_size, color='black'),
        #axis.ticks = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color='gray89', linetype='dashed'),
        panel.border = element_rect(colour = 'black', fill=alpha('blue',0)),
        panel.background = element_rect(fill='gray'),
        plot.margin = margin(2, 10, 2, 2, unit='mm'),
        legend.position = 'none'
        # legend.position = c(0, 1),
        # legend.background = element_rect(fill=alpha('white', 1), color='black'),
        # legend.key = element_rect(color=NA, fill=alpha('black', 0)),
        # legend.justification = c(0, 1),
        # legend.key.size = unit(9, 'mm'),
        # legend.direction = 'vertical',
        # legend.text.align = 0,
        # legend.margin = unit(2, 'mm')
    )
# rec_dif.plot

#
# Plot relative signal
#
rel_sig.plot = ggplot(snr.dt[c0_left %in% cl_kd_plot], 
                      aes(x=p0, y=delta_r_over_delta_r0, group=factor(c0_left/Kd), color=factor(c0_left/Kd))) +
    geom_line(size=1.5) +
    scale_x_continuous(name = expression(paste(P[0], ' [mol ', m^-2, s^{-1}, ']', sep='')), trans='log10',
                       breaks = 10^x_range, expand = c(0, 0),
                       labels = sapply(x_range, function (i) as.expression(bquote(10^ .(i))))) +
    scale_y_continuous(name = expression(paste('relative signal ', Delta, R)), trans = 'log10',
                       breaks = 10^y_range, expand = c(0, 0),
                       labels = sapply(y_range, function (i) as.expression(bquote(10^ .(i))))
                       
    ) +
    coord_cartesian(ylim=c(1e-1, 3e3)) +
    scale_color_manual(
        name = expression(paste('c(-w/2) [', K[d], ']', sep='')),
        values = rainbow(length(cl_plot), end=0.8),
        labels = sapply(round(log10(cl_plot), 0), function (i) as.expression(bquote(10^ .(i)))),
        breaks = cl_plot
    ) +
    annotation_logticks(sides='bl') +
    guides(color=guide_legend(ncol=2, byrow=T, override.aes = list(size=10))) +
    theme(
        text = element_text(size=plot_text_size, color='black'),
        axis.text = element_text(size=plot_text_size, color='black'),
        axis.title = element_text(size=plot_text_size, color='black'),
        #axis.ticks = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color='gray89', linetype='dashed'),
        panel.border = element_rect(colour = 'black', fill=alpha('blue',0)),
        panel.background = element_rect(fill='gray'),
        plot.margin = margin(2, 10, 2, 2, unit='mm'),
        legend.position = 'none'
        # legend.position = c(0, 1),
        # legend.background = element_rect(fill=alpha('white', 1), color='black'),
        # legend.key = element_rect(color=NA, fill=alpha('black', 0)),
        # legend.justification = c(0, 1),
        # legend.key.size = unit(9, 'mm'),
        # legend.direction = 'vertical',
        # legend.text.align = 0,
        # legend.margin = unit(2, 'mm')
    )
rel_sig.plot

# Plot absolute noise

sigma.plot = ggplot(snr.dt[c0_left %in% cl_kd_plot], aes(x=p0, y=sigma_delta_r, group=factor(c0_left/Kd), color=factor(c0_left/Kd))) +
    geom_line(size=1.5) +
    scale_x_continuous(name = expression(paste(P[0], ' [mol ', m^-2, s^{-1}, ']', sep='')), trans='log10',
                       breaks = 10^x_range, expand = c(0, 0),
                       labels = sapply(x_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    scale_y_continuous(
        name = expression(paste('noise ', sqrt(sigma[B]^2 + sigma[paste(Delta, 'R')]^2 / I))),
        breaks = seq(0, floor(sqrt(Sb^2 + 0.25*Rt/Is)), 25)
        # expand = c(0, 0)
    ) +
    scale_color_manual(
        name = expression(atop('c(-w/2)', paste(' [', K[d], ']', sep=''))),
        values = rainbow(length(cl_plot), end=0.8),
        labels = sapply(round(log10(cl_plot), 0), function (i) as.expression(bquote(10^ .(i)))),
        breaks = cl_plot
    ) +
    # coord_cartesian(ylim=c(0, 100)) +
    annotation_logticks(sides='b') +
    guides(color=guide_legend(ncol=1, byrow=F, override.aes = list(size=10))) +
    theme(
        text = element_text(size=plot_text_size, color='black'),
        axis.text = element_text(size=plot_text_size, color='black'),
        axis.title = element_text(size=plot_text_size, color='black'),
        #axis.ticks = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color='gray89', linetype='dashed'),
        panel.border = element_rect(colour = 'black', fill=alpha('blue',0)),
        panel.background = element_rect(fill='gray'),
        plot.margin = margin(2, 10, 2, 2, unit='mm'), # t r b l
        # legend.position = 'none'
        legend.title = element_text(size=plot_text_size*0.8, color='black'),
        legend.position = c(1, 1),
        legend.background = element_rect(fill=alpha('white', 1), color='black'),
        legend.key = element_rect(color=NA, fill=alpha('black', 0)),
        legend.justification = c(1, 1),
        legend.key.size = unit(9, 'mm'),
        legend.direction = 'vertical',
        legend.text.align = 0,
        legend.margin = unit(0, 'mm')
    )
# sigma.plot




# SNR heatmap

x_range = -14:-9
y_range = -2:5
snr.plot = ggplot(snr.dt, aes(y = c0_left/Kd, x = p0, z = snr, fill = snr)) +
    geom_raster(show.legend = T) +
    stat_contour(mapping=aes(color=..level..), color='black', size=0.7, breaks=c(1, 2, 5)) +
    scale_y_continuous(
        name = expression(paste('c(-w/2) [', K[d], ']', sep='')), trans='log10', expand = c(0, 0),
        breaks = 10^y_range,
        labels = sapply(y_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    scale_x_continuous(
        name = expression(paste(P[0], ' [mol ', m^-2, s^{-1}, ']', sep='')), trans='log10', expand = c(0, 0),
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
        plot.margin = margin(5, 10, 2, 2, unit='mm'),
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
# snr.plot


# snr  horizontal slices
x_range = -14:-9
y_range = -3:1

snr_hor.plot = ggplot(snr.dt[c0_left %in% cl_kd_plot], aes(x=p0, y=snr, group=factor(c0_left/Kd), color=factor(c0_left/Kd))) +
    geom_line(size=1.5) +
    scale_x_continuous(
        name = expression(paste(P[0], ' [mol ', m^-2, s^{-1}, ']', sep='')), trans='log10',
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
        values = rainbow(length(cl_plot), end=0.8),
        labels = sapply(round(log10(cl_plot), 0), function (i) as.expression(bquote(10^ .(i)))),
        breaks = cl_plot
    ) +
    annotation_logticks(sides='bl') +
    guides(color=guide_legend(ncol=5, byrow=F, override.aes = list(size=10))) +
    theme(
        text = element_text(size=plot_text_size, color='black'),
        axis.text = element_text(size=plot_text_size, color='black'),
        axis.title = element_text(size=plot_text_size, color='black'),
        #axis.ticks = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color='gray89', linetype='dashed'),
        panel.border = element_rect(colour = 'black', fill=alpha('blue',0)),
        panel.background = element_rect(fill='gray'),
        plot.margin = margin(10, 10, 2, 2, unit='mm'),
        # legend.position = 'none'
        legend.position = c(-0.025, 1.10),
        legend.background = element_rect(fill=alpha('white', 1), color='black'),
        legend.key = element_rect(color=NA, fill=alpha('black', 0)),
        legend.justification = c(0, 1),
        legend.key.size = unit(9, 'mm'),
        legend.direction = 'vertical',
        legend.text.align = 0,
        legend.margin = unit(2, 'mm')
    )
# snr_hor.plot


# SNR ratio heatmap
x_range = -14:-9
y_range = -2:5
snr_ratio_range = -4:4

snr_ratio.plot = ggplot(snr.dt[snr_rel > 1e-4], aes(y=c0_left/Kd, x=p0, z=snr_rel, fill=snr_rel)) +
    geom_raster(show.legend = T) +
    stat_contour(mapping=aes(color=..level..), color='black', size=0.7, breaks=c(2, 10, 100)) +
    scale_y_continuous(
        name = expression(paste('c(-w/2) [', K[d], ']', sep='')), trans='log10', expand = c(0, 0),
        breaks = 10^y_range,
        labels = sapply(y_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    scale_x_continuous(
        name = expression(paste(P[0], ' [mol ', m^-2, s^{-1}, ']', sep='')), trans='log10', expand = c(0, 0),
        breaks = 10^x_range,
        labels = sapply(x_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    scale_fill_gradientn(
        name = expression(over('SNR', 'SNR'[paste('p'[0], '=0')])), 
        trans = 'log10',
        colours = rev(c('#661f19', '#8c2a22', '#e9622e', '#f7e92b', 
                        '#86c98e', '#00b2ed', '#2061ae', '#303592')),
        breaks = 10^snr_ratio_range,
        labels = sapply(snr_ratio_range, function (i) as.expression(bquote(10^ .(i)))),
        guide = 'colourbar'
    ) +
    scale_color_gradient(guide = 'none') +
    annotation_logticks(sides='bl', color='black') +
    theme(
        text = element_text(size=plot_text_size, color='black'),
        axis.text = element_text(size=plot_text_size, color='black'),
        axis.title = element_text(size=plot_text_size, color='black'),
        axis.ticks = element_blank(),
        plot.margin = margin(5, 10, 2, 2, unit='mm'),
        panel.background = element_rect(fill='#303592'),
        legend.position = c(0, 1),
        legend.text = element_text(size=legend_text_size, color='black'),
        legend.title = element_text(size=legend_text_size, color='black'),
        legend.background = element_rect(fill=alpha('black', 0)),
        legend.justification = c(0, 1),
        legend.key.size = unit(13, 'mm'),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = 'black', fill=alpha('blue',0))        
    )
snr_ratio.plot = direct.label(snr_ratio.plot, method=list('bottom.pieces', cex=2, colour='black'))
#snr_ratio.plot

# snr ratio horizontal slices
y_range = -3:4
x_range = -14:-9
#c0_plot = c0_range[1:61][seq(1, 61, 5)]

snr_ratio_hor.plot = ggplot(snr.dt[c0_over_kd %in% c0_plot], aes(x=p0, y=snr_rel, group=factor(c0_over_kd), color=factor(c0_over_kd))) +
    geom_line(size=1.5) +
    scale_x_continuous(
        name = expression(paste(P[0], ' [mol ', m^-2, s^{-1}, ']', sep='')), trans='log10',
        breaks = 10^x_range, expand = c(0, 0),
        labels = sapply(x_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    scale_y_continuous(
        name = expression(paste('SNR/', 'SNR'[paste(p[0], '=0')])), 
        trans = 'log10', 
        breaks = 10^y_range,
        labels = sapply(y_range, function (i) as.expression(bquote(10^ .(i)))),
        expand = c(0, 0)
    ) +
    coord_cartesian(ylim=c(1e-1, 1e5)) +
    scale_color_manual(
        name = expression(paste(c[0], ' [', K[d], ']', sep='')),
        values = rainbow(length(c0_plot), end=0.8),
        labels = sapply(round(log10(c0_plot), 0), function (i) as.expression(bquote(10^ .(i)))),
        breaks = c0_plot
    ) +
    annotation_logticks(sides='bl') +
    guides(color=guide_legend(ncol=3, byrow=F, override.aes = list(size=10))) +
    theme(
        text = element_text(size=plot_text_size, color='black'),
        axis.text = element_text(size=plot_text_size, color='black'),
        axis.title = element_text(size=plot_text_size, color='black'),
        #axis.ticks = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color='gray89', linetype='dashed'),
        panel.border = element_rect(colour = 'black', fill=alpha('blue',0)),
        panel.background = element_rect(fill='gray'),
        plot.margin = margin(10, 10, 2, 2, unit='mm'),
        legend.position = 'none'
        # legend.position = c(0, 1),
        # legend.background = element_rect(fill=alpha('white', 1), color='black'),
        # legend.key = element_rect(color=NA, fill=alpha('black', 0)),
        # legend.justification = c(0, 1),
        # legend.key.size = unit(9, 'mm'),
        # legend.direction = 'vertical',
        # legend.text.align = 0,
        # legend.margin = unit(2, 'mm')
    )
#snr_ratio_hor.plot

#
# This plots the entire figure
#
pdf('../figures/fig3_pde_flux_plots.pdf', useDingbats=F, width=14, height=11.7)
multiplot(c_all.plot, rec_dif.plot, sigma.plot, snr_ratio.plot, 
          cols=2, annotate=T, vjust=0, hjust=-0.20, t=-2)
dev.off()
# Export 14 by 17.5 inch for 2x3 plots, 14x11.7 for 2x2 plots


