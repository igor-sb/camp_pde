#
# Plots Fixed PDE secretion flux model, where dcdx = const.
#

# Clear all variables
rm(list=ls(all=TRUE))

# Load plotting packages
library(ggplot2)
library(data.table)
library(directlabels)

# Set working directory
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# Load multiplot for ggplot
source('include_multiplot_2.R')

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


# Plotting parameters
plot_text_size = 22
legend_text_size = 18
heat_legend_size = 10
legend_key_size = unit(2, 'pt')
legend_key_ht = unit(1, 'pt')
slice_legend_pos = c(-0.09, 1.30)
# camp_left = expression(paste('c(-w)'))

# Load the data from COMSOL/MATLAB
input_files = c('../data/matlab_output_pde_flux_k2_13300_s-1_Km_10_uM_fixed_dcdx_0.24_mol_per_m3.txt',
                '../data/matlab_output_pde_flux_k2_13300_s-1_Km_10_uM_fixed_dcdx_0.12_mol_per_m3.txt',
                '../data/matlab_output_pde_flux_k2_13300_s-1_Km_10_uM_fixed_dcdx_0.06_mol_per_m3.txt')

p = function (c) c/(c+Kd)
erf = function (x) 2 * pnorm(x * sqrt(2)) - 1
# Error function can be approximated with erf(x/sqrt(2)) \approx sqrt(2/pi)*x
# test this. 
erf_rel_error = function (x) {
    erf1 = erf(x/sqrt(2))
    erf2 = sqrt(2/pi)*x
    return(abs(erf1-erf2)/(erf1+erf2)/2)
}


camp.dt = rbindlist(lapply(input_files, function (file) {
    dt = fread(file, sep='\t', select=c(1, 2, 5, 6, 7))
    # Select only subset of points for c0_Kd
    dt[, c0_mean_over_kd := strsplit(strsplit(file, '_')[[1]][12], 'K')[[1]][1]]
    return(dt)
}))
camp.dt[, dcdx_bulk := factor(round((c0_left/w - c0_right/w), 2))]
# camp.dt[c0_mean_over_kd=='1', 10^round(log10(unique(dcdx_bulk)), 2)]
camp.dt[, c_front := c_mean + dcdx_mean*rc/2]
camp.dt[, c_back := c_mean - dcdx_mean*rc/2]
camp.dt[, delta_r := 0.5*Rt*(p(c_front) - p(c_back))]
camp.dt[, sigma_delta_r := sqrt(Sb^2 + 0.5*Rt*p(c_front)*(1-p(c_front))/Is + 0.5*Rt*p(c_back)*(1-p(c_back))/Is)]
camp.dt[, snr := delta_r/sigma_delta_r]
camp.dt[, ci := erf(snr/sqrt(2))*(1 - 0.2*erf(snr/sqrt(2)))]
camp.dt[, c0 := 10^round(log10(0.5*c0_left + 0.5*c0_right), 2)]
# camp.dt[, c0_left_over_kd := c0_left/Kd]
# There are some duplicate entries. Remove the duplicates.
camp.dt = camp.dt[, if (.N>1) unique(.SD) else .SD, by=c('c0_left', 'p0')]
camp.dt[, dcdx_mean_label := paste0('paste("|", nabla, "c|"[app], " = ",', dcdx_bulk, ', " nM ", mu, m^{-1}, "")')]

#
#          SNR HEAT MAP
#
# Make a similar plot but now with SNR values instead of the 
# difference of occupied receptors.

x_range = -14:-9
y_range = -5:-1
# Fun fact: Gradient is numerically equal in (mol/m^3)/m and (nM/um).

snr.plot = ggplot(camp.dt, aes(x = p0, y = c0, fill = snr, z = snr)) +
    facet_wrap(~ dcdx_mean_label, scales = 'free', labeller = label_parsed) +
    geom_raster(show.legend = T) + 
    geom_contour(mapping = aes(color = ..level..), color = 'black', size = 0.7, breaks = c(0.1, 0.5, 1)) +
    scale_x_continuous(
        name = expression(paste('PDE secretion flux ', P[0], ' [mol ', m^-2, s^{-1}, ']'), sep=''),
        breaks = 10^x_range, expand = c(0,0), trans = 'log10',
        labels = sapply(x_range, function (i) as.expression(bquote(10^ .(i))))
    ) + 
    scale_y_continuous(
        name = expression(paste(c['0,app'], ' [', K[d], ']'), sep=''),
        breaks = 10^y_range, expand = c(0,0), trans = 'log10',
        labels = sapply(y_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    annotation_logticks(sides = 'bl', color = 'white') +
    scale_fill_gradientn(
        name = 'SNR', trans = 'sqrt', guide = 'colourbar',
        colours = rev(c('#661f19', '#8c2a22', '#e9622e', '#f7e92b', 
                        '#86c98e', '#00b2ed', '#2061ae', '#303592')),
        limits = c(0, camp.dt[, max(snr)]),
        breaks = c(0, 0.1, 0.5, 1, 2)
    ) +
    scale_color_gradient(guide = 'none') +
    coord_cartesian(ylim=c()) +
    theme(
        text = element_text(size=plot_text_size, color='black'),
        axis.text = element_text(size=plot_text_size, color='black'),
        axis.title = element_text(size=plot_text_size, color='black'),
        axis.ticks = element_blank(),
        plot.margin = margin(10, 10, 2, 2, unit='mm'),
        legend.position = c(0, 1),
        legend.text = element_text(size=legend_text_size, color='white'),
        legend.title = element_text(size=legend_text_size, color='white'),
        legend.background = element_rect(fill=alpha('black', 0)),
        legend.justification = c(0, 1),
        legend.key.size = unit(heat_legend_size, 'mm'),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = 'black', fill=alpha('blue',0))
    )
snr.plot = direct.label(snr.plot, method=list('top.pieces', cex=2, colour='black'))
#snr.plot

# snr vertical slices

x_range = -5:1
y_range = -4:1
p0_range = camp.dt[, sort(unique(p0))][1:41][seq(1, 41, length.out = 12)]

snr_v.plot <- ggplot(camp.dt[p0 %in% p0_range], aes(x=c0, y=snr, group=factor(p0), color=factor(p0))) +
    facet_wrap(~ dcdx_mean_label, scales = 'free', labeller = label_parsed) +
    geom_line(size=1.5) + 
    scale_y_continuous(name = 'SNR', 
                       breaks = 10^y_range,
                       labels = sapply(y_range, function (i) as.expression(bquote(10^ .(i)))),
                       trans = 'log10') + 
    coord_cartesian(ylim=c(2e-4, 1e0)) +
    scale_x_continuous(
        name = expression(paste(c['0,app'], ' [', K[d], ']'), sep=''),
        breaks = 10^x_range, expand = c(0,0), trans = 'log10',
        labels = sapply(x_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    annotation_logticks(color='black', sides='bl') +
    guides(color=guide_legend(nrow=1, byrow=F, override.aes = list(size=10))) +
    scale_color_manual(
        name = expression(paste(P[0], ' [mol ', m^-2, s^-1, ']')),
        values = rainbow(length(p0_range), end=0.8),
        labels = sapply(round(log10(p0_range), 1), function (i) as.expression(bquote(10^ .(i))))
    ) +
    theme(
        strip.text = element_blank(),
        text = element_text(size=plot_text_size, color='black'),
        axis.text = element_text(size=plot_text_size, color='black'),
        axis.title = element_text(size=plot_text_size, color='black'),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color='gray89', linetype='dashed'),
        panel.border = element_rect(colour = 'black', fill=alpha('blue',0)),
        panel.background = element_rect(fill='gray'),
        plot.margin = margin(30, 10, 2, 2, unit='mm'), # t r b l
        legend.text = element_text(size=legend_text_size, color='black'),
        legend.position = slice_legend_pos,
        legend.background = element_blank(),
        legend.justification = c(0, 1),
        legend.key = element_blank(), 
        legend.key.size = legend_key_size,
        legend.key.height = legend_key_ht,
        legend.direction = 'vertical',
        legend.text.align = 0
    )
#snr_v.plot


# snr horizontal slices


c0_values = camp.dt[, sort(unique(c0))][seq(1, 42, length.out=12)]
x_range = -14:-9
y_range = -4:1

snr_h.plot = ggplot(camp.dt[c0 %in% c0_values], 
                    aes(x=p0, y=snr, group=factor(c0), color=factor(c0))) +
    facet_wrap(~ dcdx_mean_label, scales = 'free', labeller = label_parsed) +
    geom_line(size=1.5) + 
    scale_y_continuous(
        name = 'SNR', 
        breaks = 10^y_range, expand = c(0,0), trans='log10',
        labels = sapply(y_range, function (i) as.expression(bquote(10^ .(i))))
    ) + 
    scale_x_continuous(
        name = expression(paste(P[0],' [mol ', m^-2, s^-1, ']'), sep=''),
        breaks = 10^x_range, expand = c(0,0), trans = 'log10',
        labels = sapply(x_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    coord_cartesian(ylim=c(2e-4, 3e0)) +
    annotation_logticks(color='black', sides='bl') +
    scale_color_manual(
        name = expression(paste(c['0,app'], ' [', K[d], ']')),
        values = rainbow(length(c0_values), end=0.8),
        labels = sapply(round(log10(c0_values), 1), function (i) as.expression(bquote(10^ .(i)))),
        guide = 'colourbar'
    ) +
    guides(color=guide_legend(nrow=1, byrow=F, override.aes = list(size=10))) +
    theme(
        strip.text = element_blank(),
        text = element_text(size=plot_text_size, color='black'),
        axis.text = element_text(size=plot_text_size, color='black'),
        axis.title = element_text(size=plot_text_size, color='black'),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color='gray89', linetype='dashed'),
        panel.border = element_rect(colour = 'black', fill=alpha('blue',0)),
        panel.background = element_rect(fill='gray'),
        plot.margin = margin(30, 10, 2, 2, unit='mm'),
        legend.text = element_text(size=legend_text_size, color='black'),
        legend.position = slice_legend_pos,
        legend.background = element_blank(),
        legend.justification = c(0, 1),
        legend.key = element_blank(), 
        legend.key.size = legend_key_size,
        legend.key.height = legend_key_ht,
        legend.direction = 'vertical',
        # legend.title.align = -4e7,
        legend.text.align = 0
    )
# snr_h.plot


# ci horizontal slices

c0_values = camp.dt[, sort(unique(c0))][seq(1, 42, length.out=12)]
x_range = -14:-9
y_range = -4:1

ci_h.plot = ggplot(camp.dt[c0 %in% c0_values], 
                    aes(x=p0, y=ci, group=factor(c0), color=factor(c0))) +
    facet_wrap(~ dcdx_mean_label, scales = 'free', labeller = label_parsed) +
    geom_line(size=1.5) + 
    scale_y_continuous(
        name = 'CI', 
        breaks = 10^y_range, expand = c(0,0), trans='log10',
        labels = sapply(y_range, function (i) as.expression(bquote(10^ .(i))))
    ) + 
    scale_x_continuous(
        name = expression(paste(P[0],' [mol ', m^-2, s^-1, ']'), sep=''),
        breaks = 10^x_range, expand = c(0,0), trans = 'log10',
        labels = sapply(x_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    coord_cartesian(ylim=c(2e-4, 3e0)) +
    annotation_logticks(color='black', sides='bl') +
    scale_color_manual(
        name = expression(paste(c['0,app'], ' [', K[d], ']')),
        values = rainbow(length(c0_values), end=0.8),
        labels = sapply(round(log10(c0_values), 1), function (i) as.expression(bquote(10^ .(i)))),
        guide = 'colourbar'
    ) +
    guides(color=guide_legend(nrow=1, byrow=F, override.aes = list(size=10))) +
    theme(
        strip.text = element_blank(),
        text = element_text(size=plot_text_size, color='black'),
        axis.text = element_text(size=plot_text_size, color='black'),
        axis.title = element_text(size=plot_text_size, color='black'),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color='gray89', linetype='dashed'),
        panel.border = element_rect(colour = 'black', fill=alpha('blue',0)),
        panel.background = element_rect(fill='gray'),
        plot.margin = margin(30, 10, 2, 2, unit='mm'),
        legend.text = element_text(size=legend_text_size, color='black'),
        legend.position = slice_legend_pos,
        legend.background = element_blank(),
        legend.justification = c(0, 1),
        legend.key = element_blank(), 
        legend.key.size = legend_key_size,
        legend.key.height = legend_key_ht,
        legend.direction = 'vertical',
        # legend.title.align = -4e7,
        legend.text.align = 0
    )
#ci_h.plot

# snr vertical slices

x_range = -5:1
y_range = -4:1
p0_range = camp.dt[, sort(unique(p0))][1:41][seq(1, 41, length.out = 12)]

ci_v.plot <- ggplot(camp.dt[p0 %in% p0_range], aes(x=c0, y=ci, group=factor(p0), color=factor(p0))) +
    facet_wrap(~ dcdx_mean_label, scales = 'free', labeller = label_parsed) +
    geom_line(size=1.5) + 
    scale_y_continuous(name = 'CI', 
                       breaks = 10^y_range,
                       labels = sapply(y_range, function (i) as.expression(bquote(10^ .(i)))),
                       trans = 'log10') + 
    coord_cartesian(ylim=c(2e-4, 1e0)) +
    scale_x_continuous(
        name = expression(paste(c['0,app'], ' [', K[d], ']'), sep=''),
        breaks = 10^x_range, expand = c(0,0), trans = 'log10',
        labels = sapply(x_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    annotation_logticks(color='black', sides='bl') +
    guides(color=guide_legend(nrow=1, byrow=F, override.aes = list(size=10))) +
    scale_color_manual(
        name = expression(paste(P[0], ' [mol ', m^-2, s^-1, ']')),
        values = rainbow(length(p0_range), end=0.8),
        labels = sapply(round(log10(p0_range), 1), function (i) as.expression(bquote(10^ .(i))))
    ) +
    theme(
        strip.text = element_blank(),
        text = element_text(size=plot_text_size, color='black'),
        axis.text = element_text(size=plot_text_size, color='black'),
        axis.title = element_text(size=plot_text_size, color='black'),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color='gray89', linetype='dashed'),
        panel.border = element_rect(colour = 'black', fill=alpha('blue',0)),
        panel.background = element_rect(fill='gray'),
        plot.margin = margin(30, 10, 2, 2, unit='mm'), # t r b l
        legend.text = element_text(size=legend_text_size, color='black'),
        legend.position = slice_legend_pos,
        legend.background = element_blank(),
        legend.justification = c(0, 1),
        legend.key = element_blank(), 
        legend.key.size = legend_key_size,
        legend.key.height = legend_key_ht,
        legend.direction = 'vertical',
        legend.text.align = 0
    )
# ci_v.plot


# Generate entire figure for constant c0_mean
pdf('../figures/si_pde_flux_fixed_dcdx_snr_plots.pdf', useDingbats = F, width=14, height=23.3)
multiplot(snr.plot, snr_h.plot, snr_v.plot, ci_v.plot,
          cols=1, annotate=T, vjust=0, t=-35, t1=0, hjust=-0.07, new_lines=3)
dev.off()