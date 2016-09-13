#
# Plots for the difference in the number of occupied
# receptors vs PDE secretion rate for the COMSOL
# model with the cell in the middle of the microfluidic
# device.
#
# We vary:
#   p0 - PDE secretion rate 
#   c0_L - cAMP concentration on the left side
#
# Author: Igor Segota
# Date:   12/07/2013

# To do:
# ------
# Make "vertical" line graphs (CI vs bulk mean conc.) for various PDE
# secretion rates to show how the dynamic range of cAMP detection 
# is improved when you add PDE.

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
source('include_multiplot.R')


# Load Model parameters
params.data <- read.table('../data/model_parameters.txt', header=T, sep='\t', row.names=1)
Km = params.data['Km',]$value # 10 uM in mol/m^3, the data file says 1uM but it looks like a typo.
Kd = params.data['Kd',]$value # (in units mol/m^3) 30 nM Kd for cAMP
Rt = params.data['Rt',]$value # total number of cAR1 receptors per cell
rc = params.data['rc',]$value # cell radius (m)
Is = params.data['Is',]$value # sampling fold (1.8 from van Haastert, Biophys J., 93, 2007
Sb = params.data['Sb',]$value # non-receptor noise, (25 from van Haastert ^)
w  = params.data['w', ]$value # width of our microfluidic device in meters (1mm wide in x dir)
Dc = params.data['Dc',]$value # [m^2/s] cAMP diffusion
k2 = params.data['k2',]$value # s^-1

# Plotting parameters
plot_text_size = 25
legend_text_size = 20
heat_legend_size = 10
legend_key_size = unit(2, 'pt')
legend_key_ht = unit(1, 'pt')
camp_left = expression(paste('c(-w)'))

# Load data from MATLAB/COMSOL
input_files = c('../data/matlab_output_pde_flux_k2_13300_s-1_Km_10_uM_fine_1.txt',
                '../data/matlab_output_pde_flux_k2_13300_s-1_Km_10_uM_fine_2.txt')

p = function (c) c/(c+Kd)
erf = function (x) 2 * pnorm(x * sqrt(2)) - 1
ci_f = function (snr) snr/(snr + 1)

camp.dt = rbindlist(lapply(input_files, fread, sep='\t', select=c(1, 2, 5, 6, 7)))
camp.dt[, c_front := c_mean + dcdx_mean*rc/2]
camp.dt[, c_back := c_mean - dcdx_mean*rc/2]
camp.dt[, delta_r := 0.5*Rt*(p(c_front) - p(c_back))]
camp.dt[, sigma_delta_r := sqrt(Sb^2 + 0.5*Rt*p(c_front)*(1-p(c_front))/Is + 0.5*Rt*p(c_back)*(1-p(c_back))/Is)]
camp.dt[, snr := delta_r/sigma_delta_r]
#camp.dt[, ci := erf(snr/sqrt(2))*(1 - 0.2*erf(snr/sqrt(2)))]
camp.dt[, ci := ci_f(snr)]
camp.dt[, c0_left_over_kd := c0_left/Kd]
# There are some duplicate entries. Remove the duplicates.
camp.dt = camp.dt[, if (.N>1) unique(.SD) else .SD, by=c("c0_left", "p0")]


#
#          SNR HEAT MAP
#
# Make a similar plot but now with SNR values instead of the 
# difference of occupied receptors.

x_range = -14:-9
y_range = -2:5
# plot_text_size = 22
# Fun fact: Gradient is numerically equal in (mol/m^3)/m and (nM/um).

snr.plot = ggplot(camp.dt, aes(x = p0, y = c0_left_over_kd, fill = snr, z = snr)) +
    geom_raster(show.legend = T, interpolate=T) + 
    geom_contour(mapping = aes(color = ..level..), color = "black", size = 0.7, breaks = c(1,2,5,10)) +
    scale_x_continuous(
        name = expression(paste("PDE secretion flux ", P[0], " [mol ", m^-2, s^{-1}, "]"), sep=""),
        breaks = 10^x_range, expand = c(0,0), trans = "log10",
        labels = sapply(x_range, function (i) as.expression(bquote(10^ .(i))))
    ) + 
    scale_y_continuous(
        name = expression(paste("cAMP conc. c(-w/2) [", K[d], "]"), sep=""),
        breaks = 10^y_range, expand = c(0,0), trans = "log10",
        labels = sapply(y_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    annotation_logticks(sides = "bl", color = "white") +
    scale_fill_gradientn(
        name = "SNR", trans = "sqrt", guide = 'colourbar',
        colours = rev(c("#661f19", "#8c2a22", "#e9622e", "#f7e92b", 
                        "#86c98e", "#00b2ed", "#2061ae", "#303592")),
        breaks = c(0, 1, 4, 9, 16, 25), limits=c(0, camp.dt[, max(snr)])
    ) +
    scale_color_gradient(guide = 'none') +
    theme(
        text = element_text(size=plot_text_size, color="black"),
        axis.text = element_text(size=plot_text_size, color="black"),
        axis.title = element_text(size=plot_text_size, color="black"),
        axis.ticks = element_blank(),
        plot.margin = margin(10, 10, 2, 2, unit="mm"),
        legend.position = c(1, 0.1),
        legend.text = element_text(size=legend_text_size, color="white"),
        legend.title = element_text(size=legend_text_size, color="white"),
        legend.background = element_rect(fill=alpha("black", 0)),
        legend.justification = c(1, 0),
        legend.key.size = unit(heat_legend_size, "mm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=alpha("blue",0))
    )
snr.plot = direct.label(snr.plot, method=list("bottom.pieces", cex=2, colour="black"))
# snr.plot


# snr vertical slices

x_range = -2:4
y_range = -4:1
p0_range = camp.dt[, sort(unique(p0))][1:41][seq(1, 41, length.out = 12)]

snr_v.plot <- ggplot(camp.dt[p0 %in% p0_range], aes(x=c0_left_over_kd, y=snr, group=factor(p0), color=factor(p0))) +
    geom_line(size=1.5) + 
    scale_y_continuous(name = "SNR", 
                       breaks = 10^y_range,
                       labels = sapply(y_range, function (i) as.expression(bquote(10^ .(i)))),
                       trans = "log10") + 
    coord_cartesian(ylim=c(1e-4, 3e1)) +
    scale_x_continuous(
        expression(paste("c(-w/2) [", K[d], "]")),
        breaks = 10^x_range, expand = c(0,0), trans = "log10",
        labels = sapply(x_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    annotation_logticks(color="black", sides="b") +
    guides(color=guide_legend(ncol=1, byrow=F, override.aes = list(size=10))) +
    scale_color_manual(
        name = expression(paste(P[0], " [mol ", m^-2, s^-1, "]")),
        values = rainbow(length(p0_range), end=0.8),
        labels = sapply(round(log10(p0_range), 1), function (i) as.expression(bquote(10^ .(i))))
    ) +
    theme(
        text = element_text(size=plot_text_size, color="black"),
        axis.text = element_text(size=plot_text_size, color="black"),
        axis.title = element_text(size=plot_text_size, color="black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color="gray89", linetype="dashed"),
        panel.border = element_rect(colour = "black", fill=alpha("blue",0)),
        panel.background = element_rect(fill="gray"),
        plot.margin = margin(10, 40, 2, 4, unit="mm"), # t r b l
        legend.text = element_text(size=legend_text_size, color='black'),
        legend.position = c(1, 1.14),
        legend.background = element_blank(),
        legend.justification = c(0, 1),
        legend.key = element_blank(), 
        legend.key.size = legend_key_size,
        legend.key.height = legend_key_ht,
        legend.direction = "vertical",
        legend.title.align = 100,
        legend.text.align = 0
    )
# snr_v.plot


# snr horizontal slices


c0_values = camp.dt[, sort(unique(c0_left_over_kd))][seq(1, 36, length.out=12)]
x_range = -14:-10
y_range = -4:1

snr_h.plot = ggplot(camp.dt[c0_left_over_kd %in% c0_values], 
                    aes(x=p0, y=snr, group=factor(c0_left_over_kd), color=factor(c0_left_over_kd))) +
    geom_line(size=1.5) + 
    scale_y_continuous(
        name = "SNR", 
        breaks = 10^y_range, expand = c(0,0), trans="log10",
        labels = sapply(y_range, function (i) as.expression(bquote(10^ .(i))))
    ) + 
    scale_x_continuous(
        name = expression(paste(P[0]," [mol ", m^-2, s^-1, "]"), sep=""),
        breaks = 10^x_range, expand = c(0,0), trans = "log10",
        labels = sapply(x_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    coord_cartesian(ylim=c(1e-4, 3e1)) +
    annotation_logticks(color="black", sides="bl") +
    scale_color_manual(
        name = expression(paste('c(-w/2) [', K[d], ']')),
        values = rainbow(length(c0_values), end=0.8),
        labels = sapply(round(log10(c0_values), 1), function (i) as.expression(bquote(10^ .(i)))),
        guide = 'colourbar'
    ) +
    guides(color=guide_legend(ncol=1, byrow=F, override.aes = list(size=10))) +
    theme(
        text = element_text(size=plot_text_size, color="black"),
        axis.text = element_text(size=plot_text_size, color="black"),
        axis.title = element_text(size=plot_text_size, color="black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color="gray89", linetype="dashed"),
        panel.border = element_rect(colour = "black", fill=alpha("blue",0)),
        panel.background = element_rect(fill="gray"),
        plot.margin = margin(10, 40, 2, 2, unit="mm"),
        legend.text = element_text(size=legend_text_size, color='black'),
        legend.position = c(1, 1.13),
        legend.background = element_blank(),
        legend.justification = c(0, 1),
        legend.key = element_blank(), 
        legend.key.size = legend_key_size,
        legend.key.height = legend_key_ht,
        legend.direction = "vertical",
        #legend.title = element_text(hjust=1, margin=margin(r=200, unit='mm')),
        legend.title.align = -4e7,
        legend.text.align = 0
    )
#snr_h.plot


# FOr vertical slices of the heat map plot, for each PDE secretion rate P0,
# we find the max of SNR wrt to c0_left_over_Kd (that's in opt2.dt) and calculate
# the width around this peak where SNR > 1.

p0s = camp.dt[p0 < 4e-11 , unique(p0)]
fwhms = rep(0, length(p0s))
for (i in seq(1, length(p0s))) {
    snr_max = camp.dt[p0==p0s[i], max(snr)]
    c0_max = camp.dt[p0==p0s[i] & snr==snr_max, c0_left_over_kd]
    half_max_snr_dif = camp.dt[p0==p0s[i], min(abs(snr-1))]
    c0_half_max = camp.dt[p0==p0s[i] & abs(snr-1)==half_max_snr_dif, c0_left_over_kd]
    fwhms[i] = 10^(2*abs(log10(c0_max) - log10(c0_half_max)))
}

width.dt = data.table(p0=p0s, width=fwhms)


# Plot parameters
x_range = -14:-10
y_range = 0:5

fwhm.plot = ggplot(width.dt, aes(x=p0, y=width)) +
    geom_smooth(se=F, color="black", size=1.5, span=2) + theme_bw() +
    scale_x_continuous(
        name = expression(paste(P[0], " [mol ", m^-2, s^-1, "] "), sep=""),
        breaks = 10^x_range, trans = "log10",
        labels = sapply(x_range, function (i) as.expression(bquote(10^ .(i))))
    ) + 
    scale_y_continuous(
        name = expression(paste("c(-w/2) range [", K[d], "]", sep="")),
        breaks = 10^y_range, trans = "log10", 
        labels = sapply(y_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    annotation_logticks(sides="bl") +
    # coord_cartesian(ylim=c(10, 2e5)) +
    theme(
        text = element_text(size=plot_text_size, color="black"),
        axis.text = element_text(size=plot_text_size, color="black"),
        axis.title = element_text(size=plot_text_size, color="black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype="dashed", color="gray", size=0.5),
        panel.border = element_rect(colour = "black", fill=alpha("blue",0)),
        panel.background = element_rect(fill="white"),
        plot.margin = margin(10, 40, 1, 2, unit="mm"), # t r b l
        plot.background = element_blank(), 
        legend.position = c(1, 1),
        legend.background = element_blank(), 
        legend.key = element_rect(color=NA, fill=alpha("black", 0)),
        legend.justification = c(1, 1),
        legend.key.size = unit(11, "mm"),
        legend.direction = "vertical",
        legend.text.align = 1,
        legend.text = element_text(size=legend_text_size),
        legend.title = element_blank()        
    )
#fwhm.plot


# ci heatmap

x_range = -14:-9
y_range = -2:5
# Fun fact: Gradient is numerically equal in (mol/m^3)/m and (nM/um).

ci_heat.plot = ggplot(camp.dt, aes(x = p0, y = c0_left_over_kd, fill = ci)) +
    geom_raster(interpolate = T) + 
    scale_x_continuous(
        name = expression(paste(P[0], " [mol ", m^-2, s^-1, "]"), sep=""),
        breaks = 10^x_range, expand = c(0,0), trans = "log10",
        labels = sapply(x_range, function (i) as.expression(bquote(10^ .(i))))
    ) + 
    scale_y_continuous(
        name = expression(paste("c(-w/2) [", K[d], "]"), sep=""),
        breaks = 10^y_range, expand = c(0,0), trans = "log10",
        labels = sapply(y_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    annotation_logticks(sides = "bl", color = "white") +
    scale_fill_gradientn(
        name = "CI", 
        colours = rev(c("#661f19", "#8c2a22", "#e9622e", "#f7e92b", 
                        "#86c98e", "#00b2ed", "#2061ae", "#303592")),
        breaks = seq(0, 1, 0.2), limits = c(0, 1)
    ) +
    theme(
        text = element_text(size=plot_text_size, color="black"),
        axis.text = element_text(size=plot_text_size, color="black"),
        axis.title = element_text(size=plot_text_size, color="black"),
        axis.ticks = element_blank(),
        plot.margin = margin(10, 10, 2, 2, unit="mm"),
        legend.position = c(1, 0.1),
        legend.text = element_text(size=legend_text_size, color="white"),
        legend.title = element_text(size=legend_text_size, color="white"),
        #legend.background = element_rect(fill=alpha("black", 0), color='white'),
        legend.background = element_blank(),
        legend.justification = c(1, 0),
        legend.key.size = unit(heat_legend_size, "mm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=alpha("blue", 0))
    )
#ci_heat.plot

# CI vertical slices


x_range = -2:4
ci_v.plot = ggplot(camp.dt[p0 %in% p0_range], aes(x=c0_left_over_kd, y=ci, group=factor(p0), color=factor(p0))) +
    geom_line(size=1.5) + 
    scale_y_continuous(
        name = "CI", 
        breaks = seq(0, 1, 0.2)
    ) + 
    scale_x_continuous(
        expression(paste("c(-w/2) [", K[d], "]")),
        breaks = 10^x_range, expand = c(0,0), trans = "log10",
        labels = sapply(x_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    annotation_logticks(color="black", sides="b") +
    scale_color_manual(
        name = expression(paste('c(-w/2) [', K[d], ']')),
        values = rainbow(length(c0_values), end=0.8),
        labels = sapply(round(log10(c0_values), 1), function (i) as.expression(bquote(10^ .(i)))),
        guide = 'colourbar'
    ) +
    guides(color=guide_legend(ncol=1, byrow=F, override.aes = list(size=10))) +
    theme(
        text = element_text(size=plot_text_size, color="black"),
        axis.text = element_text(size=plot_text_size, color="black"),
        axis.title = element_text(size=plot_text_size, color="black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color="gray89", linetype="dashed"),
        panel.border = element_rect(colour = "black", fill=alpha("blue",0)),
        panel.background = element_rect(fill="gray"),
        plot.margin = margin(10, 40, 2, 2, unit="mm"),
        legend.text = element_text(size=legend_text_size, color='black'),
        legend.position = c(1, 1.13),
        legend.background = element_blank(),
        legend.justification = c(0, 1),
        legend.key = element_blank(), 
        legend.key.size = legend_key_size,
        legend.key.height = legend_key_ht,
        legend.direction = "vertical",
        #legend.title = element_text(hjust=1, margin=margin(r=200, unit='mm')),
        legend.title.align = -4e7,
        legend.text.align = 0
    )    
#ci_v.plot


# This generates figure 2
pdf('../figures/fig2_pde_flux_plots.pdf', useDingbats = F, width=14, height=17.5)
multiplot(snr.plot, snr_h.plot, snr_v.plot, fwhm.plot, ci_heat.plot, ci_v.plot, 
          cols=2, annotate=T, vjust=0, hjust=-0.22, t=-5)
dev.off()