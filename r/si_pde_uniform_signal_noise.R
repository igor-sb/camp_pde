#
#   Uniform PDE model, analytical solution. How does PDE change gradient, mean concentration?
#

# Clear all variables
rm(list=ls(all=T))

library(ggplot2)
library(data.table)
library(directlabels)

# Set working directory
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

source('include_multiplot.R')

plot_text_size = 25
legend_text_size = 19
snr_max = 20

# c0_values  = 10^seq(-2, 5, 0.5)

# x  = 0.5*w*sqrt(k2*p0/(Dc*Km))
rc = 5e-6      # cell radius (m)
Kd = 3e-5     # (in units mol/m^3) 30 nM Kd for cAMP
Rt = 70000    # total number of cAR1 receptors per cell
Is = 1.4      # sampling fold (1.8 from van Haastert, Biophys J., 93, 2007
Sb = 73       # non-receptor noise, (25 from van Haastert ^)
w  = 1e-3     # width of our microfluidic device in meters (1mm wide in x dir)
Km = 0.01     # 10 uM in mol/m^3, the data file says 1uM but it looks like a typo.
k2 = 13300    # s^-1
Dc = 4.44e-10 # [m^2/s] cAMP diffusion

L = function (p0) sqrt((Km*Dc)/(k2*p0))
p = function (c) c/(c+Kd)
cx = function (x, p0, cL) cL*sinh(0.5*w/L(p0) - x/L(p0))/sinh(w/L(p0))
cx_lin = function (x, p0, cL) 0.5*cL/cosh(0.5*w/L(p0)) - x*0.5*cL/(L(p0)*sinh(0.5*w/L(p0)))

dcdx_ratio_fun = function (p0) 0.5*w/(L(p0)*sinh(0.5*w/L(p0)))
c0_ratio_fun = function (p0) 1/cosh(0.5*w/L(p0))
dcdx_over_c0_fun = function (p0) dcdx_ratio_fun(p0)/c0_ratio_fun(p0)
snr_rel_f = function (p0, c) 0.5*sqrt(cosh(0.5*w/L(p0)))*w*(c+1)/(sinh(0.5*w/L(p0))*(1 + c/cosh(0.5*w/L(p0)))*L(p0))

# Generate the data using analytical equations
c0s = 10^seq(-2, 5, 0.1)*Kd
p0s = 10^seq(-11, -6, length.out = 100)

snr.dt = data.table(expand.grid(c0s, p0s))

names(snr.dt) = c('c0', 'p0')
snr.dt[, c_front := cx_lin(-0.5*rc, p0, c0*2), by=names(snr.dt)]
snr.dt[, c_back := cx_lin(0.5*rc, p0, c0*2), by=names(snr.dt)]
snr.dt[, delta_r := 0.5*Rt*(p(c_front) - p(c_back))]
snr.dt[, sigma_delta_r := sqrt(Sb^2 + 0.5*Rt*p(c_front)*(1-p(c_front))/Is + 0.5*Rt*p(c_back)*(1-p(c_back))/Is)]
snr.dt[, snr := delta_r/sigma_delta_r]
snr.dt[, c0_over_kd := c0/Kd]
snr.dt[, snr_rel := snr/.SD[p0==1e-11, snr]]
snr.dt[, snr_rel_2 := snr_rel_f(p0, c0_over_kd)]

dcdx.dt = data.table(p0=p0s, dcdx_ratio=dcdx_ratio_fun(p0s))
c0.dt = data.table(p0=p0s, c0_ratio=c0_ratio_fun(p0s))
dcdxoverc0 = data.table(p0=p0s, dcdx_over_c0=dcdx_over_c0_fun(p0s))
c_all.dt = Reduce(function (dt1, dt2) merge(dt1, dt2, by='p0'), list(dcdx.dt, c0.dt, dcdxoverc0))
c_all.dt = melt(c_all.dt, id=1, measure=2:4)

# Plot the gradient , relative gradient and mean concentration as a function of PDE concentration
x_range = (-11:-6)+6
y_range = -3:2
c_all.plot = ggplot(c_all.dt, aes(x=p0*1e6, y=value, group=variable, linetype=variable)) +
    geom_line(mapping = aes(), size = 1.5) +
    scale_x_continuous(
        name = expression(paste(p[0], ' [nM]', sep='')),
        trans = 'sqrt', 
        expand = c(0, 0)
        # breaks = 10^x_range,
        # labels = sapply(x_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    scale_y_continuous(
        name = '', trans = 'log10',
        breaks = 10^y_range,
        labels = sapply(y_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    scale_linetype_manual(
        values = c('solid', 'dashed', 'dotted'),
        labels = c( expression(
                        over( paste('|', nabla, 'c|'), 
                              paste('|', nabla, 'c|'[paste(p[0], '=0')])
                        )
                    ), 
                    expression(
                        over( paste(c[0], ' '),
                              c[paste('0, ', p[0], '=0')]
                        )
                    ),
                    expression(
                        over( paste('|', nabla, 'c|/', c[0]),
                              paste('(|', nabla, 'c|/', c[0], ')'[paste(p[0], '=0')])
                        )
                    )
                  )
    ) +
    coord_cartesian(ylim=c(1e-3, 1e2), xlim=c(0, 0.15)) +
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
        plot.margin = margin(5, 10, 2, 2, unit='mm'),
        plot.background = element_blank(), 
        legend.position = c(0, 0),
        # legend.background = element_rect(fill=alpha('white', 1)),
        legend.background = element_blank(), 
        legend.key = element_rect(color=NA, fill=alpha('black', 0)),
        legend.justification = c(0, 0),
        legend.key.size = unit(11, 'mm'),
        legend.key.width = unit(20, 'mm'),
        legend.key.height = unit(25, 'mm'),
        legend.direction = 'vertical',
        legend.text.align = 0,
        legend.text = element_text(size=1.1*plot_text_size),
        legend.title = element_blank()
    )      

# c_all.plot




# Absolute receptor difference 

c0_plot = snr.dt[, sort(unique(c0_over_kd))]
c0_plot = c0_plot[seq(1, length(c0_plot), 5)]
y_range = -1:4
x_range = -5:0

rec_dif.plot = ggplot(snr.dt[c0_over_kd %in% c0_plot], 
        aes(x=p0*1e6, y=delta_r, group=factor(c0_over_kd), color=factor(c0_over_kd))) +
    geom_line(size=1.5) +
    scale_x_continuous(name = expression(paste(p[0], ' [nM]', sep='')), trans='log10',
                       breaks = 10^x_range, expand = c(0, 0),
                       labels = sapply(x_range, function (i) as.expression(bquote(10^ .(i))))) +
    scale_y_continuous(name = expression(paste('signal ', Delta, R)), trans = 'log10',
                       breaks = 10^y_range, expand = c(0, 0),
                       labels = sapply(y_range, function (i) as.expression(bquote(10^ .(i))))
                       
    ) +
    coord_cartesian(ylim=c(1e-1, 3e3)) +
    scale_color_manual(
        name = expression(paste(c[0], ' [', K[d], ']', sep='')),
        values = rainbow(length(c0_plot), end=0.8),
        labels = sapply(round(log10(c0_plot), 1), function (i) as.expression(bquote(10^ .(i)))),
        breaks = c0_plot
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



# Plot absolute noise

sigma.plot = ggplot(snr.dt[c0_over_kd %in% c0_plot], aes(x=p0*1e6, y=sigma_delta_r, group=factor(c0_over_kd), color=factor(c0_over_kd))) +
    geom_line(size=1.5) +
    scale_x_continuous(name = expression(paste(p[0], ' [nM]', sep='')), trans='log10',
                       breaks = 10^x_range, expand = c(0, 0),
                       labels = sapply(x_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    scale_y_continuous(
        name = expression(paste('noise ', sqrt(sigma[B]^2 + sigma[paste(Delta, 'R')]^2 / I))),
        breaks = seq(0, floor(sqrt(Sb^2 + 0.25*Rt/Is)), 25)
        # expand = c(0, 0)
    ) +
    scale_color_manual(
        name = expression(paste(c[0], ' [', K[d], ']', sep='')),
        values = rainbow(length(c0_plot), end=0.8),
        labels = sapply(round(log10(c0_plot[6:13]), 1), function (i) as.expression(bquote(10^ .(i)))),
        breaks = c0_plot[6:13]
    ) +
    # coord_cartesian(ylim=c(0, 100)) +
    annotation_logticks(sides='b') +
    guides(color=guide_legend(ncol=2, byrow=F, override.aes = list(size=10))) +
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
# sigma.plot




# SNR heatmap

x_range = -5:0
y_range = -2:5
snr.plot = ggplot(snr.dt, aes(y = c0_over_kd, x = p0*1e6, z = snr, fill = snr)) +
    geom_raster(show.legend = T) +
    stat_contour(mapping=aes(color=..level..), color='black', size=0.7, breaks=c(1, 2, 5)) +
    scale_y_continuous(
        name = expression(paste(c[0], ' [', K[d], ']', sep='')), trans='log10', expand = c(0, 0),
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
y_range = -3:1

snr_hor.plot = ggplot(snr.dt[c0_over_kd %in% c0_plot], aes(x=p0*1e6, y=snr, group=factor(c0_over_kd), color=factor(c0_over_kd))) +
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
        name = expression(paste(c[0], ' [', K[d], ']', sep='')),
        values = rainbow(length(c0_plot), end=0.8),
        labels = sapply(round(log10(c0_plot), 1), function (i) as.expression(bquote(10^ .(i)))),
        breaks = c0_plot
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
x_range = -5:0
y_range = -2:5
snr_ratio_range = -4:4

snr_ratio.plot = ggplot(snr.dt[snr_rel > 1e-4], aes(y=c0_over_kd, x=p0*1e6, z=snr_rel, fill=snr_rel)) +
    geom_raster(show.legend = T) +
    stat_contour(mapping=aes(color=..level..), color='black', size=0.7, breaks=c(2, 10, 100)) +
    scale_y_continuous(
        name = expression(paste(c[0], ' [', K[d], ']', sep='')), trans='log10', expand = c(0, 0),
        breaks = 10^y_range,
        labels = sapply(y_range, function (i) as.expression(bquote(10^ .(i))))
    ) +
    scale_x_continuous(
        name = expression(paste(p[0], ' [nM]')), trans='log10', expand = c(0, 0),
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
        plot.margin = margin(10, 10, 2, 2, unit='mm'),
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
#c0_plot = c0_range[1:61][seq(1, 61, 5)]

snr_ratio_hor.plot = ggplot(snr.dt[c0_over_kd %in% c0_plot], aes(x=p0*1e6, y=snr_rel, group=factor(c0_over_kd), color=factor(c0_over_kd))) +
    geom_line(size=1.5) +
    scale_x_continuous(
        name = expression(paste(p[0], ' [nM]', sep='')), trans='log10',
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
        labels = sapply(round(log10(c0_plot), 1), function (i) as.expression(bquote(10^ .(i)))),
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

# This generates figure 3
pdf('../figures/si_pde_uniform_signal_noise_snr.pdf', useDingbats=F, width=14, height=17.5)
multiplot(snr.plot, snr_hor.plot, c_all.plot, rec_dif.plot, sigma.plot, snr_ratio.plot, 
          cols=2, annotate=T, vjust=0, hjust=-0.20, t=-5)
dev.off()

# Export 14 by 17.5 inch
# SNR ratio horizontal slices (paramterize by c0

