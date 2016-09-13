#
#  Plots from van Haastert Postma Biophys J 2007, CI(SNR) and cost minimization
#

# Clear all variables
rm(list=ls(all=T))

# Load libraries
library(data.table)
library(ggplot2)

# van Haastert / Postma parameters
Sb = 25
Is = 1.8
Rt_vp = 40000
Kd_vp = 100        # nM
rc_vp = 5 # um

# Ours parameters
# Load model parameters and import into variables
params.data <- read.table('../data/model_parameters.txt', header=T, sep='\t', row.names=1)
Kd <- params.data['Kd',]$value
Rt <- params.data['Rt',]$value
rc <- params.data['rc',]$value

# Plot parameters
plot_text_size = 22

# Set working directory
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# Include multiplot
source('include_multiplot.R')

# Load clean data from Excel after extracting numerical values from their plots
ci.dt = fread('../data/ci_dcdx_c0_from_experiments.txt')

# First repeat their Figure 1 a and b to see if everything checks out 
# fig1.plot = ggplot(ci.dt[experiment %like% 'pipette'], 
#                    aes(x=distance_from_pipette_um, y=CI, group=factor(cAMP_pipette_uM), shape=factor(cAMP_pipette_uM))) +
#     geom_line(size=0.7) + 
#     geom_point(size=4) +
#     scale_shape_manual(values = c('100'=19, '10'=18, '1'=17, '0.1'=15)) +
#     scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2))
# fig1.plot
# 
# fig2.plot = ggplot(ci.dt[!(experiment %like% 'pipette') & experiment %like% 'van haastert'], 
#                    aes(x=cAMP_nM, y=CI, group=factor(experiment), shape=factor(experiment))) +
#     geom_line(size=0.7) + 
#     geom_point(size=4) +
#     scale_shape_manual(values = c('van haastert postma bridge'=19, 'van haastert postma population 650 um'=17, 'van haastert postma population 250 um'=15)) +
#     scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2)) +
#     scale_x_continuous(trans="log10", breaks=10^(0:5), limits=c(1, 1e5))
# fig2.plot

# Everything looks correct. Now calculate their SNR, try to replicate their plots.
erf = function (x) 2 * pnorm(x * sqrt(2)) - 1


ci.dt[, snr_vp := (0.25*Rt_vp*2*rc_vp*dcdx_nM_um/Kd_vp)/sqrt(Sb^2 + (c0_nM*Rt_vp)/(Kd_vp*Is))]
ci.dt[, Ppos_vp := erf(snr_vp/sqrt(2))]
ci.dt[, experiment_finer := experiment]
ci.dt[experiment %like% 'pipette', experiment_finer := paste0(experiment, ' ', as.character(cAMP_pipette_uM), ' uM')]



p = function (c) c/(c+Kd)
dR = function (cf, cb) 0.5*Rt*(p(cf) - p(cb))
noise = function (cf, cb) sqrt(Sb^2 + 0.5*Rt*p(cf)*(1-p(cf))/Is + 0.5*Rt*p(cb)*(1-p(cb))/Is)
log_inv = function (x) log(x/(1-x))

p_vp = function (c) c/(c+Kd_vp)
snr_vp_f = function (c_front, c_back) 0.5*Rt_vp*(p_vp(c_front) - p_vp(c_back))/sqrt(Sb^2 + 0.5*Rt_vp*p_vp(c_front)*(1-p_vp(c_front))/Is + 0.5*Rt_vp*p_vp(c_back)*(1-p_vp(c_back))/Is)
snr_f = function (c_front, c_back) dR(c_front, c_back) / noise(c_front, c_back)

ci.dt[experiment %like% 'van Haastert', snr_vp_exact := snr_vp_f(c0_nM + 0.5*rc_vp*dcdx_nM_um, c0_nM - 0.5*rc_vp*dcdx_nM_um)]
ci.dt[experiment %like% 'van Haastert', Ppos_vp_exact := erf(snr_vp_exact/sqrt(2))]
ci.dt[, snr := snr_f(c0_nM*1e-6 + 0.5*rc*dcdx_nM_um, c0_nM*1e-6 - 0.5*rc*dcdx_nM_um)]
ci.dt[, log_snr_sqrt2 := log10(snr/sqrt(2))]
ci.dt[, Ppos := erf(snr/sqrt(2))]
setorder(ci.dt, 'experiment_finer')

# CI vs Ppos using van Haastert's own parameters and using exact SNR. not his approximation.
ci_vs_ppos.plot = ggplot(ci.dt[experiment %like% 'van Haastert'], 
    aes(x=Ppos_vp_exact, y=CI, group=experiment_finer, 
        shape=gsub('van Haastert, Postma; ', '', experiment_finer, fixed=T))) +
    geom_point(size=4) +
    scale_shape_manual(values = c(8, 18, 15, 17, 19, 1, 2)) +
    scale_x_continuous(name='Ppos', limits=c(0,1), breaks=seq(0,1,0.2)) +
    scale_y_continuous(name='CI', limits=c(0,1), breaks=seq(0,1,0.2)) +
    theme(
        text = element_text(size=plot_text_size, color='black'),
        legend.text = element_text(size=plot_text_size*0.7),
        axis.text = element_text(color='black'),
        plot.margin = margin(10, 10, 2, 2, unit="mm"),
        legend.title = element_blank()
        #legend.position = c(0.5, 1),
        #legend.justification = c(0.5, 1)
    )
pdf('van_haastert_postma_replot_ci_ppos.pdf', width=6.9, height=4.5, useDingbats = F)
ci_vs_ppos.plot
dev.off()

# This RESEMBLES their plot but 
# 1. population assay points are a bit off. The equation they gave for gradient 
#    produces data that is not consistent with what they have even in the table.
# 2. Ppos > 0.95 is now included so the data is not cherry picked.


# ci_vs_snr.plot

#
#  To find CI vs Ppos best fit. Calculate Ppos for each data point c0 dcdx, for a combination of parameters Is, Sb
#  For specific Is,Sb, calculate sum( (CI-Ppos)^2 ) where sum is over all points (c0, dcdx) and then show this cost
#  as a function of Sb, Is and find the minimum.
#

Sb_range = seq(0, 300, 1)
Is_range = seq(1, 10, 0.2)

ci.sdt = ci.dt[year > 1989 & !(experiment_finer %like% 'population' | experiment_finer %like% 'bridge' | experiment_finer %like% 'Fuller')]

pfit.dt = data.table(expand.grid(Sb_range, Is_range))
names(pfit.dt) = c('Sb', 'Is')

snr_custom = function (c0, dcdx, Sb_, Is_) ( 0.5*Rt*(p(c0*1e-6 + 0.5*rc*dcdx) - p(c0*1e-6 - 0.5*rc*dcdx)) ) /
    sqrt(Sb_^2 + 0.5*Rt*p(c0*1e-6 + 0.5*rc*dcdx)*(1-p(c0*1e-6 + 0.5*rc*dcdx))/Is_ + 0.5*Rt*p(c0*1e-6 - 0.5*rc*dcdx)*(1-p(c0*1e-6 - 0.5*rc*dcdx))/Is_)

mm = function (x) x/(x+1)

pfit.dt[, erf_cost := ci.sdt[, sum( ( CI - erf(snr_custom(c0_nM, dcdx_nM_um, Sb, Is)/sqrt(2)) )^2,  na.rm=T)], by=c('Sb', 'Is')]
pfit.dt[, mm_cost := ci.sdt[, sum( ( CI - mm(snr_custom(c0_nM, dcdx_nM_um, Sb, Is)/sqrt(2)) )^2,  na.rm=T)], by=c('Sb', 'Is')]

# Find Sb and Is that minimize the fit CI(SNR) = snr/(snr + 1)
Sb = pfit.dt[mm_cost==min(mm_cost)]$Sb
Is = pfit.dt[mm_cost==min(mm_cost)]$Is
ci.sdt[, snr := snr_f(c0_nM*1e-6 + 0.5*rc*dcdx_nM_um, c0_nM*1e-6 - 0.5*rc*dcdx_nM_um)]
ci.sdt[, Ppos := erf(snr/sqrt(2))]
ci.sdt[, Ppos2 := Ppos^2]
ci.sdt[, Ppos3 := Ppos^3]
# ci.sdt[, CI_fit_mm := snr/(snr + summary(nls(formula = CI ~ snr/(a + snr), start = c(a=1)))$parameters[1])]
cubic_fit_params = t(ci.sdt[, summary(lm(CI ~ Ppos + Ppos2 + Ppos3))$coefficients[,1]])
quad_fit_params = t(ci.sdt[, summary(lm(CI ~ Ppos + Ppos2))$coefficients[,1]])
setorder(ci.sdt, 'experiment_finer')
ci_fit.dt = data.table(Ppos=seq(0,1, length.out=50))
ci_fit.dt$CI_fit_cubic = as.vector(cubic_fit_params %*% t(as.matrix(ci_fit.dt[, list(Ones=1, Ppos=Ppos, Ppos2=Ppos^2, Ppos3=Ppos^3)])))
ci_fit.dt$CI_fit_quad = as.vector(quad_fit_params %*% t(as.matrix(ci_fit.dt[, list(Ones=1, Ppos=Ppos, Ppos2=Ppos^2)])))

# Plot best fit values
# Is_values = pfit.dt[, unique(Is)][seq(1, 46, 2)]
Is_values = c(1, 1.4, pfit.dt[, unique(Is)][seq(6, 46, 5)])
sbis.plot = ggplot(pfit.dt[Is %in% Is_values], aes(x=Sb, y=mm_cost, color=factor(Is))) +
    geom_line(size=1.5) +
    scale_x_continuous(
        name = expression(sigma[B]),
        breaks = seq(0, 150, 50),
        expand = c(0, 0)
    ) +
    scale_y_continuous(
        name = expression(paste('S(', sigma[B], ', ', I, ')')),
        breaks = seq(0.2, 1, 0.2), expand = c(0,0)
    ) +
    scale_color_manual(
        name = 'I',
        values = rainbow(length(Is_values), end=0.8),
        labels = sprintf('%.1f', Is_values)
    ) +
    coord_cartesian(ylim=c(0.2, 1), xlim=c(0, 200)) +
    guides(color=guide_legend(ncol=7, byrow=F, override.aes = list(size=10))) +
    theme(
        text = element_text(size=plot_text_size*1.2, color='black'),
        axis.text = element_text(color='black'),
        # axis.ticks.length = unit(1.5, 'mm'),
        # axis.ticks = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size=plot_text_size*0.75),
        legend.position = c(-0.025, -0.01),
        legend.justification = c(0, 0),
        legend.background = element_blank(),
        plot.margin = unit(c(0.6, 0.5, 0.5, 0.5), "lines"),
        panel.border = element_rect(colour='black', fill=NA),
        panel.background = element_rect(fill='gray'),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color='gray89', linetype='dashed')
    )

pdf('../figures/si_ci_snr_cost_min.pdf', width=6, height=6, useDingbats=F)
sbis.plot
dev.off()
# Export 6 by 6

sbis_heat.plot = ggplot(pfit.dt, aes(x=Sb, y=Is, fill=mm_cost)) +
    geom_raster() +
    scale_x_continuous(
        name = expression(paste(sigma[B])),
        expand = c(0, 0)
    ) +
    scale_y_continuous(
        name = 'I', breaks = 1:10,
        expand = c(0, 0)
    ) +
    scale_fill_gradientn(
        name = 'cost', limits = c(0.3, 1),
        breaks = c(0.1, 1, 10),
        colours = rev(c("#661f19", "#8c2a22", "#e9622e", "#f7e92b", 
                        "#86c98e", "#00b2ed", "#2061ae", "#303592"))
    )
# sbis_heat.plot

# Compare this to R's nonlinear fit nls()
# ci.sdt[, cf := c0_nM*1e-6 + 0.5*rc*dcdx_nM_um]
# ci.sdt[, cb := c0_nM*1e-6 - 0.5*rc*dcdx_nM_um]
# snr_2 = function (cf, cb, Sb_, Is_) 0.5*Rt*(cf/(cf+Kd) - cb/(cb+Kd)) / sqrt(Sb_^2 + 0.5*Rt*cf*Kd/(Is_*(cf+Kd)^2) + 0.5*Rt*cb*Kd/(Is_*(cb+Kd)^2))
# tmp = ci.sdt[, nls(CI ~ snr_2(cf, cb, sigmaB, is) / 
#                        (snr_2(cf, cb, sigmaB, is) + 1),
#                    start=c(is=1.5, sigmaB=50))]
# This throws out an error. 

# Export table to LaTeX
ci.sdt[, experiment_format := paste0(experiment, '(', year, ')')]
# Run this then copy/paste to Texniccenter
# xtable(ci.sdt[, .(c0_nM, dcdx_nM_um, CI, experiment_format)])

# CI vs SNR fits
ci_snr_fit.dt = data.table(snr=10^seq(-5, 3, length.out=50))
snr_km = ci.sdt[, summary(nls(formula = CI ~ snr/(a + snr), start = c(a=1)))]$parameters[1]
ci_snr_fit.dt[, ci_fit_mm := snr/(snr + 1)]

ci.sdt[, experiment_format := paste(gsub('; pipette', '\n', experiment, fixed=T), year, sep='')]
# ci.sdt[, experiment_format := factor(experiment_format, levels=rev(unique(experiment_format)))]



x_range = -4:2
ci_vs_snr.plot = ggplot(ci.sdt,
    aes(x=snr, y=CI, group=experiment_format, color=experiment_format)
    ) +
    geom_line(data=ci_snr_fit.dt, mapping=aes(x=snr, y=ci_fit_mm), color='black', inherit.aes=F, size=1.1) +
    geom_point(size=4) +
    coord_cartesian(xlim=c(1e-4, 1e2)) +
    annotation_logticks(sides='b') +
    scale_x_continuous(name='SNR', trans='log10', breaks=10^x_range,
                       labels = sapply(x_range, function (i) as.expression(bquote(10^ .(i))))) +
    scale_y_continuous(name='CI', limits=c(0,1), breaks=seq(0,1,0.2)) +
    scale_color_discrete() +
    theme(
        text = element_text(size=plot_text_size*1.2, color='black'),
        axis.text = element_text(size=plot_text_size, color='black'),
        legend.text = element_text(size=plot_text_size*0.75, color='black', vjust=1),
        legend.title = element_blank(),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.key.height = unit(11, 'mm'),
        panel.grid.minor = element_blank()
    )


ci_vs_snr.plot2 = ggplot(
        ci.sdt,
        aes(x=mm(snr), y=CI, group=experiment_format, color=experiment_format)
    ) +
    geom_line(data=ci_snr_fit.dt, mapping=aes(x=mm(snr), y=ci_fit_mm), color='black', inherit.aes=F, size=1.1) +
    geom_point(size=4) +
    coord_cartesian(xlim=c(0, 1)) +
    scale_x_continuous(name='SNR/(SNR + 1)', breaks=seq(0,1,0.2)) +
    scale_y_continuous(name='CI', limits=c(0,1), breaks=seq(0,1,0.2)) +
    scale_color_discrete() +
    theme(
        text = element_text(size=plot_text_size, color='black'),
        axis.text = element_text(size=plot_text_size, color='black'),
        legend.text = element_text(size=plot_text_size*0.8, color='black'),
        legend.title = element_blank(),
        # legend.position = c(0, 1),
        legend.position = 'none',
        legend.justification = c(0, 1),
        legend.key.height = unit(11, 'mm'),
        panel.grid.minor = element_blank()
    )
#ci_vs_snr.plot2

pdf('../figures/si_ci_vs_snr_mm_fit.pdf', width=11, height=5.5, useDingbats=F)
multiplot(ci_vs_snr.plot, ci_vs_snr.plot2, cols=2, annotate=T, vjust=0, hjust=-0.19, t=-2.5)
dev.off()
# Export 7.1 by 5.5

# CI vs Ppos fits
ci_vs_ppos.plot = ggplot(ci.sdt,
                         aes(x=Ppos, y=CI, group=factor(experiment_finer), color=factor(experiment_finer))) +
    geom_point(size=4) +
    geom_line(data = ci_fit.dt, mapping = aes(x=Ppos, y=CI_fit_cubic), color='red', inherit.aes = F) +
    scale_x_continuous(name='Ppos', limits=c(0,1), breaks=seq(0,1,0.2)) +
    scale_y_continuous(name='CI', limits=c(0,1), breaks=seq(0,1,0.2))
# ci_vs_ppos.plot

# Show all experiments CS vs Ppos
ci_vs_ppos.plot = ggplot(ci.dt[year > 2000 & !(experiment_finer %like% c('population'))], 
                         aes(x=Ppos, y=CI, group=factor(experiment_finer), 
    color=factor(experiment_finer), shape=factor(gradient_type))) +
    geom_point(size=4) +
    scale_x_continuous(name='Ppos', limits=c(0,1), breaks=seq(0,1,0.2)) +
    scale_y_continuous(name='CI', limits=c(0,1), breaks=seq(0,1,0.2))
# ci_vs_ppos.plot


