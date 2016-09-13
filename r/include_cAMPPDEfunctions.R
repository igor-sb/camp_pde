charLength <- function(Km, Dc, k2, P0) {
  # Calculate characteristic length for uniform cAMP-PDE model
  # 
  # Km = Michaelis menten constant for cAMP-PDE binding
  # Dc = cAMP diffusion coef.
  # k2 = 5'AMP synthesis rate
  # P0 = PDE concentration 
  return(sqrt(Km*Dc/(k2*P0)))
}




cAMPConc <- function(cL, w, L, r, x) {
  #
  # Calculates average of cAMP concentration over all combinations
  # of theta and phi which are just vectors.
  # --------------------------------------------------------------
  # cL = concentration at the left end
  # w  = channel width
  # L  = characteristic length scale (depends on PDE conc.)
  #    = sqrt(Km*Dc/(k2*P0))
  #    Km = 0.01 (mol/m^3)
  #    Dc = 4.44e-10 (m^2/s)
  #    k2 = 13300 (s^-1)
  # r  = cell radius; 5e-6 (5um)
  # x <- r*sin(theta)*cos(phi)
  
  return(mean(cL*sinh((w/(2*L) - x/L))/(sinh(w/L))))
}


cAMPConcAvg <- function(cL, w, L, r, theta, phi) {
  #
  # Calculates average of cAMP concentration over all combinations
  # of theta and phi which are just vectors.
  # --------------------------------------------------------------
  # cL = concentration at the left end
  # w  = channel width
  # L  = characteristic length scale (depends on PDE conc.)
  #    = sqrt(Km*Dc/(k2*P0))
  #    Km = 0.01 (mol/m^3)
  #    Dc = 4.44e-10 (m^2/s)
  #    k2 = 13300 (s^-1)
  # r  = cell radius; 5e-6 (5um)
  # theta,phi = angles at which to evaluate cAMP conc
  #
  # Note we do not want to average concentrations. Instead, we want
  # to average the c/(c+Kd) to average the receptor numbers between
  # front and the back of the cell.
  #
  grid <- expand.grid(theta, phi)
  colnames(grid) <- c("theta", "phi")
  grid$x <- r*sin(grid$theta)*cos(grid$phi)
  grid$cAMPConc <- cL*sinh((w/(2*L) - grid$x/L))/(sinh(w/L))
  # x <- r*sin(theta)*cos(phi)
  
  return(mean(grid$cAMPConc))
}

probRecOcc <- function(cL, w, L, r, Kd, x) {
  #
  # Calculate the probability of receptor occupancy in the case of
  # constant and uniform PDE at position x, where x=0 denotes the
  # position in the middle of the device.
  # --------------------------------------------------------------
  # cL = concentration at the left end
  # w  = channel width
  # L  = characteristic length scale (depends on PDE conc.)
  #    = sqrt(Km*Dc/(k2*P0))
  #    Km = 0.01 (mol/m^3)
  #    Dc = 4.44e-10 (m^2/s)
  #    k2 = 13300 (s^-1)
  # Kd = cAMP-receptor dissociation constant
  # r  = cell radius; 5e-6 (5um)
  # x <- r*sin(theta)*cos(phi)
  conc <- cL*sinh((w/(2*L) - x/L))/(sinh(w/L))
  return(conc/(conc+Kd))
}

probRecOccAvg <- function(cL, w, L, r, Kd, theta, phi) {
  #
  # ---------------------------------------------------------------
  # Calculates average receptor occupancy over different positions
  # on a sphere, in the case of uniform PDE model:
  #   < c(theta, phi) / ( c(theta, phi) + Kd ) >
  # where the average is taken over different {theta, phi}
  # ---------------------------------------------------------------
  #
  # cL = concentration at the left end
  # w  = channel width
  # L  = characteristic length scale (depends on PDE conc.)
  #    = sqrt(Km*Dc/(k2*P0))
  #    Km = 0.01 (mol/m^3)
  #    Dc = 4.44e-10 (m^2/s)
  #    k2 = 13300 (s^-1)
  # r  = cell radius; 5e-6 (5um)
  # Kd = dissociation constant for cAMP-cAMPreceptor binding
  # theta,phi = angles at which to evaluate cAMP conc
  #
  # ---------------------------------------------------------------
  # Note we do not want to average concentrations. Instead, we want
  # to average the c/(c+Kd) to average the receptor numbers between
  # front and the back of the cell.
  #
  grid <- expand.grid(theta, phi)
  colnames(grid) <- c("theta", "phi")
  grid$x <- r*sin(grid$theta)*cos(grid$phi)
  grid$cAMPConc <- cL*sinh((w/(2*L) - grid$x/L))/(sinh(w/L))
  grid$probRecOcc <- grid$cAMPConc/(grid$cAMPConc + Kd)
  # x <- r*sin(theta)*cos(phi)
  
  return(mean(grid$probRecOcc))
}



SNR <- function(pf, pb, cMean, Kd, I, Sb, Rt) {
  #
  # pf = prob of occ rec at the front
  # pb = prob of occ rec at the back
  # 
  # Kd = dissociation const.
  # I  = fold change
  # Sb = nonreceptor noise (sigma_B)
  #
  
  # Receptor numbers front and back
  Rf <- 0.5*Rt*pf
  Rb <- 0.5*Rt*pb
  
  numerator <- Rf - Rb
  
  # sigma_R
  # can we calculate this with just pf,pb? I don't think so, since we need
  # cf and cb. (or in this case we can just use cMean since the sigmaR^2 )
  # depends only on cMean in the first order (Eq.7 in manuscript)
  SigmaR2 <- Rt*cMean*Kd/(cMean+Kd)^2
  
  denum <- sqrt(SigmaR2/I + Sb^2)
  
  return(numerator/denum)
}

lambda0 <- function(j0, r0, Dc, Kd) {
  return( (j0*r0^2)/(Dc*Kd) )
}

lambda <- function(j0, r0, Dc, Kd, Km, k2, P0) {
  L <- charLength(Km, Dc, k2, P0)
  return( 
    lambda0(j0, r0, Dc, Kd)*(L/(r0 + L))*exp(r0/L) 
  )
}

cAMPConcSpherical <- function(r, j0, r0, Dc) {
  return(
    (j0*r0^2)/(Dc*r)
    # (lambda*Kd*exp(-r/L))/r
  )
}

cAMPGradSpherical <- function(r, lambda, Kd, L) {
  return(
    cAMPConcSpherical(r, lambda, Kd, L)*(1/r + 1/L)
  )
}

ExponentialLengthUm <- function(P0, Km, Dc, k2) {
  # Exponential length in micrometers
  # P0 is the input PDE concentration in nM
  return( 
    sqrt((Km*Dc)/(k2*P0*1e-6))*1e6
  )
}

SNRSpherical3DExact <- function(r, P0) {
  #
  # For specific, (r, p0) give back the single value of SNR
  #
  # r = distance to the center where the measuring cell would be
  # r0 = cell diameter
  # Kd = cAMP-cAR1 dissociation const.
  lambda <- lambda(j0, r0, Dc, Kd, Km, k2, P0)
  L <- charLength(Km, Dc, k2, P0)
  
  pb <- cAMPConcSpherical(r+r0/2, lambda, Kd, L)/(cAMPConcSpherical(r+r0/2, lambda, Kd, L) + Kd)
  pf <- cAMPConcSpherical(r-r0/2, lambda, Kd, L)/(cAMPConcSpherical(r-r0/2, lambda, Kd, L) + Kd)
  deltaR <- Rt*(pf - pb)/2
  
  SigmaR2 = Rt*pf*(1-pf) + Rt*pb*(1-pb)
  denum <- sqrt(SigmaR2/Is + Sb^2)
  
  return(deltaR/denum)
}

SNRShallowGrad <- function(Nrec, Kd, deltaC, cMean) {
  return( 0.5*sqrt(Kd*Nrec)*deltaC/(sqrt(cMean)*(cMean + Kd)) )
}
deltaCSpherical <- function(lambda, Kd, r0, r, L) {
  return( lambda*Kd*r0*(r+L)*exp(-r/L)/(L*r^2) )
}
cSpherical <- function(lambda, Kd, r, L) {
  return( lambda*Kd*exp(-r/L)/r )
}
cPointSourceInKd <- function(lambda, r) {
  return(lambda/r)
}



SigmaRSq <- function(Nrec, c, Kd) {
  return( Nrec*c*Kd/(c+Kd)^2 )
}

RMean <- function(Nrec, c, Kd) {
  return( Nrec*c/(c+Kd) )
}

SNRExact <- function(Rf, Rb, SigmaR, SigmaB, I) {
  return( (Rf-Rb)/sqrt((1/I)*SigmaR^2 + SigmaB^2) )
}

SNRSpherical3DExactDataFrame <- function(r, P0, r0, Kd, Km, Dc, k2, j0, Nrec, Sb, I) {
  #
  # For arrays r and p0 give back the data fram
  # r(um) | P0 [nM] | SNR
  # ...   | ...     | ...
  #
  # r = distance to the center where the measuring cell would be
  # r0 = cell diameter
  # Kd = cAMP-cAR1 dissociation const.
}

test <- function(r) {
  return(r*Dc)
}


theme_basic <- function(plotTextSize=14, base_family="Helvetica") {
  theme(
    text = element_text(size=plotTextSize, color="black"),
    axis.text = element_text(size=plotTextSize, color="black"),
    axis.title = element_text(size=plotTextSize, color="black"),
    legend.text = element_text(size=plotTextSize, color="black"),
    legend.title = element_text(size=plotTextSize, color="black"),
    #legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(colour = "black", fill=alpha("blue",0))   
  )
}

