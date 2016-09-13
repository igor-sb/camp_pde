%
% FEM struct for constant PDE flux
%

% Use Chemical Engineering application mode, Diffusion
fem.appl.module     = 'CHEM';       % Select Chem.Eng. module
fem.appl.mode.class = 'Diffusion';  % Mass diffusion class from Chem.Eng. module
fem.appl.name       = 'chdi';
fem.appl.dim        = {'c', 'p';};
fem.appl.gporder    = 4;
fem.appl.cporder    = 2;

% fem.appl.sshape     = 2; 

fem.appl.assignsuffix   = '_chdi';
fem.appl.prop.analysis  = 'static'; 

% Reaction and diffusion terms:
fem.appl.equ.D      = {'Dc'; 'Dp';};
fem.appl.equ.name   = {'mainDomain'};
fem.appl.equ.R      = {'-k2*c*p/KM';};
fem.appl.equ.ind    = 1;

% Boundary conditions
fem.appl.bnd.c0     = {'c0L', 'c0R', 0, 0};
fem.appl.bnd.name   = {'leftWall', 'rightWall', 'cellBoundary', 'reflectiveBoundaries'};
fem.appl.bnd.type   = { ...
    {'C';'C';}, ...
    {'C';'C';}, ...
    {'N';'N';}, ...
    {'N';'N';}
    };
fem.appl.bnd.N = { ...
    0, 0, {0;'p0';}, 0
    };
fem.appl.bnd.ind = [1,4,4,4,4,3,3,3,3,2];  

% Constants
% These numbers will get rewritten in the for loop below
% fem.const = { ...
%     'k2','5000 [1/s]', ...
%     'KM','10 [umol/l]', ...
%     'p0','1.5e-13*1e2 [mol/(m^2*s)]', ...
%     'c0L','40 [nmol/l]', ...
%     'c0R','0 [nmol/l]', ...
%     'Dc','4.44e-10 [m^2/s]', ...
%     'Dp','7e-11 [m^2/s]' ...
% ;};

% Other stuff
fem.frame   = {'ref'};
fem.shape   = {'shlag(2,''c'')','shlag(2,''p'')';};
fem.border  = 1;
fem.outform = 'general';
fem.form    = 'general';

% Units
fem.units.basesystem = 'SI';
fem.units.vardim = {
    '3Dxx_p_chdi','diffusion coefficient','3Dxy_p_chdi','diffusion coefficient','3Dxz_c_chdi',...
    'diffusion coefficient','3D_c_chdi','diffusion coefficient','2kc_p_chdi','mass transfer coefficient',...
    '3Dyz_c_chdi','diffusion coefficient','3Dxy_c_chdi','diffusion coefficient','3Dyx_c_chdi',...
    'diffusion coefficient','3Dts_p_chdi',[0,0,0,0,0,0,0,0],'2c0_p_chdi','concentration','3R_c_chdi',...
    [-3,0,-1,0,0,0,1,0],'2ny_chdi',[0,0,0,0,0,0,0,0],'2c0_c_chdi','concentration','2kc_c_chdi',...
    'mass transfer coefficient','2Dbnd_p_chdi','diffusion coefficient','3Dxz_p_chdi',...
    'diffusion coefficient','2N_c_chdi','molar flux','3Dzz_c_chdi','diffusion coefficient',...
    '3Dzx_c_chdi','diffusion coefficient','3Dyx_p_chdi','diffusion coefficient','3Dzy_p_chdi',...
    'diffusion coefficient','2N_p_chdi','molar flux','2nx_chdi',[0,0,0,0,0,0,0,0],'2d_chdi','length',...
    '2nz_chdi',[0,0,0,0,0,0,0,0],'2cb_c_chdi','concentration','3Dyy_p_chdi','diffusion coefficient',...
    '3Dxx_c_chdi','diffusion coefficient','3Dts_c_chdi',[0,0,0,0,0,0,0,0],'3Dzy_c_chdi',...
    'diffusion coefficient','3Dzx_p_chdi','diffusion coefficient','3Dyz_p_chdi',...
    'diffusion coefficient','2Dbnd_c_chdi','diffusion coefficient','3Dyy_c_chdi',...
    'diffusion coefficient','3c','concentration','3Dzz_p_chdi','diffusion coefficient','3D_p_chdi',...
    'diffusion coefficient','3p','concentration','2cb_p_chdi','concentration','3R_p_chdi',...
    [-3,0,-1,0,0,0,1,0];};


% Equation forms for COMSOL
fem.equ.f{1,1} = {'R_c_chdi'; 0;};
fem.equ.c{1,1}{1,1} = {'-diff(-Dxx_c_chdi*cx-Dxy_c_chdi*cy-Dxz_c_chdi*cz,cx)','-diff(-Dxx_c_chdi*cx-Dxy_c_chdi*cy-Dxz_c_chdi*cz,cy)','-diff(-Dxx_c_chdi*cx-Dxy_c_chdi*cy-Dxz_c_chdi*cz,cz)';'-diff(-Dyx_c_chdi*cx-Dyy_c_chdi*cy-Dyz_c_chdi*cz,cx)','-diff(-Dyx_c_chdi*cx-Dyy_c_chdi*cy-Dyz_c_chdi*cz,cy)','-diff(-Dyx_c_chdi*cx-Dyy_c_chdi*cy-Dyz_c_chdi*cz,cz)';'-diff(-Dzx_c_chdi*cx-Dzy_c_chdi*cy-Dzz_c_chdi*cz,cx)','-diff(-Dzx_c_chdi*cx-Dzy_c_chdi*cy-Dzz_c_chdi*cz,cy)','-diff(-Dzx_c_chdi*cx-Dzy_c_chdi*cy-Dzz_c_chdi*cz,cz)';};
fem.equ.c{1,1}{1,2} = {'-diff(-Dxx_c_chdi*cx-Dxy_c_chdi*cy-Dxz_c_chdi*cz,px)','-diff(-Dxx_c_chdi*cx-Dxy_c_chdi*cy-Dxz_c_chdi*cz,py)','-diff(-Dxx_c_chdi*cx-Dxy_c_chdi*cy-Dxz_c_chdi*cz,pz)';'-diff(-Dyx_c_chdi*cx-Dyy_c_chdi*cy-Dyz_c_chdi*cz,px)','-diff(-Dyx_c_chdi*cx-Dyy_c_chdi*cy-Dyz_c_chdi*cz,py)','-diff(-Dyx_c_chdi*cx-Dyy_c_chdi*cy-Dyz_c_chdi*cz,pz)';'-diff(-Dzx_c_chdi*cx-Dzy_c_chdi*cy-Dzz_c_chdi*cz,px)','-diff(-Dzx_c_chdi*cx-Dzy_c_chdi*cy-Dzz_c_chdi*cz,py)','-diff(-Dzx_c_chdi*cx-Dzy_c_chdi*cy-Dzz_c_chdi*cz,pz)';};
fem.equ.c{1,1}{2,1} = {'-diff(-Dxx_p_chdi*px-Dxy_p_chdi*py-Dxz_p_chdi*pz,cx)','-diff(-Dxx_p_chdi*px-Dxy_p_chdi*py-Dxz_p_chdi*pz,cy)','-diff(-Dxx_p_chdi*px-Dxy_p_chdi*py-Dxz_p_chdi*pz,cz)';'-diff(-Dyx_p_chdi*px-Dyy_p_chdi*py-Dyz_p_chdi*pz,cx)','-diff(-Dyx_p_chdi*px-Dyy_p_chdi*py-Dyz_p_chdi*pz,cy)','-diff(-Dyx_p_chdi*px-Dyy_p_chdi*py-Dyz_p_chdi*pz,cz)';'-diff(-Dzx_p_chdi*px-Dzy_p_chdi*py-Dzz_p_chdi*pz,cx)','-diff(-Dzx_p_chdi*px-Dzy_p_chdi*py-Dzz_p_chdi*pz,cy)','-diff(-Dzx_p_chdi*px-Dzy_p_chdi*py-Dzz_p_chdi*pz,cz)';};
fem.equ.c{1,1}{2,2} = {'-diff(-Dxx_p_chdi*px-Dxy_p_chdi*py-Dxz_p_chdi*pz,px)','-diff(-Dxx_p_chdi*px-Dxy_p_chdi*py-Dxz_p_chdi*pz,py)','-diff(-Dxx_p_chdi*px-Dxy_p_chdi*py-Dxz_p_chdi*pz,pz)';'-diff(-Dyx_p_chdi*px-Dyy_p_chdi*py-Dyz_p_chdi*pz,px)','-diff(-Dyx_p_chdi*px-Dyy_p_chdi*py-Dyz_p_chdi*pz,py)','-diff(-Dyx_p_chdi*px-Dyy_p_chdi*py-Dyz_p_chdi*pz,pz)';'-diff(-Dzx_p_chdi*px-Dzy_p_chdi*py-Dzz_p_chdi*pz,px)','-diff(-Dzx_p_chdi*px-Dzy_p_chdi*py-Dzz_p_chdi*pz,py)','-diff(-Dzx_p_chdi*px-Dzy_p_chdi*py-Dzz_p_chdi*pz,pz)';};

fem.equ.a{1,1}      = {'-diff(R_c_chdi,c)','-diff(R_c_chdi,p)';0,0;};

fem.equ.ga{1,1}{1,1} = {'-Dxx_c_chdi*cx-Dxy_c_chdi*cy-Dxz_c_chdi*cz';'-Dyx_c_chdi*cx-Dyy_c_chdi*cy-Dyz_c_chdi*cz';'-Dzx_c_chdi*cx-Dzy_c_chdi*cy-Dzz_c_chdi*cz';};
fem.equ.ga{1,1}{2,1} = {'-Dxx_p_chdi*px-Dxy_p_chdi*py-Dxz_p_chdi*pz';'-Dyx_p_chdi*px-Dyy_p_chdi*py-Dyz_p_chdi*pz';'-Dzx_p_chdi*px-Dzy_p_chdi*py-Dzz_p_chdi*pz';};

fem.equ.al{1,1}{1,1} = {'-diff(-Dxx_c_chdi*cx-Dxy_c_chdi*cy-Dxz_c_chdi*cz,c)';'-diff(-Dyx_c_chdi*cx-Dyy_c_chdi*cy-Dyz_c_chdi*cz,c)';'-diff(-Dzx_c_chdi*cx-Dzy_c_chdi*cy-Dzz_c_chdi*cz,c)';};
fem.equ.al{1,1}{1,2} = {'-diff(-Dxx_c_chdi*cx-Dxy_c_chdi*cy-Dxz_c_chdi*cz,p)';'-diff(-Dyx_c_chdi*cx-Dyy_c_chdi*cy-Dyz_c_chdi*cz,p)';'-diff(-Dzx_c_chdi*cx-Dzy_c_chdi*cy-Dzz_c_chdi*cz,p)';};
fem.equ.al{1,1}{2,1} = {'-diff(-Dxx_p_chdi*px-Dxy_p_chdi*py-Dxz_p_chdi*pz,c)';'-diff(-Dyx_p_chdi*px-Dyy_p_chdi*py-Dyz_p_chdi*pz,c)';'-diff(-Dzx_p_chdi*px-Dzy_p_chdi*py-Dzz_p_chdi*pz,c)';};
fem.equ.al{1,1}{2,2} = {'-diff(-Dxx_p_chdi*px-Dxy_p_chdi*py-Dxz_p_chdi*pz,p)';'-diff(-Dyx_p_chdi*px-Dyy_p_chdi*py-Dyz_p_chdi*pz,p)';'-diff(-Dzx_p_chdi*px-Dzy_p_chdi*py-Dzz_p_chdi*pz,p)';};

fem.equ.be{1,1}{1,1} = {'-diff(R_c_chdi,cx)';'-diff(R_c_chdi,cy)';'-diff(R_c_chdi,cz)';};
fem.equ.be{1,1}{1,2} = {'-diff(R_c_chdi,px)';'-diff(R_c_chdi,py)';'-diff(R_c_chdi,pz)';};
fem.equ.be{1,1}{2,1} = {0;0;0;};
fem.equ.be{1,1}{2,2} = {0;0;0;};

fem.equ.ind = 1;
fem.equ.dim = {'c'; 'p';};
fem.equ.var = {'grad_c_x_chdi','cx','dflux_c_x_chdi','-Dxx_c_chdi*cx-Dxy_c_chdi*cy-Dxz_c_chdi*cz','grad_c_y_chdi','cy','dflux_c_y_chdi','-Dyx_c_chdi*cx-Dyy_c_chdi*cy-Dyz_c_chdi*cz','grad_c_z_chdi','cz','dflux_c_z_chdi','-Dzx_c_chdi*cx-Dzy_c_chdi*cy-Dzz_c_chdi*cz','grad_c_chdi','sqrt(grad_c_x_chdi^2+grad_c_y_chdi^2+grad_c_z_chdi^2)','dflux_c_chdi','sqrt(dflux_c_x_chdi^2+dflux_c_y_chdi^2+dflux_c_z_chdi^2)','grad_p_x_chdi','px','dflux_p_x_chdi','-Dxx_p_chdi*px-Dxy_p_chdi*py-Dxz_p_chdi*pz','grad_p_y_chdi','py','dflux_p_y_chdi','-Dyx_p_chdi*px-Dyy_p_chdi*py-Dyz_p_chdi*pz','grad_p_z_chdi','pz','dflux_p_z_chdi','-Dzx_p_chdi*px-Dzy_p_chdi*py-Dzz_p_chdi*pz','grad_p_chdi','sqrt(grad_p_x_chdi^2+grad_p_y_chdi^2+grad_p_z_chdi^2)','dflux_p_chdi','sqrt(dflux_p_x_chdi^2+dflux_p_y_chdi^2+dflux_p_z_chdi^2)','D_c_chdi','Dc','Dxx_c_chdi','Dc','Dyx_c_chdi',0,'Dzx_c_chdi',0,'Dxy_c_chdi',0,'Dyy_c_chdi','Dc','Dzy_c_chdi',0,'Dxz_c_chdi',0,'Dyz_c_chdi',0,'Dzz_c_chdi','Dc','D_p_chdi','Dp','Dxx_p_chdi','Dp','Dyx_p_chdi',0,'Dzx_p_chdi',0,'Dxy_p_chdi',0,'Dyy_p_chdi','Dp','Dzy_p_chdi',0,'Dxz_p_chdi',0,'Dyz_p_chdi',0,'Dzz_p_chdi','Dp','Dts_c_chdi',1,'R_c_chdi','-k2*c*p/KM','Dts_p_chdi',1,'R_p_chdi',0;};

% Boundary conditions
fem.bnd.shape = [1;2;];
fem.bnd.g = {0,0,{0;'N_p_chdi';},0};

fem.bnd.h{1,1} = {'-diff(-c+c0_c_chdi,c)','-diff(-c+c0_c_chdi,p)';'-diff(-p,c)','-diff(-p,p)';};
fem.bnd.h{1,2} = 0;
fem.bnd.h{1,3} = 0;
fem.bnd.h{1,4} = {'-diff(-c+c0_c_chdi,c)','-diff(-c+c0_c_chdi,p)';'-diff(-p,c)','-diff(-p,p)';};

fem.bnd.r{1,1} = {'-c+c0_c_chdi';'-p';};
fem.bnd.r{1,2} = 0;
fem.bnd.r{1,3} = 0;
fem.bnd.r{1,4} = {'-c+c0_c_chdi';'-p';};

fem.bnd.q = {0, 0, {0,0;'-diff(N_p_chdi,c)','-diff(N_p_chdi,p)';}, 0};

fem.bnd.ind   = [1,2,2,2,2,3,3,3,3,4];
fem.bnd.var   = {'ndflux_c_chdi','nx_chdi*dflux_c_x_chdi+ny_chdi*dflux_c_y_chdi+nz_chdi*dflux_c_z_chdi','ndflux_p_chdi','nx_chdi*dflux_p_x_chdi+ny_chdi*dflux_p_y_chdi+nz_chdi*dflux_p_z_chdi','d_chdi',1,'nx_chdi','nx','ny_chdi','ny','nz_chdi','nz','N_c_chdi',0,'kc_c_chdi',0,'cb_c_chdi',0,'c0_c_chdi',{'c0L',0,0,'c0R';},'Dbnd_c_chdi',0,'N_p_chdi',{0,0,'p0',0;},'kc_p_chdi',0,'cb_p_chdi',0,'c0_p_chdi',0,'Dbnd_p_chdi',0;};

fem.edg.shape = [1;2;];
fem.edg.ind   = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;];

fem.pnt.shape = [1;2;];
fem.pnt.ind   = [1,1,1,1,1,1,1,1,1,1,1,1,1;];

% fem.draw.p.objs = {};
% fem.draw.p.name = {};
% fem.draw.c.objs = {};
% fem.draw.c.name = {};
% fem.draw.f.objs = {};
% fem.draw.f.name = {};
% fem.draw.s.objs = {<1x1 solid3>;};
% fem.draw.s.name = {};