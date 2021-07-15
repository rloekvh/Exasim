function avField = getavfield2d_wEoS(udg,qdg,vdg,wdg,param)

% wdg = [p, mu, kappa, gamma, a, dT_dri, dT_d(re)];

p = wdg(1);
% mu = wdg(2);
% kappa = wdg(3);
gam = wdg(4); %only used inside equation of state for ideal gas.
% a = wdg(5); %speed of sound only used in boundary conditions
% dT_dri = wdg(6);
% dT_dre = wdg(7);


gam1 = gam - 1.0;
% artificial viscosity
avbulk = param(12);
porder = param(15);

% regularization parameters
alpha = 1.0e3;
rmin = 1.0e-2;
Hmin = 1.0e-4;

% regularization parameters for the bulk viscosity
kb = avbulk;      % kb = 1.5 for Ma<6
sb0   = param(16);
sbmax = param(17) / sqrt(gam*gam - 1.0);
sbmin = 0.0;

% mesh size
hm = vdg(2);

% Get base variables
r = udg(1);
ru = udg(2);
rv = udg(3);
rE = udg(4);
rx = qdg(1);
rux = qdg(2);
rvx = qdg(3);
ry = qdg(5);
ruy = qdg(6);
rvy = qdg(7);

% Regularization of density 
%rreg = r0 + lmax(r-r0,alpha);
r = rmin + lmax(r-rmin,alpha);
rinv = 1./r;
uv = ru.*rinv;
vv = rv.*rinv;
E = rE.*rinv;
q = 0.5*(uv.*uv+vv.*vv);
% H = gam*E - gam1*q;                    % Enthalpy (Critical Speed of Sound ???)
H = E + p*rinv;
H = Hmin + lmax(H-Hmin,alpha);         % Regularized Enthalpy

% Critical speed of Sound %TODO: huh? 
c_star = sqrt((2.*gam1.*H) ./ (gam+1));

% Computing derivatives for the sensors
ux = (rux - rx.*uv).*rinv;
vx = (rvx - rx.*vv).*rinv;
uy = (ruy - ry.*uv).*rinv;
vy = (rvy - ry.*vv).*rinv;
div_v = - (ux + vy);
vort = - (vx - uy);
vort = sqrt(vort.*vort);

% limit  divergence and vorticity
sigm = 1e4;
div_v = limiting(div_v,-sigm,sigm,alpha,-sigm);
vort = limiting(vort,0.0,sigm,alpha,0.0);

% Dilatation Sensor sb
DucrosRatio = div_v.*div_v ./ (div_v.*div_v + vort.*vort + 1.0e-16);
% DucrosRatio = 1.0;
sb = - (hm./porder) .* (div_v./c_star) .* DucrosRatio;
sb = limiting(sb,sbmin,sbmax,alpha,sb0);
% Artificial Bulk viscosity
avb = r.*(kb.*hm./(porder)) .* sqrt(uv.*uv + vv.*vv + c_star.*c_star) .* sb;

% Assign artificial viscosities
avField(1) = avb;  %  bulk









