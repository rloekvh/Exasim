function f = getfluxav2d_wEoS(udg,qdg,vdg,wdg,param,eta)

% pde.physicsparam = [gam Re Pr Minf rinf ruinf rvinf rEinf Tinf Tref Twall avb avk avs pde.porder sb0 sb1 rscale vscale pscale];
%                      1  2  3   4    5     6     7      8    9   10    11   12  13  14     15     16  17   18     19      20     
% gam = param(1);
% gam1 = gam - 1.0;
% Re = param(2);
% Pr = param(3);
% Minf = param(4);
% Tref = param(10);
% muRef = 1/Re;
% Tinf = 1/(gam*gam1*Minf^2);
c23 = 2.0/3.0;

% regularization parameters
alpha = 1.0e3;
rmin = 1.0e-2;
pmin = 1.0e-3;

avb = vdg(1); % Bulk    viscosity
avr = 0;      % density    viscosity
avk = param(13)*vdg(1); % Thermal viscosity
avs = param(14)*vdg(1); % Shear   viscosity

r = udg(1);
ru = udg(2);
rv = udg(3);
rE = udg(4);
dr_dx = qdg(1);
rux = qdg(2);
rvx = qdg(3);
drE_dx = qdg(4);
dr_dy = qdg(5);
ruy = qdg(6);
rvy = qdg(7);
drE_dy = qdg(8);
% wdg = [p, mu, kappa, gamma, a, dT_dri, dT_d(re)];

p = wdg(1);
mu = wdg(2);
kappa = wdg(3);
% gam = wdg(4); %only used inside equation of state for ideal gas.
% a = wdg(5); %speed of sound only used in boundary conditions
dT_dri = wdg(6);
dT_dre = wdg(7);

% gam1 = gam - 1.0;

r = rmin + lmax(r-rmin,alpha);
% Density sensor
dr = atan(alpha*(r - rmin))/pi + (alpha*(r - rmin))/(pi*(alpha^2*(r - rmin)^2 + 1)) + 1/2;
%dr=1;
dr_dx = dr_dx*dr;
dr_dy = dr_dy*dr;
r1 = 1/r;
uv = ru*r1;
vv = rv*r1;
E = rE*r1;
q = 0.5*(uv*uv+vv*vv);
% p = gam1*(rE-r*q);
% Regularization of pressure p (cannot be smaller than pmin)
p = pmin + lmax(p-pmin,alpha);
% Pressure sensor
dp = atan(alpha*(p - pmin))/pi + (alpha*(p - pmin))/(pi*(alpha^2*(p - pmin)^2 + 1)) + 1/2;
%dp=1;
% Total enthalpy
h = E+p*r1;
% Inviscid fluxes
fi = [ru, ru*uv+p, rv*uv, ru*h, ...
        rv, ru*vv, rv*vv+p, rv*h];
ux = (rux - dr_dx*uv)*r1;
vx = (rvx - dr_dx*vv)*r1;
dq_dx = uv*ux + vv*vx;
% px = gam1*(drE_dx - dr_dx*q - r*dq_dx);
% px = px*dp;
% Tx = 1/gam1*(px*r - p*dr_dx)*r1^2;
uy = (ruy - dr_dy*uv)*r1;
vy = (rvy - dr_dy*vv)*r1;
dq_dy = uv*uy + vv*vy;

% re = rE - r q/2
dre_dr = -q/2;
dre_dq = -r/2;
dre_drE = 1;
dre_dx = dre_dr * dr_dx + dre_dq * dq_dx + dre_drE * drE_dx;
dre_dy = dre_dr * dr_dy + dre_dq * dq_dy + dre_drE * drE_dy;

Tx = dT_dri * dr_dx + dT_dre * dre_dx;
Ty = dT_dri * dr_dy + dT_dre * dre_dy;

% py = gam1*(drE_dy - dr_dy*q - r*dq_dy);
% py = py*dp;
% Ty = 1/gam1*(py*r - p*dr_dy)*r1^2;
% Adding Artificial viscosities
% T = p/(gam1*r);
% 
% Tphys = Tref/Tinf * T;
% mu = getViscosity(muRef,Tref,Tphys,1);
mu = mu + avs;
fc = kappa;
% Viscous fluxes with artificial viscosities
txx = (mu)*c23*(2*ux - vy) + (avb)*(ux+vy);
txy = (mu)*(uy + vx);
tyy = (mu)*c23*(2*vy - ux) + (avb)*(ux+vy);
fv = [avr*dr_dx, txx, txy, uv*txx + vv*txy + (fc+avk)*Tx, ...
      avr*dr_dy, txy, tyy, uv*txy + vv*tyy + (fc+avk)*Ty];
f = fi+fv;

