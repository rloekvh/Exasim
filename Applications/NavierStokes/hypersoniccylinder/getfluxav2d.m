function f = getfluxav2d(udg,qdg,vdg,wdg,param,eta)

% pde.physicsparam = [gam Re Pr Minf rinf ruinf rvinf rEinf Tinf Tref Twall avb avk avs pde.porder sb0 sb1 rscale vscale pscale];
%                      1  2  3   4    5     6     7      8    9   10    11   12  13  14     15     16  17   18     19      20     
% not all physics params should be used
% pde.externalparam = [porder_interp, A1d(:), T_alpha(:), P_alpha(:), mu_alpha(:), kappa_alpha(:), gamma_alpha(:), a_alpha(:), a_rho, b_rho, a_e, b_e];

porder_interp = 12;
nbasis = (porder_interp+1)^2; % num basis functions = (porder_interp+1)^2
A1d = reshape(eta(2:(nbasis+1)), [porder_interp+1, porder_interp+1]);
T_alpha = eta(nbasis+2:2*nbasis+1);
% P_alpha = eta(2*nbasis+2:3*nbasis+1);
% mu_alpha = eta(3*nbasis+2:4*nbasis+1);
% kappa_alpha = eta(4*nbasis+2:5*nbasis+1);
% gamma_alpha = eta(5*nbasis+2:6*nbasis+1);
% a_alpha = eta(6*nbasis+2:7*nbasis+1);
% a_rho = eta(7*nbasis+2);
% b_rho = eta(7*nbasis+3);
% a_e = eta(7*nbasis+4);
% b_e = eta(7*nbasis+5);
a_rho = 0.0005;
b_rho = 0.2;
a_e = 0.005;
b_e = 0.5;


rscale= param(18);

gam = param(1);
gam1 = gam - 1.0;
Re = param(2);
Pr = param(3);
Minf = param(4);
Tref = param(10);
muRef = 1/Re;
Tinf = 1/(gam*gam1*Minf^2);
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
rx = qdg(1);
rux = qdg(2);
rvx = qdg(3);
rEx = qdg(4);
ry = qdg(5);
ruy = qdg(6);
rvy = qdg(7);
rEy = qdg(8);
% wdg = [p, mu, kappa, gamma, a, dT_dri, dT_d(re)];

% T = wdg(1);
% p = wdg(2);
% mu
% rpows = r.^((w0:21)');
% A = rand(22*22,1);
% A = reshape(A, [22 22]);
% Regularization of rho (cannot be smaller than rmin)
r = rmin + lmax(r-rmin,alpha);
% Density sensor
dr = atan(alpha*(r - rmin))/pi + (alpha*(r - rmin))/(pi*(alpha^2*(r - rmin)^2 + 1)) + 1/2;
%dr=1;
rx = rx*dr;
ry = ry*dr;
r1 = 1/r;
uv = ru*r1;
vv = rv*r1;
E = rE*r1;
q = 0.5*(uv*uv+vv*vv);
p = gam1*(rE-r*q);
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
ux = (rux - rx*uv)*r1;
vx = (rvx - rx*vv)*r1;
qx = uv*ux + vv*vx;
px = gam1*(rEx - rx*q - r*qx);
px = px*dp;
Tx = 1/gam1*(px*r - p*rx)*r1^2;
uy = (ruy - ry*uv)*r1;
vy = (rvy - ry*vv)*r1;
qy = uv*uy + vv*vy;
py = gam1*(rEy - ry*q - r*qy);
py = py*dp;
Ty = 1/gam1*(py*r - p*ry)*r1^2;
% Adding Artificial viscosities
T = p/(gam1*r);

% NEED TO BRING INPUTS TO [0,1] DOMAIN!!!

% xi = (r*rscale-a_rho)/(b_rho-a_rho);
% eta = (E-a_e)/(b_e - a_e);
% T = evalExpansionMatrix(porder_interp,xi,eta,A1d,T_alpha);
% rhopows = (2*xi-1).^((porder_interp:-1:0)');
% epows = (2*eta-1).^((porder_interp:-1:0)');
%     
% psirho_alpha = A1d*rhopows;
% psie_alpha = A1d*epows;
% F = 0;
% for ieta = 0:porder_interp
%     psie_curr = psie_alpha(ieta+1);
% %         psirho_curr = psirho_alpha(i+1);
%     for ixi = 0:porder_interp
%         psirho_curr = psirho_alpha(ixi+1);
%         % disp(ixi+(p+1)*ieta+1)
% %             psie_curr = psie_alpha(j+1);
%         F = F + T_alpha(ixi+(porder_interp+1)*ieta+1)*psirho_curr*psie_curr;
%     end
% end

Tphys = Tref/Tinf * T;
mu = getViscosity(muRef,Tref,Tphys,1);
mu = mu + avs;
% mu = T;
fc = mu*gam/(Pr);
% Viscous fluxes with artificial viscosities
txx = (mu)*c23*(2*ux - vy) + (avb)*(ux+vy);
txy = (mu)*(uy + vx);
tyy = (mu)*c23*(2*vy - ux) + (avb)*(ux+vy);
fv = [avr*rx, txx, txy, uv*txx + vv*txy + (fc+avk)*Tx, ...
      avr*ry, txy, tyy, uv*txy + vv*tyy + (fc+avk)*Ty];
f = fi+fv;

