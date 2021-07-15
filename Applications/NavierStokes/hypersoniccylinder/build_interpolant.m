function [A1d,T_alpha,P_alpha,mu_alpha, kappa_alpha,gamma_alpha,a_alpha] = build_interpolant(p,a_rho, b_rho,a_e, b_e,currdat,param)
% returns structures needed to contstruct high order interpolants of data
% OUT:
%   - A1d: ((p+1) x (p+1))
%     for a 1d basis function psi_k = a_k^p x^p + ... + a_k^0 x^0, each row of
%     A1d is a vector [a_k^p, a_k^{p-1}, ..., a_k^0] so that psi_k(x) can
%     be evaluated as A1d * (x).^(p:-1:0)
%   - X_alpha: (p+1)^2 
%           for an output X, the coefficients such that 
%           X(r,e) \approx \sum_{k=0}^p X_alpha psi_k(r,e)
%   

% pde.physicsparam = [gam Re Pr Minf rinf ruinf rvinf rEinf Tinf Tref Twall avb avk avs pde.porder sb0 sb1 rscale vscale pscale];
                %    1  2  3   4    5     6     7      8    9   10    11   12  13  14     15     16  17   18     19      20     
Re = param(3);
Tinf = param(9);
Tref = param(10);
rscale = param(18);
vscale = param(19);
pscale = param(20);

muRef = 1/Re;
                
L_rho = b_rho - a_rho;
L_e = b_e - a_e;
% n_alpha = (p+1)^2;
[xg1d,wg1d] = gaussquad(2*p, 1, 1);
[xg2d,wg2d] = gaussquad(2*p, 2, 1);
intPsiPsi1d = wg1d'*koornwinder(xg1d,p).^2;
% intPsiPsi2d = intPsiPsi1d.^2*L_rho*L_e;


% mutoutfile = dlmread('/Users/loekvh/Desktop/mit/research/Mutationpp/examples/c++/equilibrium_air/example.txt');
% A = load('example.mat','mutoutfile');
% currdat = mutoutfile;
Tphys = currdat(:,3);
Pphys = currdat(:,4);
muphys = currdat(:,5);
kappa = currdat(:,6);
gamma = currdat(:,7);
a = currdat(:,8);
% Tphys = reshape(currdat(:,3),sqrt(size(currdat,1)),sqrt(size(currdat,1)));
T = Tphys*Tinf/Tref;
P = Pphys/pscale;
mu = muphys/(muRef*50);
basisxquad = tensorproduct(xg2d,p);
calcCoeffs = @(X) (X.*wg2d)'*basisxquad;
T_alpha = calcCoeffs(T);
P_alpha = calcCoeffs(P);
gamma_alpha = calcCoeffs(gamma);
mu_alpha = calcCoeffs(mu);
kappa_alpha = calcCoeffs(kappa);
a_alpha = calcCoeffs(a);

A1d = koornwinderCoeffs1d(p);
end

