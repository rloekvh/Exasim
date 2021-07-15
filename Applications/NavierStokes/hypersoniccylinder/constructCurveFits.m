%% Curve fits of thermo and transport properties with orthogonal polynomials

a_rho = 0.0005;
b_rho = 0.2;
a_e = 0.005;
b_e = 0.5;

L_rho = b_rho - a_rho;
L_e = b_e - a_e;

%% step 1: calculate < psi_a, psi_a > 
porder = 12;
n_alpha = (porder+1)^2;

% get quadrature nodes
[xg1d,wg1d] = gaussquad(2*porder, 1, 1);
[xg2d,wg2d] = gaussquad(2*porder, 2, 1);
% evaluate at quad nodes

%rho integral
intPsiPsi1d = wg1d'*koornwinder(xg1d,porder).^2;
%e integral
intPsiPsi2d = intPsiPsi1d.^2*L_rho*L_e;
%what poylnomials are included in each index

%% Calculate <T, psi_a>
%okay so map quad points to interval first...
savefile='/Users/loekvh/Desktop/mit/research/Mutationpp/examples/c++/equilibrium_air/matlabinput.mat';
rhoquadpts = a_rho + xg1d*L_rho; %DO NOT NONDIMENZIONALIZE
equadpts = (a_e + xg1d*L_e)*vscale^2;       %DO DIMENSIONALIZE
[xx,yy] = ndgrid(rhoquadpts,equadpts);
figure(1); scatter(xx(:),yy(:));
format long g
% disp([rhoquadpts, equadpts])
A = [rhoquadpts, equadpts];
save(savefile,'A')

%% RUN MUTATION FIRST
% could get ugly here...hopefully i can evaluate the state at appropriate
% quadrature nodes
% mutationdata_outquad;
% currdat = quadtest7;
mutoutfile = dlmread('/Users/loekvh/Desktop/mit/research/Mutationpp/examples/c++/equilibrium_air/example.txt');
% A = load('example.mat','mutoutfile');
currdat = mutoutfile;
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
basisxquad = tensorproduct(xg2d,porder);
calcCoeffs = @(X) (X.*wg2d)'*basisxquad;
T_alpha = calcCoeffs(T);
P_alpha = calcCoeffs(P);
gamma_alpha = calcCoeffs(gamma);
mu_alpha = calcCoeffs(mu);
kappa_alpha = calcCoeffs(kappa);
a_alpha = calcCoeffs(a);
% c_alpha = zeros(n_alpha,1);
% for i = 1:porder+1
%     for j = 1:porder+1
%     % integrate <T, psi_a>
%         c_alpha((i-1)*nalpha+j) = 
%     end
% end
%and then just do a division

%% visualize evaluation
rhotst = linspace(a_rho,b_rho)';
etst = linspace(a_e,b_e)';
ntst = length(rhotst);

xi = (rhotst-a_rho)/(b_rho-a_rho);
eta = (etst-a_e)/(b_e - a_e);

[xx,yy] = ndgrid(xi,eta);
basisxtest = tensorproduct([xx(:),yy(:)],porder);
evalExpansion = @(X_alpha) sum(basisxtest*diag(X_alpha),2);

% Tphys = currdat(:,3);
% Pphys = currdat(:,4);
% mu = currdat(:,5);
% kappa = currdat(:,6);
% gamma = currdat(:,7);
% a = currdat(:,8);
Teval = evalExpansion(T_alpha);
Teval = Teval*Tref/Tinf;
Teval = reshape(Teval,[ntst, ntst]);

Peval = evalExpansion(P_alpha);
Peval = Peval*pscale;
Peval = reshape(Peval,[ntst,ntst]);

mueval = evalExpansion(mu_alpha);
mueval = mueval*(muRef*50); %no scaling? 
mueval = reshape(mueval,[ntst,ntst]);

kappaeval = evalExpansion(kappa_alpha);
% kappaeval = kappaeval;
kappaeval = reshape(kappaeval,[ntst,ntst]);

aeval = evalExpansion(a_alpha);
% kappaeval = kappaeval;
aeval = reshape(aeval,[ntst,ntst]);

gammaeval = evalExpansion(gamma_alpha);
gammaeval = reshape(gammaeval,[ntst,ntst]);


% surf(rhotst,etst*vscale^2,Teval')
% surf(rhotst,etst*vscale^2,Teval')

rhonum = rhotst;
enum = etst*vscale^2;

figure(4)
subplot(2,3,1);
surf(rhonum,enum,Teval','EdgeColor','none')
title("T")
xlabel("\rho")
ylabel("e");
subplot(2,3,2);
surf(rhonum,enum,Peval','EdgeColor','none')
title("P")
xlabel("\rho")
ylabel("e");
subplot(2,3,3);
surf(rhonum,enum,mueval','EdgeColor','none')
title("\mu")
xlabel("\rho")
ylabel("e");

subplot(2,3,4);
surf(rhonum,enum,kappaeval','EdgeColor','none')
title("\kappa")
xlabel("\rho")
ylabel("e");
subplot(2,3,5);
surf(rhonum,enum,gammaeval','EdgeColor','none')
title("\gamma")
xlabel("\rho")
ylabel("e");
subplot(2,3,6);
surf(rhonum,enum,aeval','EdgeColor','none')
title("a")
xlabel("\rho")
ylabel("e");



%% check errors
mutationdata_out;
currdat = outhighres3;
mutationdata_in; 
rhotst = A(:,1); %A here is input for high res run
etst = A(:,2)/(vscale^2); %NEED TO RESCALE!
ntst = length(rhotst);

xi = (rhotst-a_rho)/(b_rho-a_rho);
eta = (etst-a_e)/(b_e - a_e);


[xx,yy] = ndgrid(xi,eta);
basisxtest = tensorproduct([xx(:),yy(:)],porder);
evalExpansion = @(X_alpha) sum(basisxtest*diag(X_alpha),2);

Tmut = reshape(currdat(:,3),[ntst ntst])';
Pmut = reshape(currdat(:,4),[ntst ntst])';
mumut = reshape(currdat(:,5), [ntst ntst])';
kappamut = reshape(currdat(:,6), [ntst ntst])';
gammamut = reshape(currdat(:,7), [ntst ntst])';
amut = reshape(currdat(:,8), [ntst ntst])';

Teval = evalExpansion(T_alpha);
Teval = Teval*Tref/Tinf;
Teval = reshape(Teval,[ntst, ntst])';

Peval = evalExpansion(P_alpha);
Peval = Peval*pscale;
Peval = reshape(Peval,[ntst,ntst])';

mueval = evalExpansion(mu_alpha);
mueval = mueval*(muRef*50); %no scaling? 
mueval = reshape(mueval,[ntst,ntst])';

kappaeval = evalExpansion(kappa_alpha);
% kappaeval = kappaeval;
kappaeval = reshape(kappaeval,[ntst,ntst])';

aeval = evalExpansion(a_alpha);
% kappaeval = kappaeval;
aeval = reshape(aeval,[ntst,ntst])';

gammaeval = evalExpansion(gamma_alpha);
gammaeval = reshape(gammaeval,[ntst,ntst])';


% surf(rhotst,etst*vscale^2,Teval')
% surf(rhotst,etst*vscale^2,Teval')

rhonum = rhotst;
enum = etst*vscale^2;

avis = a_rho;
bvis = 0.01;

figure(4)
subplot(2,3,1);
surf(rhonum,enum,abs((Teval-Tmut))./Tmut,'EdgeColor','none')
title("T")
xlabel("\rho")
ylabel("e");
% xlim([avis, bvis])

subplot(2,3,2);
surf(rhonum,enum,abs((Peval-Pmut))./Pmut,'EdgeColor','none')
title("P")
xlabel("\rho")
ylabel("e");
% xlim([avis, bvis])

subplot(2,3,3);
surf(rhonum,enum,abs((mueval-mumut))./mumut,'EdgeColor','none')
title("\mu")
xlabel("\rho")
ylabel("e");
% xlim([avis, bvis])

subplot(2,3,4);
surf(rhonum,enum,abs((kappaeval-kappamut))./kappamut,'EdgeColor','none')
title("\kappa")
xlabel("\rho")
ylabel("e");
% xlim([avis, bvis])

subplot(2,3,5);
surf(rhonum,enum,abs((gammaeval-gammamut))./gammamut,'EdgeColor','none')
title("\gamma")
xlabel("\rho")
ylabel("e");
% xlim([avis, bvis])

subplot(2,3,6);
surf(rhonum,enum,abs((aeval-amut))./amut,'EdgeColor','none')
title("a")
xlabel("\rho")
ylabel("e");
% xlim([avis, bvis])
