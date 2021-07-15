% mutationdata_out
% 
% % currdat = outhighres;
% currdat = quadtest;
% rho = reshape(currdat(:,1),sqrt(size(currdat,1)),sqrt(size(currdat,1)));
% e = reshape(currdat(:,2),sqrt(size(currdat,1)),sqrt(size(currdat,1)));
% T = reshape(currdat(:,3),sqrt(size(currdat,1)),sqrt(size(currdat,1)));
% P = reshape(currdat(:,4),sqrt(size(currdat,1)),sqrt(size(currdat,1)));
% mu = reshape(currdat(:,5),sqrt(size(currdat,1)),sqrt(size(currdat,1)));
% kappa = reshape(currdat(:,6),sqrt(size(currdat,1)),sqrt(size(currdat,1)));
% gamma = reshape(currdat(:,7),sqrt(size(currdat,1)),sqrt(size(currdat,1)));
% a = reshape(currdat(:,8),sqrt(size(currdat,1)),sqrt(size(currdat,1)));
% figure(2)
% subplot(2,3,1);
% surf(rho,e,T)
% title("T")
% xlabel("\rho")
% ylabel("e");
% subplot(2,3,2);
% surf(rho,e,P)
% title("P")
% xlabel("\rho")
% ylabel("e");
% subplot(2,3,3);
% surf(rho,e,mu)
% title("\mu")
% xlabel("\rho")
% ylabel("e");
% subplot(2,3,4);
% surf(rho,e,kappa)
% title("\kappa")
% xlabel("\rho")
% ylabel("e");
% subplot(2,3,5);
% surf(rho,e,gamma)
% title("\gamma")
% xlabel("\rho")
% ylabel("e");
% subplot(2,3,6);
% surf(rho,e,a)
% title("a")
% xlabel("\rho")
% ylabel("e");
% 
% 
% %% ideal gas version
% rscale = 0.001;
% vscale = 5000;
% pscale = rscale*vscale^2;
% rho = reshape(currdat(:,1),sqrt(size(currdat,1)),sqrt(size(currdat,1)));
% e = reshape(currdat(:,2),sqrt(size(currdat,1)),sqrt(size(currdat,1)));
% 
% gam = 1.4*ones(size(rho));
% P = (gam-1).*rho.*(e);
% % T = Minf^2*gam.*(P/pscale)./(rho/rscale)*Tref;
% T = P./(287*rho);
% % T = Tref/Tinf * T;
% a = sqrt(gam.*P./rho);
% muRef = 1/Re;
% mu = getViscosity(muRef,Tref,T,1);
% kappa = mu.*gam/(Pr)*Tref/Tinf;
% figure(3)
% subplot(2,3,1);
% surf(rho,e,T)
% title("T")
% xlabel("\rho")
% ylabel("e");
% subplot(2,3,2);
% surf(rho,e,P)
% title("P")
% xlabel("\rho")
% ylabel("e");
% subplot(2,3,3);
% surf(rho,e,mu)
% title("\mu")
% xlabel("\rho")
% ylabel("e");
% 
% subplot(2,3,4);
% surf(rho,e,kappa)
% title("\kappa")
% xlabel("\rho")
% ylabel("e");
% subplot(2,3,5);
% surf(rho,e,gam)
% title("\gamma")
% xlabel("\rho")
% ylabel("e");
% subplot(2,3,6);
% surf(rho,e,a)
% title("a")
% xlabel("\rho")
% ylabel("e");

mutationdata_out

currdat = outhighres3;

rho = reshape(currdat(:,1),sqrt(size(currdat,1)),sqrt(size(currdat,1)));
e = reshape(currdat(:,2),sqrt(size(currdat,1)),sqrt(size(currdat,1)));
T = reshape(currdat(:,3),sqrt(size(currdat,1)),sqrt(size(currdat,1)));
P = reshape(currdat(:,4),sqrt(size(currdat,1)),sqrt(size(currdat,1)));
mu = reshape(currdat(:,5),sqrt(size(currdat,1)),sqrt(size(currdat,1)));
kappa = reshape(currdat(:,6),sqrt(size(currdat,1)),sqrt(size(currdat,1)));
gamma = reshape(currdat(:,7),sqrt(size(currdat,1)),sqrt(size(currdat,1)));
a = reshape(currdat(:,8),sqrt(size(currdat,1)),sqrt(size(currdat,1)));
figure(2)
subplot(2,3,1);
surf(rho,e,T)
title("T")
xlabel("\rho")
ylabel("e");
subplot(2,3,2);
surf(rho,e,P)
title("P")
xlabel("\rho")
ylabel("e");
subplot(2,3,3);
surf(rho,e,mu)
title("\mu")
xlabel("\rho")
ylabel("e");
subplot(2,3,4);
surf(rho,e,kappa)
title("\kappa")
xlabel("\rho")
ylabel("e");
subplot(2,3,5);
surf(rho,e,gamma)
title("\gamma")
xlabel("\rho")
ylabel("e");
subplot(2,3,6);
surf(rho,e,a)
title("a")
xlabel("\rho")
ylabel("e");


%% ideal gas version
rscale = 0.001;
vscale = 5000;
pscale = rscale*vscale^2;
rho = reshape(currdat(:,1),sqrt(size(currdat,1)),sqrt(size(currdat,1)));
e = reshape(currdat(:,2),sqrt(size(currdat,1)),sqrt(size(currdat,1)));

% OKAY SO WHAT DOES "NONDIMENSIONALIZATION" LOOK LIKE IN PRACTICE
% Code variables || Mutation Uses
% rho                  rho * rscale^2
% e                    e * vscale^2 
% P                    Pphys = P * pscale
% T                    Tphys = T * Tref/Tinf
% a                    a * vscale^2 (I'd suspect?)
% mu                   same mu ? 
%                        why does the code seem to not worry about
%                        dimensionality of mu??? The way this is set up I
%                        find confusing...
% kappa                kappaPhys = kappa * Tref/Tinf 
%                       (seems like it starts dimensionally correct but
%                       in the code it includes the temperature scaling)

gamma = 1.4*ones(size(rho));
P = (gamma-1).*rho.*(e);
T = Minf^2*gam.*(P/pscale)./(rho/rscale)*Tref;
% T = P./(287*rho);
% T = Tref/Tinf * T;
a = sqrt(gamma.*P./rho);
muRef = 1/Re;
mu = getViscosity(muRef,Tref,T,1);
codeKappa = mu*gamma/(Pr);

kappa = mu.*gamma/(Pr)*Tref/Tinf;
figure(3)
subplot(2,3,1);
surf(rho,e,T)
title("T")
xlabel("\rho")
ylabel("e");
subplot(2,3,2);
surf(rho,e,P)
title("P")
xlabel("\rho")
ylabel("e");
subplot(2,3,3);
surf(rho,e,mu)
title("\mu")
xlabel("\rho")
ylabel("e");

subplot(2,3,4);
surf(rho,e,kappa)
title("\kappa")
xlabel("\rho")
ylabel("e");
subplot(2,3,5);
surf(rho,e,gamma)
title("\gamma")
xlabel("\rho")
ylabel("e");
subplot(2,3,6);
surf(rho,e,a)
title("a")
xlabel("\rho")
ylabel("e");


%% nondim version
gammanum = gamma;
Pnum = P/pscale;
rhonum = rho;
enum = e/(vscale^2);
% T = Minf^2*gam.*(P/pscale)./(rho/rscale)*Tref;
% T = P./(287*rho);
Tnum = T*Tinf/Tref;
% T = Tref/Tinf * T;
% a = sqrt(gamma.*P./rho);
anum = a/vscale^2;
% muRef = 1/Re;
% mu = getViscosity(muRef,Tref,T,1);
munum = mu;
codeKappa = mu*gamma/(Pr);
kappanum = kappa*Tinf/Tref;
figure(3)
subplot(2,3,1);
surf(rhonum,enum,Tnum)
title("T")
xlabel("\rho")
ylabel("e");
subplot(2,3,2);
surf(rhonum,enum,Pnum)
title("P")
xlabel("\rho")
ylabel("e");
subplot(2,3,3);
surf(rhonum,enum,munum)
title("\mu")
xlabel("\rho")
ylabel("e");

subplot(2,3,4);
surf(rhonum,enum,kappanum)
title("\kappa")
xlabel("\rho")
ylabel("e");
subplot(2,3,5);
surf(rhonum,enum,gammanum)
title("\gamma")
xlabel("\rho")
ylabel("e");
subplot(2,3,6);
surf(rhonum,enum,anum)
title("a")
xlabel("\rho")
ylabel("e");

