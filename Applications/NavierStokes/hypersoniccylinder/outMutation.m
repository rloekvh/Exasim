mutationdata_out

currdat = outhighres;

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

gam = 1.4*ones(size(rho));
P = (gam-1).*rho.*(e);
% T = Minf^2*gam.*(P/pscale)./(rho/rscale)*Tref;
T = P./(287*rho);
% T = Tref/Tinf * T;
a = sqrt(gam.*P./rho);
muRef = 1/Re;
mu = getViscosity(muRef,Tref,T,1);
kappa = mu.*gam/(Pr)*Tref/Tinf;
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
surf(rho,e,gam)
title("\gamma")
xlabel("\rho")
ylabel("e");
subplot(2,3,6);
surf(rho,e,a)
title("a")
xlabel("\rho")
ylabel("e");
