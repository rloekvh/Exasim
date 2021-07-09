%% grab solution
ti = randi(50);
Uout = getsolution(['dataout/out_t' num2str(pde.soltime(ti))],dmd,master.npe);
rho = Uout(:,1,:);
uv = Uout(:,2,:)./Uout(:,1,:);
vv = Uout(:,3,:)./Uout(:,1,:);
rhoE = Uout(:,4,:);
T = eulereval(Uout,'t',gam,Minf);
P = eulereval(Uout,'p',gam,Minf);
H = eulereval(Uout,'h',gam,Minf);

rhoe = rhoE - rho.*(uv.^2 + vv.^2)/2;
%% check different x
i = randi(90000);
%%
% rscale = 0.001;
rscale=.001;
vscale = 5000;
pscale = rscale*vscale^2;
disp([rho(i)*rscale, rhoe(i)*pscale])
% gamnew = gam;
%% adjust gamma?
% gamnew =   1.3926;
T = eulereval(Uout,'t',gam,Minf);
P = eulereval(Uout,'p',gam,Minf);

Ttest = (P(i)/((gam-1)*rho(i)))*Tref/Tinf;
H = eulereval(Uout,'h',1.37344,Minf);
disp("T: " + string(T(i)*Tref));
disp("Ttest: " + string(Ttest));
disp("P: " + string(P(i)*pscale));