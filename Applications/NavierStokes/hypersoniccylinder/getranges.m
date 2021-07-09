%% params
nt = length(pde.soltime);

densities = zeros(nt,2);
energies = zeros(nt,2);
rhoprevmin = 10;
rhoprevmax = 0;
rhoeprevmax = 0;
rhoeprevmin = 1e8;
eprevmax= 0; eprevmin = 1e8;
rscale = 0.001;
vscale = 5000;
pscale = rscale * vscale^2;
%% get ranges
for ti = 1:nt
    Uout = getsolution(['dataout/out_t' num2str(pde.soltime(ti))],dmd,master.npe);
    rho = Uout(:,1,:);
    uv = Uout(:,2,:)./Uout(:,1,:);
    vv = Uout(:,3,:)./Uout(:,1,:);
    rhoE = Uout(:,4,:);

    rhoe = rhoE - rho.*(uv.^2 + vv.^2)/2;
    e = rhoE./rho - (uv.^2 + vv.^2)/2;
%     rhoe = rhoE - 0.5*(Uout(:,2,:).*uv + Uout(:,3,:).*vv);
    densities(ti,1) = min(rho,[],'all');
    densities(ti,2) = max(rho,[],'all');
    
    rhoprevmin = min(rhoprevmin, min(rho*rscale,[],'all'));
    rhoprevmax = max(rhoprevmax, max(rho*rscale, [],'all'));
   
    
    energies(ti,1) = min(rhoe,[],'all');
    energies(ti,2) = max(rhoe,[],'all');
    
    rhoeprevmin = min(rhoeprevmin, min(rhoe*pscale, [],'all'));
    rhoeprevmax = max(rhoeprevmax, max(rhoe*pscale, [],'all'));
    eprevmin = min(eprevmin, min(e*vscale^2, [],'all'));
    eprevmax = max(eprevmax, max(e*vscale^2, [],'all'));
    
    p = (gam-1)*rhoe;
    T = Minf^2 * gam * p(:)./rho(:) * Tref;
    rhot = rscale*rho(:); et = (rhoe(:)./rho(:))*vscale^2;
    rhoet = rhoe(:)*pscale;    
    disp(sum(rhot(:)<0.00085)./length(rhot(:)))

    
    hold on
    subplot(1,2,2)
    scatter(rhot(1:end),et(1:end),'.')
    hold on
    subplot(1,2,1)
     scatter(rhot(1:5:end),rhoet(1:5:end),'.')
%     figure(4)
%     scatter(rhot(1:5:end),T(1:5:end),'.')
    
    waitforbuttonpress;
end
subplot(1,2,1)
grid on
subplot(1,2,2)
grid on
% figure(1);
% plot(1:nt,densities(:,1), 1:nt,densities(:,2))
% figure(2);
% plot(1:nt,energies(:,1), 1:nt,energies(:,2))
