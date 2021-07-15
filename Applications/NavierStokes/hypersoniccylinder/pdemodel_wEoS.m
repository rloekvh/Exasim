function pde = pdemodel_wEoS
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.initu = @initu;
pde.avfield = @avfield;
pde.sourcew = @sourcew;
pde.initw = @initw;
end

function m = mass(u, q, w, v, x, t, mu, eta)
m = sym([1.0; 1.0; 1.0; 1.0]); 
end

function f = flux(u, q, w, v, x, t, mu, eta)
    f = getfluxav2d_wEoS(u,q,v,w,mu,eta);
    f = reshape(f,[4,2]);        
end

function f = avfield(u, q, w, v, x, t, mu, eta)
    f = getavfield2d_wEoS(u,q,v,w,mu);
end

function s = source(u, q, w, v, x, t, mu, eta)
    s = [sym(0.0); sym(0.0); sym(0.0); sym(0.0)];
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    f = flux(uhat, q, w, v, x, t, mu, eta);
    fi = f(:,1)*n(1) + f(:,2)*n(2) + tau*(u-uhat); % numerical flux at freestream boundary
    
    % adiabatic wall
    faw = fi;
    faw(1) = 0.0;   % zero velocity 
    faw(end) = 0.0; % adiabatic wall -> zero heat flux
    
    % Flux Thermal Wall
    ftw = fi;
    ftw(1) = 0.0;
    
    % freestream, adiabatic wall, isothermal wall, adiabatic slip wall, supersonic inflow, supersonic outflow
    fb = [fi faw ftw faw fi fi]; 
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)

    % gam = mu(1);
    % gam1 = gam - 1.0;
    % p = wdg(1);
    % mu = wdg(2);
    % kappa = wdg(3);
    gam = w(4); %only used inside equation of state for ideal gas.
    a = w(5); %speed of sound only used in boundary conditions
    % dT_dri = wdg(6);
    % dT_dre = wdg(7);
    gam1 = gam - 1.0;

    r = u(1);
    ru = u(2);
    rv = u(3);
    rE = u(4);
    nx = n(1);
    ny = n(2);
    
    rinv = 1/r;
    uv = ru*rinv;
    vv = rv*rinv;
    E = rE*rinv;
    p = gam1*(rE-r*0.5*(uv*uv+vv*vv));
    h = E+p*rinv;

    % a = sqrt(gam*p*rinv);
    
    run = ru*nx + rv*ny;
    rut = -ru*ny + rv*nx;
    un = run/r;
    ut = rut/r;
    
    K = [ 1 , 1 , 0 , 1 ;...
          un-a , un , 0 , un+a ;...
          ut , ut , 1 , ut ;...
          h - un*a , (1/2)*(un^2 + ut^2) , ut , h+un*a ];
    Kinv = (gam1/(2*a^2))*[ h + (a/gam1)*(un-a) , -(un+a/gam1) , -ut , 1 ;...
                            -2*h + (4/gam1)*a^2 , 2*un , 2*ut , -2 ;...
                            -2*(ut*a^2)/gam1 , 0 , 2*(a^2)/gam1 , 0 ;...
                            h - a*(un+a)/gam1 , -un+a/gam1 , -ut , 1 ];
    T = [ 1 , 0 , 0 , 0;...
          0 , nx , ny , 0;...
          0 , -ny , nx , 0;...
          0 , 0 , 0 , 1];
    Tinv = [ 1 , 0 , 0 , 0;...
             0 , nx ,-ny , 0;...
             0 , ny , nx , 0;...
             0 , 0 , 0 , 1];
    Lambda = [ tanh(1e3*(un-a)) , 0 , 0 , 0 ;...
                     0 , tanh(1e3*(un)) , 0 , 0 ;...
                     0 , 0 , tanh(1e3*(un)) , 0 ;...
                     0 , 0 , 0 , tanh(1e3*(un+a)) ];
    E = (K * Lambda * Kinv);
    An = (Tinv * E * T);
    
    % freestream boundary condition
    uinf = sym(mu(5:8)); % freestream flow
    uinf = uinf(:);
    u = u(:);          % state variables 
    
    % freestream
    ui = 0.5*((u+uinf) + An*(u-uinf));  % Riemann solution
    
    % adiabatic wall
    uaw = u;
    uaw(2:3) = 0; % zero velocity

    % Isothermal Wall
    Tinf = mu(9);
    Tref = mu(10);
    Twall = mu(11); %ouch... do we need to edit Tinf here? 
    TisoW = Twall/Tref * Tinf;
    utw = u(:);
    utw(2:3) = 0;
    utw(4) = u(1)*TisoW;
    
    % Slip wall
    usw = u;
    usw(2) = u(2) - nx * (u(2)*nx + u(3)*ny);
    usw(3) = u(3) - ny * (u(2)*nx + u(3)*ny);
    
    % freestream, adiabatic wall, isothermal wall, adiabatic slip wall, supersonic inflow, supersonic outflow
    ub = [ui uaw utw usw uinf u]; 
end

function u0 = initu(x, mu, eta)
    u0 = sym(mu(5:8)); % freestream flow   
end

function w0 = initw(x, mu, eta)
    w0 = sourcew(sym(mu(5:8)), [],[],[],x,0,mu,eta);
end

function sw = sourcew(u, q, w, v, x, t, mu, eta)
    sw = fluidmodel(u, t, mu, eta, 0);
end


function fluxquantities = fluidmodel(u,t,mu,eta,modelflag)
    % Returns nondimensional quantities needed for flux and eventually source terms
    % fluxquantities: given conservative variables, return thermo and
    % transport quantities
    
%     modelFlag 0 = ideal gas
%               1 = calorically imperfect gas
%     wdg = [p, mu, kappa, gamma, a, dT_dri, dT_d(re)];
    if modelflag == 0
%         fluxquantities = zeros(1,7);
        Re = mu(2);
        Pr = mu(3);
        Minf = mu(4);
        Tref = mu(10);
        % should Tinf be from here?
        
        r  = u(1);
        ru = u(2);
        rv = u(3);
        rE = u(4);
    
        rinv = 1/r;
        u = ru*rinv;
        v = rv*rinv;
        E = rE*rinv;
        
        gamma = 1.4;
        gammaminus1 = gamma-1;
        Ekinetic = 0.5*r*(u*u + v*v);
        
        % Ideal gas pressure
        re = rE - Ekinetic;
        
        p = gammaminus1 * (rE - Ekinetic);
        
        % Sutherland's law for viscosity
        muRef = 1/Re;
        Tinf = 1/(gamma*gammaminus1*Minf^2);
        T = p/(gammaminus1*r);
        Tphys = Tref/Tinf * T;
        mu = getViscosity(muRef,Tref,Tphys,1);
        
        % Prandtl number for conductivity
        kappa = mu*gamma/(Pr);
        
        % speed of sound
        a = sqrt(gamma*p*rinv);
        
        % Temperature derivatives
        rinv2 = rinv^2;
        dT_dri = -re*rinv2;
        dT_dre = rinv;
        
        fluxquantities = [p, mu, kappa, gamma, a, dT_dri, dT_dre];
    else
        error("Model not implemented");
    end

end
