function [F, Fxi, Feta] = evalExpansionMatrix(p,xi,eta,A1d,c_alpha)
%     A1d = koornwinderCoeffs1d(p);
    rhopows = (2*xi-1).^((p:-1:0)');
    epows = (2*eta-1).^((p:-1:0)');
    
    psirho_alpha = A1d*rhopows;
    psie_alpha = A1d*epows;
    F = 0;
    Fxi = 0;
    Feta = 0;
    for ieta = 0:p
        psie_curr = psie_alpha(ieta+1);
%         psirho_curr = psirho_alpha(i+1);
        for ixi = 0:p
            psirho_curr = psirho_alpha(ixi+1);
            disp(ixi+(p+1)*ieta+1)
%             psie_curr = psie_alpha(j+1);
            F = F + c_alpha(ixi+(p+1)*ieta+1)*psirho_curr*psie_curr;
        end
    end
end

