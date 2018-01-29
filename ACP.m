function [XU,VAPU, VEPU] = ACP(X)
%--------------------------------
% Calcul ACP
%--------------------------------
% Recherche des valeurs et vecteurs propres
[VEPU, VAPU] = eig(X); 
VAPU         = diag(VAPU);        
%
% Ordonnancement des valeurs et vecteurs propres
[VAPU,s] = sort(VAPU, 'descend');
%VAPU     = VAPU(s); 
VEPU     = VEPU(:,s);
VEPU=sign(VEPU(1,1)).*VEPU;
%
% Nouvelles Coordonnées (Composantes principales)
XU = VEPU * diag(sqrt(VAPU));
%XU=sign(XU(1,1)).*XU;
end

