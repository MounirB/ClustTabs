function [r] = prod_hs(A,B,D)
    % Definition de produit scalaire d'Hilbert Schmidt
    r = trace(A*D*B*D);
end