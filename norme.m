function [An]= norme(A,D)
    %--------------------------------
    % Definition de norme
    %--------------------------------
    if nargin < 2
        n= size(A,2);
        D =1/n * eye(n);
    end

    An= sqrt(prod_hs(A,A,D));
end