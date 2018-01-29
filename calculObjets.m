function [ O ] = calculObjets(X, M, norm)
    %-------------------------------------------------------------------------------
    % Definition des objets
    %-------------------------------------------------------------------------------
    nbObjets = size(X,3);
    for t = 1:nbObjets
        if ~norm
           O(:,:,t) = X(:,:,t)*M*X(:,:,t)';
        else
            O(:,:,t) = X(:,:,t)*M*X(:,:,t)';
            O(:,:,t) = O(:,:,t)./norme(O(:,:,t));
        end     
    end
end

