function [ O ] = calculObjets_bisFromCells(Xcells)
    %-------------------------------------------------------------------------------
    % Definition des objets
    %-------------------------------------------------------------------------------
    nbObjets = size(Xcells,2);
    for t = 1:nbObjets
        O(:,:,t) = cell2mat(Xcells(t))*cell2mat(Xcells(t))';  
    end
end


