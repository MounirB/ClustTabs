function [ RVdist ] = distanceRV( O, centre, D )
    RVdist = prod_hs(O, centre, D) / (norme(O,D)*norme(centre,D));
end

