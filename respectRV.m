function [ respect ] = respectRV( W, roundedRV, cible, D)
respect = true;
    for i=1:cible-1
        if(~isequal(0.1*round(distanceRV(W(:,:,cible), W(:,:,i), D)*10), roundedRV(cible, i)))
            respect = false;
            break;
        end;
    end;

end

