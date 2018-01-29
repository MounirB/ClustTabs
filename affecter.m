function [ index_affect ] = affecter(W, centres_Wk, D, s, s_centres)
%% Fonction d'affectation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables d'entrée :
% W = Conteneur des t objets
% centres_Wk = Conteneur des k centres Wk
% D = Métrique des poids des variables
% s = Vecteur des normes des objets W
% s_centres = Vecteur des normes des centres centres_Wk
%
% Variables de sortie : 
% index_affect = Index d'affectation des objets
%
% Usage:
% [ index_affect ] = affecter(W, centres_Wk, D, s, s_centres)
%
% Auteur: Mounir Bendali-Braham
% Date de création : Avril 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nbObjets = size(W, 3);
nbCentres = size(centres_Wk, 3);
index_affect = ones(1, nbObjets);

%% Calcul de la matrice des produits scalaires d'Hilbert Schmidt S
    S = ones(nbObjets, nbCentres);
    for t = 1:nbObjets
        for k = 1:nbCentres
            S(t,k) = prod_hs(W(:,:,t), centres_Wk(:,:,k), D);
        end;
    end;
    
%% Calcul des coefficients de corrélation RV
   distRV = ones(nbObjets, nbCentres);
   for t = 1:nbObjets
        for k = 1:nbCentres
            distRV(t,k) = S(t,k)/(sqrt(s(t))*sqrt(s_centres(k)));
        end
   end
   
%% Index d'affectation
    for t = 1:nbObjets
        temp = find(distRV(t,:)==max(distRV(t,:)));
        index_affect(t) = temp(1);
    end;
end

