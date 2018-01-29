function [ RV ] = calculMatriceRV(W, D)
%% Fonction de cr�ation de la matrice RV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables d'entr�e :
% W = Conteneur des t objets
% D = M�trique des poids des individus
%
% Variables de sortie : 
% RV = Matrice des produits scalaires d'Hilbert Schmidt
%
% Usage:
% [ RV ] = calculMatriceS(W, D)
%
% Auteur: Mounir Bendali-Braham
% Date de cr�ation : Mai 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Initialisation
    % T : nombre d'objets
    T = size(W,3); 

    %% Calcul de la matrice S
    for t1 =1:T
        for t2 =1:T
            RV(t1, t2) = prod_hs(W(:,:,t1), W(:,:,t2), D)/(sqrt(prod_hs(W(:,:,t1),W(:,:,t1),D))*sqrt(prod_hs(W(:,:,t2),W(:,:,t2),D)));
        end;
    end;

end

