function [ S ] = calculMatriceS(W, D)
%% Fonction de création de la matrice S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables d'entrée :
% W = Conteneur des t objets
% D = Métrique des poids des individus
%
% Variables de sortie : 
% S = Matrice des produits scalaires d'Hilbert Schmidt
%
% Usage:
% [ S ] = calculMatriceS(W, D)
%
% Auteur: Mounir Bendali-Braham
% Date de création : Mai 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Initialisation
    % T : nombre d'objets
    T = size(W,3); 

    %% Calcul de la matrice S
    for t1 =1:T
        for t2 =1:T
            S(t1, t2) = prod_hs(W(:,:,t1), W(:,:,t2), D);
        end;
    end;


end

