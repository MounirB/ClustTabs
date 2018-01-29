function [ centres_Wk ] = minimiser(W, S, D, Delta, nbCentres, index_classif, norm)
%% Fonction de minimisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables d'entr�e :
% W = Conteneur des t objets
% S = Matrice des produits scalaires d'Hilbert Schmidt S
% D = M�trique des poids des individus
% Delta = Matrice des poids des �tudes
% index_classif = R�sultat de l'affectation
% nbCentres = Nombre de classes
% norm = Recours ou non � la normalisation des objets
%
% Variables de sortie : 
% centres_Wk = Nouveaux centro�ds de chaque classe
%
% Usage:
% [ centres_Wk ] = minimiser(W, index_classif)
%
% Auteur: Mounir Bendali-Braham
% Date de cr�ation : Avril 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pr�traitements
    [nbIndividus, nbIndividus, nbEtudes] = size(W);
    centres_Wk = ones(nbIndividus, nbIndividus, nbCentres);
    
%% Minimisation
    for k = 1:nbCentres
        % Image euclidienne des objets
        k_index_classif = find(index_classif==k);
        if(~isempty(k_index_classif))
            localDelta = Delta(k_index_classif, k_index_classif)*(nbEtudes/length(k_index_classif));
            localS = S(k_index_classif, k_index_classif);
            SDelta = sqrt(localDelta)'*localS*sqrt(localDelta);
            [~, VaP, VeP] = ACP(SDelta);
            [centres_Wk(:,:,k), ~] = compromis(W(:,:,k_index_classif), localS, ...
            localDelta, VaP, VeP, norm);
        end;
    end;
end

