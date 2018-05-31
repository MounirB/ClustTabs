function [inertieRVglobale, inertiesRV_intra, nbObjets_parClasse] = inertieIntraClasse( W, D, index_classif, centres_finaux )
%% Fonction de calcul de l'inertie intra classe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables d'entrée :
% W = Conteneur des t objets
% D = Métrique des poids des individus
% index_classif = Résultat de la classification
% centres_finaux = Conteneur des centres de gravité
%
% Variables de sortie : 
% inertie = Moyenne des inerties par classe
% inerties_intra = Vecteur contenant les inerties par classe
% nbObjets_parClasse = Vecteur contenant le nombre d'objets par classe
%
% Usage:
% [inertie, inerties_intra, nbObjets_parClasse] = inertieIntraClasse(W, index_classif, centres_finaux)
%
% Auteur: Mounir Bendali-Braham
% Date de création : Mai 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialisation
nbCentres = size(centres_finaux, 3);
nbObjets = size(W, 3);
nbObjets_parClasse = ones(1, nbCentres);
inertiesRV_intra = ones(1, nbCentres);
S = calculMatriceS(W, D);
normesW = diag(S);
S_centres = calculMatriceS(centres_finaux, D);
normes_centres = diag(S_centres);

% S_prime = ones(nbObjets, nbCentres);
for t = 1:nbObjets
    for k = 1:nbCentres
        S_prime(t,k) = ...
            prod_hs(W(:,:,t), ...
            centres_finaux(:,:,k), D);
    end;
end;
distRV = ones(nbObjets, nbCentres);
for t = 1:nbObjets
    for k = 1:nbCentres
        distRV(t,k) = S_prime(t,k)/(sqrt(normes_centres(k))*sqrt(normesW(t)));
    end;
end;

%% Programme
for c = 1:nbCentres
    objets_centre = find(index_classif==c);
    nbObjets_parClasse(c) = size(objets_centre, 2);
    disp(sum(distRV(objets_centre, c)));
    inertiesRV_intra(c) = sum(distRV(objets_centre, c))/nbObjets;
end;

inertieRVglobale = sum(inertiesRV_intra);
end
