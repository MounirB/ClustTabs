%% Programme principal de classification Kmeans des variables mfeatures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ce script fait appel à l'algorithme ClustVar 'Som_cov_2' (Clustering de variables
% tel que décrit par Vigneau et El Qannari, mais sans qu'il soit initialisé par une 
% procédure de classification hiérarchique)
%
% Auteur: Mounir Bendali-Braham
% Date de création : Septembre 2017
% Date de dernière modification : 26 Novembre 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% Importation des données mfeatures
fac = importdata('./mfeat/mfeat-fac');
fou = importdata('./mfeat/mfeat-fou');
zer = importdata('./mfeat/mfeat-zer');
mor = importdata('./mfeat/mfeat-mor');

%% Individus

Starts = 1;
Ends = 2000;
nbIndividus = 2000;

%% Chargement des données

etude(1) = {fac(Starts:Ends, 1:216)};
etude(2) = {fou(Starts:Ends, 1:76)};
etude(3) = {zer(Starts:Ends, 1:47)};
etude(4) = {mor(Starts:Ends, 1:6)};

%% Initialisation

nbEtudes = size(etude, 2);

%% Normaliser les données

for t= 1:nbEtudes
    etudeN(t) = {normaliser(cell2mat(etude(t)))};
end;

%% Création de la matrice des données dataMatrix

dataMatrix = [];

for t= 1:nbEtudes
    dataMatrix = [dataMatrix, cell2mat(etudeN(t))];
end;

tDataMatrix = dataMatrix'; % Données transposées pour classer les variables

%% Matrice de corrélations et visualisation

corr = corrcoef(dataMatrix);
imagesc(corr);
caxis([0, 1]);
colorbar;
title('Matrice des corrélations des variables');

%% Préparation de la classification

maxIter = 1000;
nbCentres = 4;
varNbReferents = [3:8]; % On fait varier le nombre de référents
nbVariations_nbReferents = length(varNbReferents); % Le nombre de variations du nombre de référents
nbExecutions = 50;

%% Préparation du calcul des performances

nbClasse1 = 216;
nbClasse2 = 76;
nbClasse3 = 47;
nbClasse4 = 6;
indice = 1;

for elt = 1:nbClasse1
    Classe_1(elt, :) = tDataMatrix(indice, :);
    indice = indice + 1;
end;

for elt = 1:nbClasse2
    Classe_2(elt, :) = tDataMatrix(indice, :);
    indice = indice + 1;
end;

for elt = 1:nbClasse3
    Classe_3(elt, :) = tDataMatrix(indice, :);
    indice = indice + 1;
end;

for elt = 1:nbClasse4
    Classe_4(elt, :) = tDataMatrix(indice, :);
    indice = indice + 1;
end;

% Labellisation des classes

startsEnds(1, :) = [1, 216];
startsEnds(2, :) = [217, 292];
startsEnds(3, :) = [293, 339];
startsEnds(4, :) = [340, 345];

%% Classification et calcul des performances

for nbRef = 1:nbVariations_nbReferents
    for init = 1:nbExecutions
        [bmus,sMap,iterations1_tauy,iterations2_tauy ,affectation ] = Som_cov_2(tDataMatrix, varNbReferents(nbRef), 1, 30, 2, 0, 0);
        
        index_classif_parInit(init, :) = bmus;
    end;
    
    indexes_parNbRef(nbRef) = {index_classif_parInit};
  
    %% Matrices de contingence avec la partition cible (MatCont)
    % Calcul des performances en se basant sur les MatCont
    
    partition_cible = [ones(1,216), 2*ones(1,76), 3*ones(1,47), 4*ones(1,6)];
    
    for numPartition = 1: nbExecutions
        for numClasse = 1: varNbReferents(nbRef)
            classe = find(index_classif_parInit(numPartition, :)==numClasse);
            for numClasse_cible = 1: nbCentres
                classe_cible = find(partition_cible == numClasse_cible);
                nbCommuns_pourCible = length(intersect(classe_cible, classe));
                nbCommuns_pourClasse = length(intersect(classe, classe_cible));
                matrice_contingence(numClasse, numClasse_cible) = min(nbCommuns_pourCible, nbCommuns_pourClasse);
            end;
        end;
        matCont_parInit(numPartition) = {matrice_contingence};
        clearvars matrice_contingence;
    end;
    matCont_parInit_parNbRef(nbRef) = {matCont_parInit};
    
    clearvars index_classif_parInit;
    
    % Calcul des quantités N11, N01, N10, N00, et ensuite des performances
    
    n11_parInit = 0;
    n10_parInit = 0;
    n01_parInit = 0;
    n00_parInit = 0;
    
    for numPartition = 1: nbExecutions
        %Calcul de n11
        sum_nObservs_commun = 0; %Nombre d'effectifs d'observations en commun
        
        matrice_contingence = cell2mat(matCont_parInit(numPartition));
        for numClasse = 1: varNbReferents(nbRef)
            for numClasse_cible = 1: nbCentres
                n_observs = matrice_contingence(numClasse, numClasse_cible);
                sum_nObservs_commun = sum_nObservs_commun + n_observs*(n_observs-1);
            end;
        end;
        n11_parInit = sum_nObservs_commun * 0.5;

        % Calcul de n01
        sum_nObservs = 0;
        marginal_freq01 = 0;
        sum_nObservs_carre = 0;
        for numClasse = 1: varNbReferents(nbRef)
            for numClasse_cible = 1: nbCentres
                n_observs = matrice_contingence(numClasse, numClasse_cible);
                sum_nObservs = sum_nObservs + n_observs;
                sum_nObservs_carre = sum_nObservs_carre + n_observs*n_observs;
            end;
            marginal_freq01 = marginal_freq01 + sum_nObservs*sum_nObservs;
            sum_nObservs = 0;
        end;
        n01_parInit = (marginal_freq01 - sum_nObservs_carre)* 0.5;
        
        % Calcul n10
        marginal_freq10 = 0;
        sum_nObservs_carre = 0;
        for numClasse_cible = 1: nbCentres
            for numClasse = 1: varNbReferents(nbRef)
                n_observs = matrice_contingence(numClasse, numClasse_cible);
                sum_nObservs = sum_nObservs + n_observs;
                sum_nObservs_carre = sum_nObservs_carre + n_observs*n_observs;
            end;
            marginal_freq10 = marginal_freq10 + sum_nObservs*sum_nObservs;
            sum_nObservs = 0;
        end;
        n10_parInit = (marginal_freq10 - sum_nObservs_carre)* 0.5;
        
        % Calcul n00
        sum_nObservs = 0;
        for numClasse_cible = 1: nbCentres
            for numClasse = 1: varNbReferents(nbRef)
                n_observs = matrice_contingence(numClasse, numClasse_cible);
                sum_nObservs = sum_nObservs + n_observs;
            end;
        end;

        n00_parInit = 0.5*(sum_nObservs.^2 + sum_nObservs_carre - (marginal_freq01 + marginal_freq10));
        
        rappel_parInit_parRapportPartitionCible(numPartition) = n11_parInit/(n11_parInit + n01_parInit);
        precision_parInit_parRapportPartitionCible(numPartition) = n11_parInit/(n11_parInit + n10_parInit);
        f_mesure_parInit_parRapportPartitionCible(numPartition) = (2 * rappel_parInit_parRapportPartitionCible(numPartition) * precision_parInit_parRapportPartitionCible(numPartition))/(rappel_parInit_parRapportPartitionCible(numPartition) + precision_parInit_parRapportPartitionCible(numPartition));
    end;
    
    rappel_parInitRapPartCible_parNbRef(nbRef) = {rappel_parInit_parRapportPartitionCible};
    precision_parInitRapPartCible_parNbRef(nbRef) = {precision_parInit_parRapportPartitionCible};
    f_mesure_parInitRapPartCible_parNbRef(nbRef) = {f_mesure_parInit_parRapportPartitionCible};
    
end;

% Affichage des performances
n_nbRef = 6;

for nbRef = 1:n_nbRef
    rappel = cell2mat(rappel_parInitRapPartCible_parNbRef(nbRef));
    precision = cell2mat(precision_parInitRapPartCible_parNbRef(nbRef));
    f_mesure = cell2mat(f_mesure_parInitRapPartCible_parNbRef(nbRef));
    
    figure;
    subplot(1, 3, 1);
    boxplot(rappel);
    xlabel('Rappel');
    ylim([0 1.1]);
    
    subplot(1, 3, 2);
    boxplot(precision);
    xlabel('Precision');
    ylim([0 1.1]);
    
    subplot(1, 3, 3);
    boxplot(f_mesure);
    xlabel('F-mesure');
    ylim([0 1.1]);
    
    moy(1) = mean(rappel);
    moy(2) = mean(precision);
    moy(3) = mean(f_mesure);
    moys = [moy(1), moy(2), moy(3)];
    
    moyennes(nbRef) = {moys};
end;
