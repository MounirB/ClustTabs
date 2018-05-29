%% Programme principal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script faisant appel à l'algorithme de clustering ClustTabs sur les données 
% multiple features
%
% Auteur: Mounir Bendali-Braham
% Date de création : Septembre 2017
% Date de dernière modification : 30 Mai 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% Importation des données mfeatures
fac = importdata('./mfeat/mfeat-fac');
fou = importdata('./mfeat/mfeat-fou');
kar = importdata('./mfeat/mfeat-kar');
mor = importdata('./mfeat/mfeat-mor');
pix = importdata('./mfeat/mfeat-pix');
zer = importdata('./mfeat/mfeat-zer');

%% Individus

Starts = 1;
Ends = 2000;
nbIndividus = 2000;

%% Découpage en études

etude(1) = {fac(Starts:Ends, 1:40)};
etude(2) = {fac(Starts:Ends, 41:80)};
etude(3) = {fac(Starts:Ends, 81:120)};
etude(4) = {fac(Starts:Ends, 121:160)};
etude(5) = {fac(Starts:Ends, 161:216)};
etude(6) = {fou(Starts:Ends, 1:30)};
etude(7) = {fou(Starts:Ends, 31:76)};
etude(8) = {zer(Starts:Ends, 1:20)};
etude(9) = {zer(Starts:Ends, 21:47)};
etude(10) = {mor(Starts:Ends, 1:6)};

%% Initialisation

nbEtudes = size(etude, 2);
D =(1/nbIndividus)*eye(nbIndividus);

%% Normaliser les données

for t= 1:nbEtudes
    etudeN(t) = {normaliser(cell2mat(etude(t)))};
end;

%% Création des objets représentatifs des études W

W = calculObjets_bisFromCells(etudeN);
RV = calculMatriceRV(W, D);

%% Visualisation de la matrice RV

imagesc(RV);
caxis([0, 1]);
colorbar;

%% Préparation de la classification

maxIter = 1000;
Delta = (1/size(W,3))*eye(size(W,3));
norm = 0;
nbCentres = 4;
varNbReferents = [3:8]; % On fait varier le nombre de référents
nbVariations_nbReferents = length(varNbReferents); % Le nombre de variations du nombre de référents
nbExecutions = 50;

%% Préparation du calcul des performances

nbClasse1 = 5;
nbClasse2 = 2;
nbClasse3 = 2;
nbClasse4 = 1;
indice = 1;

for elt = 1:nbClasse1
    Classe_1(:, :, elt) = W(:, :, indice);
    indice = indice + 1;
end;

for elt = 1:nbClasse2
    Classe_2(:, :, elt) = W(:, :, indice);
    indice = indice + 1;
end;

for elt = 1:nbClasse3
    Classe_3(:, :, elt) = W(:, :, indice);
    indice = indice + 1;
end;

for elt = 1:nbClasse4
    Classe_4(:, :, elt) = W(:, :, indice);
    indice = indice + 1;
end;

% Labellisation des classes

startsEnds(1, :) = [1, 5];
startsEnds(2, :) = [6, 7];
startsEnds(3, :) = [8, 9];
startsEnds(4, :) = [10, 10];

%% Classification et calcul des performances

for nbRef = 1:nbVariations_nbReferents
    for init = 1:nbExecutions
        [ index_classif, centres_finaux ] = classif_blocsVars_bis( W, varNbReferents(nbRef), maxIter, Delta, D, norm );
        
        [inertieRVglobale, inertiesRV_intra, nbObjets_parClasse] = inertieIntraClasse( W, D, index_classif, centres_finaux );
        
        globalRV_init(init) = inertieRVglobale;
        nbObjetsParClasse(init, :) = nbObjets_parClasse;
        inertiesRV_intraParInit(init, :) = inertiesRV_intra;
        index_classif_parInit(init, :) = index_classif;
    end;
    
    globalRV_parNbRef(nbRef) = {globalRV_init};
    nbObjets_parNbRef(nbRef) = {nbObjetsParClasse};
    inertiesRV_parNbRef(nbRef) = {inertiesRV_intraParInit};
    indexes_parNbRef(nbRef) = {index_classif_parInit};
    
    clearvars globalRV_init nbObjetsParClasse inertiesRV_intraParInit index_classif_parInit;
    
    %% Calcul des performances
    index_classif_parInit = cell2mat(indexes_parNbRef(nbRef));

    for numPartition = 1: nbExecutions
        for numClasse = 1: nbCentres
            classe_label(numClasse) = mode(index_classif_parInit(numPartition, startsEnds(numClasse, 1): startsEnds(numClasse, 2)));
        end;
        classeLabellisation_parInit(numPartition, :) = classe_label;
    end;
    
    classeLabellisation_parNbRef(nbRef) = {classeLabellisation_parInit};
    
    % Rappel et Precision

    for numPartition = 1: nbExecutions
        partition = index_classif_parInit(numPartition, :);
        classe_label = classeLabellisation_parInit(numPartition, :);
        rappel_parInit(numPartition) = 0;
        precision_parInit(numPartition) = 0;
        for numClasse = 1: nbCentres
            classe = index_classif_parInit(numPartition, startsEnds(numClasse, 1): startsEnds(numClasse, 2));
            rappel_parClasseParInit(numPartition, numClasse) = length(find(classe==classe_label(numClasse)))/length(classe);
            precision_parClasseParInit(numPartition, numClasse) = length(find(classe==classe_label(numClasse)))/ length(find(partition==classe_label(numClasse)));
            rappel_parInit(numPartition) = rappel_parInit(numPartition) + rappel_parClasseParInit(numPartition, numClasse);
            precision_parInit(numPartition) = precision_parInit(numPartition) + precision_parClasseParInit(numPartition, numClasse);
        end;
        rappel_parInit(numPartition) = rappel_parInit(numPartition)/nbCentres;
        precision_parInit(numPartition) = precision_parInit(numPartition)/nbCentres;
    end;
    
    rappel_parClasseParInitParNbRef(nbRef) = {rappel_parClasseParInit};
    precision_parClasseParInitParNbRef(nbRef) = {precision_parClasseParInit};
    rappel_parInitParNbRef(nbRef) = {rappel_parInit};
    precision_parInitParNbRef(nbRef) = {precision_parInit};
    
    rappel_moyennes_parNbRef(nbRef) = {mean(rappel_parClasseParInit)};
    precision_moyennes_parNbRef(nbRef) = {mean(precision_parClasseParInit)};

    % F-mesure par classe par Exécution

    for numPartition = 1: nbExecutions
        f_mesure_parInit(numPartition) = 0;
        for numClasse = 1: nbCentres
            f_mesure_parClasseParInit(numPartition, numClasse) = (2*precision_parClasseParInit(numPartition, numClasse)*rappel_parClasseParInit(numPartition, numClasse))/(precision_parClasseParInit(numPartition, numClasse)+rappel_parClasseParInit(numPartition, numClasse));
            f_mesure_parInit(numPartition) = f_mesure_parInit(numPartition) + f_mesure_parClasseParInit(numPartition, numClasse);
        end;
        f_mesure_parInit(numPartition) = f_mesure_parInit(numPartition)/nbCentres;
    end;
    
    f_mesure_parClasseParInitParNbRef(nbRef) = {f_mesure_parClasseParInit};
    
    f_mesure_moyennes_parNbRef(nbRef) = {mean(f_mesure_parClasseParInit)};
    
    %% Calcul du nombre de perfect matchs PMs pour les nbExecutions par indice de performance
    
    for numPartition = 1:nbExecutions
        valeurRappelMaximale = max(rappel_parInit);
        nbPM_rappel = length(find(rappel_parInit == valeurRappelMaximale));
    end;
    
    nbPM_rappel_parNbRef(nbRef) = {nbPM_rappel};
    
    for numPartition = 1:nbExecutions
        valeurPrecisionMaximale = max(precision_parInit);
        nbPM_precision = length(find(precision_parInit == valeurPrecisionMaximale));
    end;
    
    nbPM_precision_parNbRef(nbRef) = {nbPM_precision};
    
    for numPartition = 1:nbExecutions
        valeurFMesureMaximale = max(f_mesure_parInit);
        nbPM_f_mesure = length(find(f_mesure_parInit == valeurFMesureMaximale));
    end;
    
    nbPM_f_mesure_parNbRef(nbRef) = {nbPM_f_mesure};
    
    %% Matrices de contingence avec la partition cible (MatCont)
    % Calcul des performances en se basant sur les MatCont
    
    partition_cible = [1, 1, 1, 1, 1, 2, 2, 3, 3, 4];
    
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
    
    %% Affichages des performances par exécution
    nbMesures = 3; % Rappel, Précision, F-mesure
    
    mesures_parClasseParInit = {...,
        rappel_parClasseParInit,
        precision_parClasseParInit,
        f_mesure_parClasseParInit};
    
    mesures_text = {'Rappel', 'Precision', 'F-mesure'}; 
    
    figure;
    for numMesure = 1: nbMesures
        subplot(1, nbMesures, numMesure);
        mesure_parClasseParInit = cell2mat(mesures_parClasseParInit(numMesure));
        hold on;
        boxplot(mesure_parClasseParInit);
        xlabel('Classes');
        ylabel(mesures_text(numMesure));
    end;
    hold off;
end;

% After
n_nbRef = 6;

% for nbRef = 1:n_nbRef
%     f_mesure = cell2mat(f_mesure_parInitRapPartCible_parNbRef(nbRef));
%     f_mesure = f_mesure * 4;
%     f_mesure_parInitRapPartCible_parNbRef(nbRef) = {f_mesure};
% end;

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
