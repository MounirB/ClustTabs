function [ index_classif, centres_finaux ] = classif_blocsVars( X, nbCenters, maxIter, Delta, M, D, norm )
%% Fonction de classification des donn�es multidimensionnelles �volutives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables d'entr�e:
% X = Conteneur des T �tudes : donn�es multidimensionnelles �volutives
% k = Nombre de classes k
% maxIter = Nombre maximal d'it�rations
% Delta = Matrice diagonale des poids des �tudes
% D = M�trique des poids des individus
% M = M�trique des poids des variables
% norm = Bool�en soulignant le recours ou non � la norme 
%
% Variables de sortie:
% index_classif = index indiquant la classe de chaque �tude
%
% Usage:
% [ index_classif ] = classif_etudes( X, k, maxIter, noms_etudes, Delta, M, D, norm )
%
% Auteur: Mounir Bendali-Braham
% Date de cr�ation: Avril 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialisation
    [~, ~, nbEtudes] = size(X);
    [W] = calculObjets(X, M, norm);
    index_centres_init = select_init_centers( nbEtudes, nbCenters);
    centres_Wk = W(:, :, index_centres_init);
    S = calculMatriceS(W, D);
    s = diag(S);
    % Indicateurs de convergeance
    changement = true;
    iter = 0;
    index_classif = ones(1, nbEtudes);
    
%% Boucle principale
    while(iter < maxIter && changement) 
        %% Affectation
        S_centres = calculMatriceS(centres_Wk, D);
        s_centres = diag(S_centres);
        index_classif = affecter(W, centres_Wk, D, s, s_centres);
        %% Minimisation
        centres_prec = centres_Wk;
        centres_Wk = minimiser(W, S, D, Delta, nbCenters, index_classif, norm);
        %% Stabilisation
        iter = iter + 1;
        if(isequal(centres_prec, centres_Wk))
             changement = false;
        end;
    end;
    centres_finaux = centres_Wk;
    
%% Affichage
%     varetude = {'REYKJAVIK','NEYWORK', 'DAKAR'};
%     [Co,~,~,~,~,~,~,~] = statis_inter (X,M,Delta,norm,D,varetude);
%     colors = jet;
%     colorsDelim = round(length(colors)/nbCenters);
%     figure;
%     hold on;
%     for(c=1:nbCenters)
%         graph = plot(Co(find(index_classif==c), 1), Co(find(index_classif==c), 2), '.g', 'MarkerSize', 10);
%         set(graph, 'Color', colors(colorsDelim*c, :));
%     end;
end

