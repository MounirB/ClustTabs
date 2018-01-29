function [ Wcomp, alpha_t ] = compromis(W, S, Delta, VaP, VeP, norm)
%% Fonction de calcul du compromis pour la methode STATIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables d'entrée :
% W = Matrice avec les objets des t etudes
% S = Matrice des produits scalaire entre les Objets (W) affectes par le poids
%     Delta
% Delta = Matrice diagonal avec les poids (pi_i)
% norm = Recours ou non à la normalisation des objets
% VaP = Valeurs propres du matrice S
% VeP = Vecteurs propres du matrice S
%
% Variables de sortie :
% Wcomp = Matrice avec le Compromis entre les objets W
%
% Usage :
% [ Wcomp, alpha_t ] = compromis(W, S, Delta, VaP, VeP, norm)
%
% Auteurs : Larbi Mouchou, Rodrigo Andres Rivera Martinez, Mounir Bendali-Braham, Nafise Gouard
% Date de création : Mars 2017
% Edit (Avril 2017): Mounir Bendali-Braham
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Expression definitive du compromis
    [L,C,T] = size(W);
    pi_t=diag(Delta);
    SS = diag(S);
    gamma = (VeP(:,1));

    Wcomp = 0;
    alpha_t = zeros(1,T);
    if ~norm
        alpha_c = (pi_t'*sqrt(SS));
        for i = 1:T
            alpha_t(i) = ((1/sqrt(VaP(1)))*alpha_c*sqrt(pi_t(i))*gamma(i));
%             alpha_t(i) = 1/size(W,3);
            Wcomp = Wcomp + (alpha_t(i)*W(:,:,i));
        end
    else
        for i = 1:T
            alpha_t(i) = ((1/sqrt(VaP(1)))*(sqrt(pi_t(i)))*gamma(i));
%             alpha_t(i) = 1/size(W,3);
            Wcomp = Wcomp + (alpha_t(i)*W(:,:,i));
        end
    end
end


