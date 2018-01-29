function [Ac] = normaliser(A)
    mean_A = mean(A);
    std_A = std(A);
    %--------------------------------
    % Centrage des donn?es
    %--------------------------------
    UN = ones(size(A));
    Me = UN * diag(mean_A);
    Ecart_type = UN * diag(std_A);
    ec  = 1./Ecart_type;
    Ac  = (A - Me).*ec;
end

