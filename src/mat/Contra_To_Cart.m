function [Ucat_x, Ucat_y] = ...
    Contra_To_Cart(M, N, Ucont_x, Ucont_y)

    % This function convert contravariant to Catersian 
    % components of flux Ucont.
    % In this case it is just the average

    Ucat_x = zeros(M, N);
    Ucat_y = zeros(M, N);

    for ii = 1:M
        for jj = 1:N
            Ucat_x(ii, jj) = (Ucont_x(ii, jj) + Ucont_x(ii+1, jj)) / 2 ;
            Ucat_y(ii, jj) = (Ucont_y(ii, jj) + Ucont_y(ii, jj+1)) / 2 ;
        end
    end

end