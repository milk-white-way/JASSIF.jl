function [Ucat_x, Ucat_y] = ...
    Contra_To_Cart(Ucont_x, Ucont_y, M, N)

    % This function convert contravariant to Catersian 
    % components of flux Ucont.
    % In this case it is just the average

    Ucat_x = zeros(N, M);
    Ucat_y = zeros(N, M);

    for ii = 1:M
        for jj = 1:N
            Ucat_x(jj, ii) = (Ucont_x(jj, ii) + Ucont_x(jj, ii+1)) / 2 ;
            Ucat_y(jj, ii) = (Ucont_y(jj, ii) + Ucont_y(jj+1, ii)) / 2 ;
        end
    end

end