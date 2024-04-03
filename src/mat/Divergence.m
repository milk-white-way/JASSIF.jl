function [Div] = Divergence(M, N, Ucont_x, Ucont_y, dx ,dy)
    % This function does calculate the divergence at every grid node in the
    % Computational domain

    Div = zeros(M, N);

    for ii = 1:M
        for jj = 1:N
            Div(ii, jj) = (Ucont_x(ii+1, jj) - Ucont_x(ii, jj)) / dx ...
                        + (Ucont_y(ii, jj+1) - Ucont_y(ii, jj)) /dy;        
        end
    end

end