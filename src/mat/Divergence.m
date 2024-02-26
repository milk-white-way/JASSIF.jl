function [Div] = Divergence(Ucont_x, Ucont_y, M, N, dx ,dy)
    % This function does calculate the divergence at every grid node in the
    % Computational domain

    Div = zeros(N, M);

    for ii = 1:M
        for jj = 1:N
            Div(jj, ii) = (Ucont_x(jj, ii+1) - Ucont_x(jj, ii)) / dx ...
                        + (Ucont_y(jj+1, ii) - Ucont_y(jj, ii)) /dy;        
        end
    end

end