function [Ucat_x, Ucat_y] = ...
    Contra_To_Cart(Ucont_x, Ucont_y, M2, N2)
    % This function convert contravariant to Catersian 
    % components of flux Ucont.
    % In this case it is just the average

    Ucat_x = zeros(N2, M2);
    Ucat_y = zeros(N2, M2);

    for ii = 1:M2
        for jj = 1:N2
            Ucat_x(jj, ii) = (Ucont_x(jj, ii) + Ucont_x(jj, ii+1)) / 2 ;
            Ucat_y(jj, ii) = (Ucont_y(jj, ii) + Ucont_y(jj+1, ii)) / 2 ;
        end
    end
%{
    %-------------------- x- direction -----------------------------
    % Just average the mid-point to get the Cartesian components
    for i = 2:M-1  
        for j = 2:N-1
        % It is identical with volume flux of Cartesian node
        Ucat_x(i,j) = (Ucont_x(i,j) + Ucont_x(i-1,j)) / 2 ;
        end
    end


    %-------------------- y- direction -----------------------------
    % Just average the mid-point to get the Cartesian components
    for i = 2:M-1 
        for j = 2:N-1
        % It is identical with volume flux of Cartesian node
        Ucat_y(i,j) = (Ucont_y(i,j) + Ucont_y(i,j-1)) / 2 ;
        end
    end

    % For the boundary - set to zero
    % For x - contravariant
    for i = 1:M
        Ucat_x(i,N) = 0;
        Ucat_x(i,1) = 0;
    end

    for j = 1:N
        Ucat_x(1,j) = 0;
        Ucat_x(M,j) = 0;
    end
    % For y- contravariant
    for i = 1:M
        Ucat_y(i,N) = 0;
        Ucat_y(i,1) = 0;
    end

    for j = 1:N
        Ucat_y(1,j) = 0;
        Ucat_y(M,j) = 0;
    end
%}
end