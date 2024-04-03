function [A] = Poisson_LHS_Neumann(M, N, dx, dy)
% This function does assemble Poisson LHS 
% This is Laplacian in 4 points
    % Memory allocation for A
    A = spalloc(M*N, M*N, 5);
    sp_ii = [];
    sp_jj = [];
    sp_vals = [];

    for glb_idx = 1:M*N
        % Get local indices
        [i, j] = lidx(M, glb_idx);             
        loc_idx_P = glb_idx;

        ii = zeros(1,5);
        jj = zeros(1,5);
        vals = zeros(1,5);           

        %% ---------------------- Boundaries --------------------------
        % ---------------------- Corners --------------------------
        if ( i==1 && j==1 ) 
        % South-West corner
            loc_idx_N = glidx(M, N, i, j+1);    % normal
            loc_idx_S = glidx(M, N, i, N);      % mirrored
            loc_idx_E = glidx(M, N, i+1, j);    % normal
            loc_idx_W = glidx(M, N, M, j);      % mirrored

        elseif ( i==1 && j==N )
        % North-West corner
            loc_idx_N = glidx(M, N, i, 1);      % mirrored
            loc_idx_S = glidx(M, N, i, j-1);    % normal
            loc_idx_E = glidx(M, N, i+1, j);    % normal
            loc_idx_W = glidx(M, N, M, j);      % mirrored
                
        elseif ( i==M && j==1 )
        % South-East corner
            loc_idx_N = glidx(M, N, i, j+1);    % normal
            loc_idx_S = glidx(M, N, i, N);      % mirrored
            loc_idx_E = glidx(M, N, 1, j);      % mirrored
            loc_idx_W = glidx(M, N, i-1, j);    % normal

        elseif ( i==M && j==N )
        % North-East corner                    
            loc_idx_N = glidx(M, N, i, 1);      % mirrored
            loc_idx_S = glidx(M, N, i, j-1);    % normal
            loc_idx_E = glidx(M, N, 1, j);      % mirrored
            loc_idx_W = glidx(M, N, i-1, j);    % normal

        % ---------------------- Edges --------------------------
        elseif ( i==1 && j~=1 && j~=N )               
        % West minus corners
            loc_idx_N = glidx(M, N, i, j+1);    % normal
            loc_idx_S = glidx(M, N, i, j-1);    % normal 
            loc_idx_E = glidx(M, N, i+1, j);    % normal
            loc_idx_W = glidx(M, N, M, j);      % mirrored

        elseif ( i==M && j~=1 && j~=N )               
        % East minus corners
            loc_idx_N = glidx(M, N, i, j+1);    % normal
            loc_idx_S = glidx(M, N, i, j-1);    % normal
            loc_idx_E = glidx(M, N, 1, j);      % mirrored
            loc_idx_W = glidx(M, N, i-1, j);    % normal

        elseif ( j==1 && i~=1 && i~=M )               
        % South minus corners
            loc_idx_N = glidx(M, N, i, j+1);    % normal
            loc_idx_S = glidx(M, N, i, N);      % mirrored
            loc_idx_E = glidx(M, N, i+1, j);    % normal
            loc_idx_W = glidx(M, N, i-1, j);    % normal

        elseif ( j==N && i~=1 && i~=M )               
        % North minus corners
            loc_idx_N = glidx(M, N, i, 1);      % mirrored
            loc_idx_S = glidx(M, N, i, j-1);    % normal
            loc_idx_E = glidx(M, N, i+1, j);    % normal
            loc_idx_W = glidx(M, N, i-1, j);    % normal
                    
        %% ---------------------- Internal points ---------------------
        else            
            % Find the neighbours
            loc_idx_N = glidx(M, N, i  , j+1);  % normal
            loc_idx_S = glidx(M, N, i  , j-1);  % normal
            loc_idx_E = glidx(M, N, i+1, j  );  % normal
            loc_idx_W = glidx(M, N, i-1, j  );  % normal
        end

        % point N
        ii(1) = loc_idx_P;
        jj(1) = loc_idx_N;
        vals(1) = 1 / dy^2;
                
        % Point E
        ii(2) = loc_idx_P;
        jj(2) = loc_idx_E;
        vals(2) = 1 / dx^2;
                
        % Point P            
        ii(3) = loc_idx_P;
        jj(3) = loc_idx_P;
        vals(3) = - 2 / dx^2 - 2 /dy^2;            
                
        % Point W
        ii(4) = loc_idx_P;
        jj(4) = loc_idx_W;
        vals(4) = 1 / dx^2;
                
        % Point S
        ii(5) = loc_idx_P;
        jj(5) = loc_idx_S;
        vals(5) = 1 / dy^2;            
            
        % Assemble the matrix 
        sp_ii = [sp_ii ii];
        sp_jj = [sp_jj jj];
        sp_vals = [sp_vals vals];
            
    end
    A = sparse(sp_ii, sp_jj, sp_vals);
end