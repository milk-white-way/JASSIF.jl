function [A] = Poisson_LHS_Neumann(M,N,dx,dy)
% This function does assemble Poisson LHS 
% This is Laplacian in 4 points
% Memory allocation for A
A = spalloc(M*N,M*N, 5);
sp_ii = [];
sp_jj = [];
sp_vals = [];

for counter = 1 :M*N
    
    % Get local indices
            [i j] = lidx(counter,M,N);             
            

            % Find the neighbours            
                indexN = glidx(i  ,j-1,M,N);                                    
                indexE = glidx(i-1,j  ,M,N);                                                
                indexP = glidx(i  ,j  ,M,N);                                    
                indexW = glidx(i+1,j  ,M,N);                                                
                indexS = glidx(i,j+1  ,M,N);                           
        
         % Inside the real domain
        if (i>=2 && i <=M-1 && j >= 2 && j<=N-1)            
            
            ii = zeros(1,5);
            jj = zeros(1,5);
            vals = zeros(1,5);           
            
            % point N
            ii(1) = indexP;
            jj(1) = indexN;
            vals(1) = 1 /dy^2;
    
            % Point E
            ii(2) = indexP;
            jj(2) = indexE;
            vals(2) = 1 / dx^2;
        
            % Point P            
            ii(3) = indexP;
            jj(3) = indexP;
            vals(3) = - 2 / dx^2 - 2 /dy^2;            
            
            % Point W
            ii(4) = indexP;
            jj(4) = indexW;
            vals(4) = 1 / dx^2;
            
            % Point S
            ii(5) = indexP;
            jj(5) = indexS;
            vals(5) = 1 / dy^2;            
            
             %%-------- Assign phi at p(2,2)  = 0 for matrix not to be
             %%singular
            
            if (i== 2 && j==2)
               vals(3) = 0;                           
            end
            
            
        else
                           
                   
            %%---------------------- Boundaries---------------------
            % i == 1
            if (i == 1 && j ~=1 && j~=N)               
                  
                    ii = zeros(1,2);
                    jj = zeros(1,2);
                    vals = zeros(1,2);
      
                    % Point P            
                    ii(1) = indexP;
                    jj(1) = indexP;
                    vals(1) = - 1;            
            
                    % Point W
                    ii(2) = indexP;
                    jj(2) = indexW;
                    vals(2) = 1;
            
            end 
            
            % i == M
            if (i == M && j~=1 && j~=N)               
                  
                    ii = zeros(1,2);
                    jj = zeros(1,2);
                    vals = zeros(1,2);
      
                    % Point P            
                    ii(1) = indexP;
                    jj(1) = indexP;
                    vals(1) = - 1;            
            
                   % Point E
                    ii(2) = indexP;
                    jj(2) = indexE;
                    vals(2) = 1;
            
            end % 
            
             %j == 1
            if (j == 1 && i ~=1 && i~=M )               
                  
                    ii = zeros(1,2);
                    jj = zeros(1,2);
                    vals = zeros(1,2);
      
                    % Point P            
                    ii(1) = indexP;
                    jj(1) = indexP;
                    vals(1) = - 1;            
            
                    % Point S
                    ii(2) = indexP;
                    jj(2) = indexS;
                    vals(2) = 1;
            
            end % 
            
            % j == N
            if (j == N && i ~=1 && i~=M )               
                  
                    ii = zeros(1,2);
                    jj = zeros(1,2);
                    vals = zeros(1,2);
      
                    % Point P            
                    ii(1) = indexP;
                    jj(1) = indexP;
                    vals(1) = - 1;            
            
                    % Point N
                    ii(2) = indexP;
                    jj(2) = indexN;
                    vals(2) = 1;
            
            end %        
            
            
             %% Default
             
             if (i==1 &&j==1 || i==1 &&j==N || i==M &&j==1||i==M &&j==N)
                 
                    ii = zeros(1,1);
                    jj = zeros(1,1);
                    vals = zeros(1,1);
      
                    % Point P            
                    ii(1) = indexP;
                    jj(1) = indexP;
                    vals(1) = 1;        
             end
        end
        
        % Assemble the matrix 
        sp_ii = [sp_ii ii];
        sp_jj = [sp_jj jj];
        sp_vals = [sp_vals vals];
    
end

A = sparse(sp_ii,sp_jj,sp_vals);
