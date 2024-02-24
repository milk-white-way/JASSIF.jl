function [A] = Poisson_LHS(M,N,dx,dy)
% This function does assemble Poisson LHS 
% This is Laplacian in 4 points
% Memory allocation for A
A = spalloc(M*N,M*N, 5);
sp_ii = [];
sp_jj = [];
sp_vals = [];

for counter = 1 :M*N
    
            [i j] = lidx(counter,M,N);    
    
            indexN = glidx(i  ,j-1,M,N);                        
            indexE = glidx(i-1,j  ,M,N);                        
            indexP = glidx(i  ,j  ,M,N);                        
            indexW = glidx(i+1,j  ,M,N);                        
            indexS = glidx(i,j+1  ,M,N);                        
            
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
            
        else
            ii = zeros(1,1);
            jj = zeros(1,1);
            vals = zeros(1,1);
    
            ii(1) = indexP;
            jj(1) = indexP;
            vals(1) = 1;
        end
        
        sp_ii = [sp_ii ii];
        sp_jj = [sp_jj jj];
        sp_vals = [sp_vals vals];
    
end

A = sparse(sp_ii,sp_jj,sp_vals);
