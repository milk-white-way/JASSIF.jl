function [P_RHS] = Poisson_RHS_Neumann(Ucont_x,Ucont_y,dx,dy,dt)

M = length(Ucont_x(:,1));
N = length(Ucont_x(1,:));

P_Div = Divergence(Ucont_x, Ucont_y,dx ,dy);

for i = 1:M
    for j = 1:N
 
      
        index = glidx(i,j,M,N);
        
        if (i >=2 && i<=M-1 && j>=2 && j<=N-1)
            
            P_RHS(index) = P_Div(i,j);            
                        
        else
            % Apply Neumann boundary condition for Poisson equation
            P_RHS(index) = 0;                           
        end       
    end    
end

% Check if solvability condition is met !
Summation = 0;
for i = 1:M
    for j = 1:N    
    index = glidx(i,j,M,N);
    
        if (i >=2 && i<=M-1 && j>=2 && j<=N-1)
            
            % Not the fix pressure location
               Summation = Summation + P_RHS(index);
               
        end
    end
end
Summation = Summation * 1.5 / dt;

if (Summation > 1e-8) 
    fprintf(' Summation of RHS = %8.5f -> Solvability of Poisson equation is invalidated \n', Summation);
end

% Scale back to time discretization scheme
% Adam-Bashforth scheme
P_RHS = P_RHS * 1.5 / dt;


