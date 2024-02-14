function [P_RHS] = Poisson_RHS(Ucont_x,Ucont_y,dx,dy,dt)

M = length(Ucont_x(:,1));
N = length(Ucont_x(1,:));

P_Div = Divergence(Ucont_x, Ucont_y,dx,dy);


for i = 1:M
    for j = 1:N
        
    x = (i-2) * dx;
    y = (j-2) * dy;

    if (i==1 || i == M || j ==1 ||j==N)
      P_Div(i,j) = -( cos(2*x) + cos(2*y))* ( exp(-4*dt) - 1)/4 ;
    end
    end
end

Summation = 0;

for i = 1:M
    for j = 1:N
        P_RHS((i-1)*M +j) = P_Div(i,j);
        
        Summation = Summation + P_Div(i,j);
    end    
end

P_RHS = P_RHS * 1.5 / dt;

Summation = Summation * 1.5 / dt;
% Check summation of RHS

fprintf(' Summation of RHS = %8.7f \n', Summation);

