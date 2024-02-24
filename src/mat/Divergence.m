function [Div] = Divergence(Ucont_x, Ucont_y,dx ,dy)
% This function does calculate the divergence at every grid node in the
% Computational domain

M = length(Ucont_x(:,1));
N = length(Ucont_x(1,:));

for i = 2:M-1
    for j=2 :N-1
        Div(i,j) = (Ucont_x(i,j) - Ucont_x(i-1,j)) / dx + (Ucont_y(i,j) - Ucont_y(i,j-1)) /dy;        
    end
end

% Zero out boundaries
for i=1:M
    Div(i,1) = 0;
    Div(i,N) = 0;
end

for j = 1:N
    Div(1,j) = 0;
    Div(M,j) = 0;
end