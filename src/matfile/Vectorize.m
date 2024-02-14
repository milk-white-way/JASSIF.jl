function [Output] = Vectorize(Ucont)

M = length(Ucont(:,1));
N = length(Ucont(1,:));
% This function Vectorize(U(i,j)) = U(index)

for i = 1:M 
    for j= 1:N
        index = glidx(i,j,M,N);
        Output(index) = Ucont(i,j);
    end
end

% Transform to column vector
Output = Output';