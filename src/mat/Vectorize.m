function [Output] = Vectorize(Ucont)

MM = length(Ucont(:,1));
NN = length(Ucont(1,:));
% This function Vectorize(U(i,j)) = U(index)
Output = zeros(MM*NN, 1);

for ii = 1:MM
    for jj= 1:NN
        glb_idx = glidx(MM, NN, ii, jj);
        Output(glb_idx) = Ucont(ii, jj);
    end
end