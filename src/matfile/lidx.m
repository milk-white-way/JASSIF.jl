function [i j] = lidx(IM,M,N)

i = floor(IM / N);
j = IM  - (N * i);
if (j ~= 0)
i = i + 1;
else
    j = N;
end

