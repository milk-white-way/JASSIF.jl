function index = glidx(i,j,M,N)
index = (i-1) *N +j;

% Check the validity of the result 
% Return 0 for null
if (i < 1 || i > M ||j <1 ||j>N) 
    index = 0;
end