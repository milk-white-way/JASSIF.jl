function index = glidx(M, N, ii, jj)
    % Tam N: Passed
    index = (jj-1)*M + ii;

    if (ii < 1 || ii > M || jj < 1 ||jj > N) 
        index = 0;
    end
end