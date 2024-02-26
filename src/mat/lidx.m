function [jj, ii] = lidx(index, M)

    ii = floor(index / M);
    jj = index - (M * ii);

    if (jj ~= 0)
        ii = ii + 1;
    else
        ii = M;
    end

end

