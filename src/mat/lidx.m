function [ii, jj] = lidx(index, M4)

    jj = floor(index / M4);
    ii = index - (M4 * jj);

    if (ii ~= 0)
        jj = jj + 1;
    else
        ii = M4;
    end

end