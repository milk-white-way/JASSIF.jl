function [ii, jj] = lidx(M, glb_idx)

    jj = floor(glb_idx / M);
    ii = glb_idx - (M * jj);

    if (ii ~= 0)
        jj = jj + 1;
    else
        ii = M;
    end

end