function [Ucat_comp_x, Ucat_comp_y, Pressure_comp] = ...
    TAM_ensemble_comp(PhysDom, HaloDom, M, N, Nghost)

    %% Ghost variables' dimension in computational non-staggered grid
    M2 = M + (2*Nghost);
    N2 = N + (2*Nghost);

    Ucat_comp_x = zeros(M2,N2);
    Ucat_comp_y = zeros(M2,N2);
    Pressure_comp = zeros(M2,N2);

    for ii = 1:M2
        for jj = 1:N2
            if (ii <= Nghost) || (ii > (M + Nghost)) || (jj <= Nghost) || (jj > (N + Nghost))
                Ucat_comp_x(jj, ii) = HaloDom.Ucat_x(jj, ii);
                Ucat_comp_y(jj, ii) = HaloDom.Ucat_y(jj, ii);
                Pressure_comp(jj, ii) = HaloDom.Pressure(jj, ii);
            else
                Ucat_comp_x(jj, ii) = PhysDom.Ucat_x(jj-Nghost, ii-Nghost);
                Ucat_comp_y(jj, ii) = PhysDom.Ucat_y(jj-Nghost, ii-Nghost);
                Pressure_comp(jj, ii) = PhysDom.Pressure(jj-Nghost, ii-Nghost);
            end
        end
    end

end