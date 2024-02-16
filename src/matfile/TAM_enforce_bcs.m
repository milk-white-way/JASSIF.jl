fprintf('\nINFO: \t Enforcing boundary conditions... \n'); 

if ENABLE_BC_PERIODIC
    for ii = 2:M2-1
            
        jj = 1; 
        Ubcs_x(ii, jj) = Ucat_phy_x(ii-1, N-1);
        Ubcs_y(ii, jj) = Ucat_phy_y(ii-1, N-1);
        Pbcs(ii, jj) = Pressure_phy(ii-1, N-1);
        if ENABLE_DEBUGGING
            fprintf('Ghost node (%d, %d) ', ii, jj);
            fprintf('mirrors inner node (%d, %d)\n', ii-1, N-1);
        end

        jj = N2; 
        Ubcs_x(ii, jj) = Ucat_phy_x(ii-1, 2);
        Ubcs_y(ii, jj) = Ucat_phy_y(ii-1, 2);
        Pbcs(ii, jj) = Pressure_phy(ii-1, 2);
        if ENABLE_DEBUGGING
            fprintf('Ghost node (%d, %d) ', ii, jj);
            fprintf('mirrors inner node (%d, %d)\n', ii-1, 2);
        end

    end

    Ucat_cal_x(:,1) = Ubcs_x(:,1);
    Ucat_cal_y(:,1) = Ubcs_y(:,1);
    Pressure_cal(:,1) = Pbcs(:,1);

    Ucat_cal_x(:,N2) = Ubcs_x(:,N2);
    Ucat_cal_y(:,N2) = Ubcs_y(:,N2);
    Pressure_cal(:,N2) = Pbcs(:,N2);

    for jj = 2:N2-1
            
        ii = 1; 
        Ubcs_x(ii, jj) = Ucat_phy_x(M-1, jj-1);
        Ubcs_y(ii, jj) = Ucat_phy_y(M-1, jj-1);
        Pbcs(ii, jj) = Pressure_phy(M-1, jj-1);
        if ENABLE_DEBUGGING
            fprintf('Ghost node (%d, %d) ', ii, jj);
            fprintf('mirrors inner node (%d, %d)\n', M-1, jj-1);
        end

        ii = M2; 
        Ubcs_x(ii, jj) = Ucat_phy_x(2, jj-1);
        Ubcs_y(ii, jj) = Ucat_phy_y(2, jj-1);
        Pbcs(ii, jj) = Pressure_phy(2, jj-1);
        if ENABLE_DEBUGGING
            fprintf('Ghost node (%d, %d) ', ii, jj);
            fprintf('mirrors inner node (%d, %d)\n', 2, jj-1);
        end

    end

    Ucat_cal_x(1,:) = Ubcs_x(1,:);
    Ucat_cal_y(1,:) = Ubcs_y(1,:);
    Pressure_cal(1,:) = Pbcs(1,:);
        
    Ucat_cal_x(M2,:) = Ubcs_x(M2,:);
    Ucat_cal_y(M2,:) = Ubcs_y(M2,:);
    Pressure_cal(M2,:) = Pbcs(M2,:);

    fprintf('INFO: \t Periodic boundary conditions are enabled!\n');
end

% Return nan to the corners
Ucat_cal_x(1,1) = nan;
Ucat_cal_y(1,1) = nan;
Pressure_cal(1,1) = nan;

Ucat_cal_x(1,N2) = nan;
Ucat_cal_y(1,N2) = nan;
Pressure_cal(1,N2) = nan;

Ucat_cal_x(M2,1) = nan;
Ucat_cal_y(M2,1) = nan;
Pressure_cal(M2,1) = nan;
    
Ucat_cal_x(M2,N2) = nan;
Ucat_cal_y(M2,N2) = nan;
Pressure_cal(M2,N2) = nan;