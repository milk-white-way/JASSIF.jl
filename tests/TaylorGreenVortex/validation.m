close all; clearvars; clc;

USE_VALIDATION = 1;

control_file_path = '<your-path-to-repo>/JASSIF.jl/tests/TaylorGreenVortex/control.m';
working_dir_path  = '<your-path-to-repo>/JASSIF.jl/tests/TaylorGreenVortex';

run(control_file_path);

if USE_VALIDATION
    %% First is load the result you want to validate
    filename = 'chk_<steptime>.mat';
    load(strcat(working_dir_path, '/Checkpoints/', filename));

    %% After this step, the below is list of noticeble loaded to the workspace
    %       + PhysDom : a struct containing the Cartesian velocity and pressure fields
    %       + FluxDom : a struct containing the fluxes
    %       + dx : x grid spacing
    %       + dy : y grid spacing
    %       + t : checkpoint time
    scale_factor = 1;

    %% Then, let us calculate the exact solution    
    % Here we want to compare two solutions on a line from (0, 0) to (1, 1)
    % Since this is analytical, I will use line representation
    M_num = L/dx; dx_num = dx;
    N_num = L/dy; dy_num = dy;
    Numel.x = [];
    Numel.y = [];
    for ii = 1:M_num
        for ji = 1:N_num
            x = (ii - 1 + 0.5)*dx_num;
            y = (ji - 1 + 0.5)*dy_num;
            if x == y
                Numel.x = [Numel.x; x];
                Numel.y = [Numel.y; y];
            end
        end
    end

    M = M_num*scale_factor;
    N = N_num*scale_factor;

    dx = L/M;
    dy = L/N;

    Exact.x = [];
    Exact.y = [];
    Exact.U = [];
    Exact.V = [];
    Exact.P = [];

    for ii = 1:M
        for ji = 1:N
            x = (ii - 1 + 0.5)*dx;
            y = (ji - 1 + 0.5)*dy;

            x_num = (ii - 1 + 0.5)*dx_num;
            y_num = (ji - 1 + 0.5)*dy_num;
            
            if x == y 
                %Uexact_x(ji, ii) =   U*( sin(2*pi*x)*cos(2*pi*y) )*exp(-8*pi^2*t);
                %Uexact_y(ji, ii) =  -U*( cos(2*pi*x)*sin(2*pi*y) )*exp(-8*pi^2*t);
                %Pressure(ji, ii) = -0.25*( cos(4*pi*x) + cos(4*pi*y) )*exp(-16*pi^2*t);

                Uexact_x =   U*( sin(2*pi*x)*cos(2*pi*y) )*exp(-8*pi^2*t);
                Uexact_y =  -U*( cos(2*pi*x)*sin(2*pi*y) )*exp(-8*pi^2*t);
                Pressure = 0.25*U*U*( cos(4*pi*x) + cos(4*pi*y) )*( exp(-8*pi^2*t) * exp(-8*pi^2*t) );

                Exact.x = [Exact.x; x];
                Exact.y = [Exact.y; y];
                Exact.U = [Exact.U; Uexact_x];
                Exact.V = [Exact.V; Uexact_y];
                Exact.P = [Exact.P; Pressure];

            end
        end
    end

    %Exact.U = Uexact_x;
    %Exact.V = Uexact_y;
    %Exact.P = Pressure;
    [Exact.X, Exact.Y] = meshgrid(Exact.x, Exact.y);

    %% Next, extract the numerical solution to our validation line
    % We will use the same line as the exact solution
    Numel.U = diag(PhysDom.Ucat_x);
    Numel.V = diag(PhysDom.Ucat_y);
    Numel.P = diag(PhysDom.Pressure);

    % Plot the exact solution
    figure();
    subplot(1, 3, 1)
        plot(Exact.x, Exact.U, 'k', 'LineWidth', 2);
        hold on;
        plot(Numel.x, Numel.U, 'r');%, ...
        %    'MarkerSize', 8, ...
        %    'MarkerFaceColor', 'red', ...
        %    'MarkerEdgeColor', 'black');
        xlabel('x = y');
        ylabel('U');
        legend('X-Velocity: Exact', ' X-Velocity: Current', 'Location', 'southoutside');
    subplot(1, 3, 2)
        plot(Exact.y, Exact.V, 'k', 'LineWidth', 2);
        hold on;
        plot(Numel.y, Numel.V, 'r');%, ...
        %    'MarkerSize', 8, ...
        %    'MarkerFaceColor', 'green', ...
        %    'MarkerEdgeColor', 'black');
        xlabel('x = y');
        ylabel('V');
        legend('Y-Velocity: Exact', 'Y-Velocity: Current', 'Location', 'southoutside');
    subplot(1, 3, 3)
        plot(Exact.x, Exact.P, 'k', 'LineWidth', 2);
        hold on;
        plot(Numel.x, Numel.P, 'r'); %, ...
        %    'MarkerSize', 8, ...
        %    'MarkerFaceColor', 'cyan', ...
        %    'MarkerEdgeColor', 'black');
        xlabel('x = y');
        ylabel('P');
        legend('Pressure: Exact', 'Pressure: Current', 'Location', 'southoutside');
      
     sgtitle(sprintf('2D Taylor-Green Vortex Solutions at t = %i s', t));
     fontsize(gcf, scale=2.0);

end