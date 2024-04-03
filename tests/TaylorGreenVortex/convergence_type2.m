close all; clearvars; clc;
%% Physical parameters
Re = 1;
L = 1;
U = 1;

%% Convergence
USE_CONVERGENCE = 1;

working_dir_path  = '<your-path-to-repo>/JASSIF.jl/tests/TaylorGreenVortex';

if USE_CONVERGENCE
    %% First is load the result you want to validate
    filename = 'chk_<steptime>.mat';
    load(strcat(working_dir_path, '/Checkpoints/', filename));

    [M, N] = size(PhysDom.Pressure);

    %% After this step, the below is list of noticeble loaded to the workspace
    %       + PhysDom : a struct containing the Cartesian velocity and pressure fields
    %       + FluxDom : a struct containing the fluxes
    %       + dx : x grid spacing
    %       + dy : y grid spacing
    %       + t : checkpoint time

    %% Then, let us calculate the exact solution
    %x_vec = [];
    %y_vec = [];

    Numel.U = PhysDom.Ucat_x;
    Numel.V = PhysDom.Ucat_y;
    Numel.P = PhysDom.Pressure;

    Exact.U = zeros(M, N);
    Exact.V = zeros(M, N);
    Exact.P = zeros(M, N);

    for ii = 1:M
        for ji = 1:N
            x = (ii - 1 + 0.5)*dx;
            y = (ji - 1 + 0.5)*dy;

            Exact.U(ii, ji) =  U*( sin(2*pi*x)*cos(2*pi*y) )*exp(-8*pi^2*t);
            Exact.V(ii, ji) =  -U*( cos(2*pi*x)*sin(2*pi*y) )*exp(-8*pi^2*t);
            Exact.P(ii, ji) = 0.25*U*U*( cos(4*pi*x) + cos(4*pi*y) )*( exp(-8*pi^2*t) * exp(-8*pi^2*t) );

            %x_vec = [x_vec; x];
            %y_vec = [y_vec; y];
        end
    end

    %% Calculate the error matrix
    Error.mat_U = abs(Exact.U - Numel.U);
    Error.mat_V = abs(Exact.V - Numel.V);
    Error.mat_P = abs(Exact.P - Numel.P);

    %% Turn the error matrix into a vector
    Error.vec_U = reshape(Error.mat_U, [M*N, 1]);
    Error.vec_V = reshape(Error.mat_V, [M*N, 1]);
    Error.vec_P = reshape(Error.mat_P, [M*N, 1]);

    l2norm_error_ucat_x = norm(Error.vec_U, 2)/(sqrt(M*N));
    l2norm_error_ucat_y = norm(Error.vec_V, 2)/(sqrt(M*N));
    l2norm_error_pressure = norm(Error.vec_P, 2)/(sqrt(M*N));

    fprintf('\nBENCHMARKING | \tL2 norm of error in ucat_x: %.3i\n', l2norm_error_ucat_x);
    fprintf('BENCHMARKING | \tL2 norm of error in ucat_y: %.3i\n', l2norm_error_ucat_y);
    fprintf('BENCHMARKING | \tL2 norm of error in pressure: %.3i\n', l2norm_error_pressure);

end