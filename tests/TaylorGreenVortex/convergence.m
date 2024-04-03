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
    filename1 = 'chk_ncell256_dt2E-6_max100.mat';
    Coarse = load(strcat(working_dir_path, '/Checkpoints/', filename1));

    filename2 = 'chk_ncell256_dt1E-6_max200.mat';
    Fine = load(strcat(working_dir_path, '/Checkpoints/', filename2));

    [M, N] = size(Fine.PhysDom.Pressure);

    %% After this step, the below is list of noticeble loaded to the workspace
    %       + PhysDom : a struct containing the Cartesian velocity and pressure fields
    %       + FluxDom : a struct containing the fluxes
    %       + dx : x grid spacing
    %       + dy : y grid spacing
    %       + t : checkpoint time

    %% Then, let us calculate the exact solution
    %x_vec = [];
    %y_vec = [];

    %% Calculate the error matrix
    Error.mat_U = abs(Coarse.PhysDom.Ucat_x - Fine.PhysDom.Ucat_x);
    Error.mat_V = abs(Coarse.PhysDom.Ucat_y - Fine.PhysDom.Ucat_y);
    Error.mat_P = abs(Coarse.PhysDom.Pressure - Fine.PhysDom.Pressure);

    %% Turn the error matrix into a vector
    Error.vec_U = reshape(Error.mat_U, [M*N, 1]);
    Error.vec_V = reshape(Error.mat_V, [M*N, 1]);
    Error.vec_P = reshape(Error.mat_P, [M*N, 1]);

    l2norm_error_ucat_x = norm(Error.vec_U, 2)/(sqrt(M*N));
    l2norm_error_ucat_y = norm(Error.vec_V, 2)/(sqrt(M*N));
    l2norm_error_pressure = norm(Error.vec_P, 2)/(sqrt(M*N));

    fprintf('\nBENCHMARKING | \tL2 norm of error in ucat_x: %.5i\n', l2norm_error_ucat_x);
    fprintf('BENCHMARKING | \tL2 norm of error in ucat_y: %.5i\n', l2norm_error_ucat_y);
    fprintf('BENCHMARKING | \tL2 norm of error in pressure: %.5i\n', l2norm_error_pressure);

end