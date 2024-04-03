close all; clearvars; clc;

USE_PRODUTION = 1;

control_file_path = '<your-path-to-repo>/JASSIF.jl/tests/TaylorGreenVortex/control.m';
working_dir_path  = '<your-path-to-repo>/JASSIF.jl/tests/TaylorGreenVortex';

run(control_file_path);

%% Here I recommend absolute paths
if USE_PRODUTION
    solver_path = '<your-path-to-repo>/JASSIF.jl/src/mat';
    solver_exec = strcat( solver_path, '/Main(control_file_path, working_dir_path)' );
    run(solver_exec);
end
