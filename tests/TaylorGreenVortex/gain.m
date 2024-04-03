clearvars; close all; clc;

%A = load('L2NormError_16x16_4E-3_50.mat');
A = load('L2NormError_16x16_2E-3_100.mat');
B = load('L2NormError_16x16_1E-3_200.mat');

fprintf('BENCHMARKING | \tAccuracy Gain in Velocity: %.3i\n', A.l2norm_error_ucat_x/B.l2norm_error_ucat_x);
fprintf('BENCHMARKING | \tAccuracy Gain in Pressure: %.3i\n', A.l2norm_error_pressure/B.l2norm_error_pressure);