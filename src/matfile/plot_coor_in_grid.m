%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: G:\RunningCode\Original_Project\coor_ucat.csv
%
% Auto-generated by MATLAB on 12-Feb-2024 10:58:54

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["x", "y"];
opts.VariableTypes = ["double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
coord_ucat = readtable("G:\RunningCode\Original_Project\coor_ucat.csv", opts);
coord_ucont_x = readtable("G:\RunningCode\Original_Project\coor_ucont_x.csv", opts);
coord_ucont_y = readtable("G:\RunningCode\Original_Project\coor_ucont_y.csv", opts);

%% Convert to output type
x_ucat = coord_ucat.x;
y_ucat = coord_ucat.y;

x_ucont_x = coord_ucont_x.x;
y_ucont_x = coord_ucont_x.y;

x_ucont_y = coord_ucont_y.x;
y_ucont_y = coord_ucont_y.y;

%text(x_ucat, y_ucat, 'c', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 10);
plot(x_ucat, y_ucat, 'p');
hold on;
plot(x_ucont_x, y_ucont_x, 'x');
plot(x_ucont_y, y_ucont_y, 'o');
grid on;

%% Clear temporary variables
clear opts tbl