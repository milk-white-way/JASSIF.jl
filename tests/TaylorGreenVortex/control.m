%% Solver's parameters
M = 16;
N = 16;
dt = 1E-5;
MAXTIME = 1000;

%% Physical parameters
Re = 1;
L = 1;
U = 1;

%% Behaviour parameters
checkpoint_freq = 100;

%% FLAGS
ENABLE_CALCULATION  = 1;
ENABLE_BC_PERIODIC  = 1;

%ENABLE_CHECKPOINT   = 1; Future work

ENABLE_VISUAL_GRID  = 0;
ENABLE_VISUAL_PLOT  = 0;

ENABLE_AMRESSIF     = 0;
ENABLE_DEBUGGING    = 0;

