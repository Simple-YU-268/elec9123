% z5615351_YU_DTB_DTD_2026.m
% Top-level driver for ELEC9123 Design Task D (Tiered Learning Taxonomy)
% Run the appropriate sub-model by uncommenting the relevant section.
% Each sub-folder contains its own Main/Analysis/Simulation triplet.

% Make the sub-model functions callable regardless of the current folder
this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'M_M_1'), ...
        fullfile(this_dir, 'M_M_1_K'), ...
        fullfile(this_dir, 'M_G_1'), ...
        fullfile(this_dir, 'M_M_m'));

%% Level 2: M/M/1
mm1_main();

%% Level 3: M/M/1/K
mm1k_main();

%% Level 4: M/G/1
mg1_main();

%% Level 5: Two-class M/M/m
mmm_main();
