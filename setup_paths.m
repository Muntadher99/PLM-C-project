function setup_paths()
%SETUP_PATHS  Add every code folder of the package to the MATLAB path.
%
%   Run once per MATLAB session from any folder:
%       run('<package root>/setup_paths.m')   or   setup_paths
%
%   The results/ and reproduced/ folders hold data only and are not added.

root = fileparts(mfilename('fullpath'));
folders = { ...
    '1_proposed_method_PLM-C'
    '2_baselines'
    '3_benchmark_protocol'
    '4_ablations'
    '5_angle_count'
    '6_compliance_witnesses'
    '7_device_level_validation'
    fullfile('7_device_level_validation', 'ideal_switch_check')
    '8_loss_analysis'
    '9_sensitivity'};
addpath(root);
for i = 1:numel(folders)
    addpath(fullfile(root, folders{i}));
end
fprintf('PLM-C package: %d code folders added to the path (root %s)\n', numel(folders), root);
end
