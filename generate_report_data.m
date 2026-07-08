%% generate_report_data.m
% Generate figures, numerical results, and LaTeX tables for the
% ELEC9123 Design Task B report.
%
% This script runs the design task file z5615351_YU_DTB_2026.m, saves the
% resulting figures as PDFs, stores the numerical arrays, and writes the
% percentage-deviation tables in LaTeX format to report/tables.tex.

clc;
close all;

designScript = 'z5615351_YU_DTB_2026.m';
if ~exist(designScript, 'file')
    error('Design task script not found: %s', designScript);
end

%% Run the design task script and capture console output
fprintf('Running %s ...\n', designScript);
consoleOutput = evalc('run(designScript)');

% The design script clears the workspace, so define output paths again here.
reportDir = fullfile(pwd, 'report');
figDir    = fullfile(reportDir, 'figs');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

fid = fopen(fullfile(reportDir, 'console_output.txt'), 'w');
fwrite(fid, consoleOutput);
fclose(fid);

%% Save figures as vector PDFs
figs = findall(0, 'Type', 'figure');
[~, order] = sort([figs.Number]);
figs = figs(order);

for i = 1:numel(figs)
    fig = figs(i);
    figFile = fullfile(figDir, sprintf('figure_%d.pdf', i));
    % Use exportgraphics for vector PDF with automatic page fitting.
    exportgraphics(fig, figFile, 'ContentType', 'vector', 'Resolution', 300);
    fprintf('Saved figure %d to %s\n', fig.Number, figFile);
end

%% Save numerical results
results.SNR_dB        = SNR_dB;
results.BER_AWGN_sim  = BER_AWGN_sim;
results.BER_AWGN_ana  = BER_AWGN_ana;
results.BER_Lap_sim   = BER_Lap_sim;
results.BER_Lap_ana   = BER_Lap_ana;
results.BER_Ray_sim   = BER_Ray_sim;
results.BER_Ray_ana   = BER_Ray_ana;
results.BER_Ric_sim   = BER_Ric_sim;
results.BER_Ric_ana   = BER_Ric_ana;
results.BER_RayRand_sim = BER_RayRand_sim;
results.BER_RayRand_ana = BER_RayRand_ana;
results.BER_RicRand_sim = BER_RicRand_sim;
results.BER_RicRand_ana = BER_RicRand_ana;
results.Pout_Ray_sim  = Pout_Ray_sim;
results.Pout_Ray_ana  = Pout_Ray_ana;
results.Pout_Ric_sim  = Pout_Ric_sim;
results.Pout_Ric_ana  = Pout_Ric_ana;
results.Pout_RayRand_sim = Pout_RayRand_sim;
results.Pout_RayRand_ana = Pout_RayRand_ana;
results.Pout_RicRand_sim = Pout_RicRand_sim;
results.Pout_RicRand_ana = Pout_RicRand_ana;
results.max_cdf_error = max_cdf_error;
results.N_samples     = N_samples;

save(fullfile(reportDir, 'results.mat'), '-struct', 'results');

%% Write LaTeX tables
pct_err = @(sim, ana) abs(sim - ana) ./ max(ana, eps) * 100;

fid = fopen(fullfile(reportDir, 'tables.tex'), 'w');

% BER percentage-deviation table
fprintf(fid, '\\begin{table}[htbp]\n\\centering\n');
fprintf(fid, '\\caption{BER Percentage Deviation (Simulated vs Analytical)}\n');
fprintf(fid, '\\label{tab:ber_dev}\n');
fprintf(fid, '\\begin{tabular}{rcccccc}\n\\toprule\n');
fprintf(fid, 'SNR (dB) & AWGN & Laplace & Rayleigh & Rician & RayRand & RicRand \\\\\n\\midrule\n');
for idx = 1:length(SNR_dB)
    fprintf(fid, '%d & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f \\\\\n', ...
        SNR_dB(idx), ...
        pct_err(BER_AWGN_sim(idx), BER_AWGN_ana(idx)), ...
        pct_err(BER_Lap_sim(idx),  BER_Lap_ana(idx)), ...
        pct_err(BER_Ray_sim(idx),  BER_Ray_ana(idx)), ...
        pct_err(BER_Ric_sim(idx),  BER_Ric_ana(idx)), ...
        pct_err(BER_RayRand_sim(idx), BER_RayRand_ana(idx)), ...
        pct_err(BER_RicRand_sim(idx), BER_RicRand_ana(idx)));
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n\n');

% Outage percentage-deviation table
fprintf(fid, '\\begin{table}[htbp]\n\\centering\n');
fprintf(fid, '\\caption{Outage Probability Percentage Deviation (Simulated vs Analytical)}\n');
fprintf(fid, '\\label{tab:outage_dev}\n');
fprintf(fid, '\\begin{tabular}{rcccc}\n\\toprule\n');
fprintf(fid, 'SNR (dB) & Rayleigh & Rician & RayRand & RicRand \\\\\n\\midrule\n');
for idx = 1:length(SNR_dB)
    fprintf(fid, '%d & %.3f & %.3f & %.3f & %.3f \\\\\n', ...
        SNR_dB(idx), ...
        pct_err(Pout_Ray_sim(idx), Pout_Ray_ana(idx)), ...
        pct_err(Pout_Ric_sim(idx), Pout_Ric_ana(idx)), ...
        pct_err(Pout_RayRand_sim(idx), Pout_RayRand_ana(idx)), ...
        pct_err(Pout_RicRand_sim(idx), Pout_RicRand_ana(idx)));
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n\n');

% Optional sample-size convergence demonstration (reproducible)
demo_SNR = 10;
demo_gb  = 10^(demo_SNR/10);
demo_ana = 0.5 * erfc(sqrt(demo_gb));
rng(42);
test_bits = 2*(randi([0 1], 5e6, 1)) - 1;
test_n    = sqrt(1/2) * (randn(5e6, 1) + 1j*randn(5e6, 1));
test_r    = sqrt(demo_gb) * test_bits + test_n;
Ns_test   = [1e5 5e5 1e6 5e6];

fprintf(fid, '\\begin{table}[htbp]\n\\centering\n');
fprintf(fid, '\\caption{AWGN BER Convergence with Sample Size at $\\gamma_b=%d$~dB (Analytical BER = %.3e)}\n', demo_SNR, demo_ana);
fprintf(fid, '\\label{tab:sample_conv}\n');
fprintf(fid, '\\begin{tabular}{rcc}\n\\toprule\n');
fprintf(fid, '$N$ & Simulated BER & Deviation (\\%%) \\\\\n\\midrule\n');
for Nt = Ns_test
    ber_test = mean(real(test_r(1:Nt)) .* test_bits(1:Nt) < 0);
    err_test = abs(ber_test - demo_ana) / demo_ana * 100;
    fprintf(fid, '%.0f & %.3e & %.3f \\\\\n', Nt, ber_test, err_test);
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n');

fclose(fid);

fprintf('Report data generated in %s\n', reportDir);
