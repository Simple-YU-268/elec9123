%% generate_report_data.m
% 为 ELEC9123 设计任务 B 的实验报告生成图片、数值结果和 LaTeX 表格。
%
% 本脚本会：
%   1) 调用主设计任务脚本 z5615351_YU_DTB_2026.m 并捕获命令行输出；
%   2) 将产生的所有 Figure 保存为矢量 PDF；
%   3) 将所有数值结果保存为 results.mat；
%   4) 生成 LaTeX 格式的百分比偏差表格，写入 report/tables.tex。

clc;        % 清空命令行窗口
clear;      % 清除工作区变量
close all;  % 关闭所有 Figure 窗口

% 主设计任务脚本的文件名
designScript = 'z5615351_YU_DTB_2026.m';
% 如果找不到该脚本，则抛出错误并停止运行
if ~exist(designScript, 'file')
    error('找不到设计任务脚本: %s', designScript);
end

%% ===================== 运行设计任务脚本并捕获控制台输出 =====================
fprintf('正在运行 %s ...\n', designScript);
% evalc 会执行 MATLAB 代码并将所有打印到命令行的内容捕获为字符串
consoleOutput = evalc('run(designScript)');

% 注意：主脚本开头执行了 clear，因此这里需要重新定义输出路径
reportDir = fullfile(pwd, 'report');       % 报告目录：当前工作目录/report
figDir    = fullfile(reportDir, 'figs');   % 图片子目录
if ~exist(figDir, 'dir')
    mkdir(figDir);  % 如果目录不存在则创建
end

% 将捕获到的控制台输出写入文本文件，便于报告中引用
fid = fopen(fullfile(reportDir, 'console_output.txt'), 'w');
fwrite(fid, consoleOutput);
fclose(fid);

%% ===================== 将 Figure 保存为矢量 PDF =====================
% findall(0, 'Type', 'figure') 找到所有已创建的 Figure 句柄
figs = findall(0, 'Type', 'figure');
% 按 Figure 编号排序，确保保存顺序与生成顺序一致
[~, order] = sort([figs.Number]);
figs = figs(order);

for i = 1:numel(figs)
    fig = figs(i);
    figFile = fullfile(figDir, sprintf('figure_%d.pdf', i));
    % exportgraphics 用于高质量导出，'ContentType', 'vector' 表示矢量 PDF
    exportgraphics(fig, figFile, 'ContentType', 'vector', 'Resolution', 300);
    fprintf('已将 Figure %d 保存至 %s\n', fig.Number, figFile);
end

%% ===================== 保存数值结果到 MAT 文件 =====================
% 将工作区中由主脚本生成的关键变量打包到 results 结构体
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

% 以结构体形式保存，方便后续用 load('results.mat') 直接读取各字段
save(fullfile(reportDir, 'results.mat'), '-struct', 'results');

%% ===================== 写入 LaTeX 表格 =====================
% 百分比偏差计算函数：防止除以 0，分母取 max(理论值, eps)
pct_err = @(sim, ana) abs(sim - ana) ./ max(ana, eps) * 100;

% 打开 tables.tex 文件用于写入
fid = fopen(fullfile(reportDir, 'tables.tex'), 'w');

% ---------- 表 1：BER 百分比偏差 ----------
fprintf(fid, '\\begin{table}[htbp]\n\\centering\n');
fprintf(fid, '\\caption{BER 仿真值与理论值的百分比偏差}\n');
fprintf(fid, '\\label{tab:ber_dev}\n');
fprintf(fid, '\\begin{tabular}{rcccccc}\n\\toprule\n');
fprintf(fid, 'SNR (dB) & AWGN & Laplace & Rayleigh & Rician & RayRand & RicRand \\\n\\midrule\n');
for idx = 1:length(SNR_dB)
    fprintf(fid, '%d & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f \\\n', ...
        SNR_dB(idx), ...
        pct_err(BER_AWGN_sim(idx), BER_AWGN_ana(idx)), ...
        pct_err(BER_Lap_sim(idx),  BER_Lap_ana(idx)), ...
        pct_err(BER_Ray_sim(idx),  BER_Ray_ana(idx)), ...
        pct_err(BER_Ric_sim(idx),  BER_Ric_ana(idx)), ...
        pct_err(BER_RayRand_sim(idx), BER_RayRand_ana(idx)), ...
        pct_err(BER_RicRand_sim(idx), BER_RicRand_ana(idx)));
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n\n');

% ---------- 表 2：Outage 概率百分比偏差 ----------
fprintf(fid, '\\begin{table}[htbp]\n\\centering\n');
fprintf(fid, '\\caption{Outage 概率仿真值与理论值的百分比偏差}\n');
fprintf(fid, '\\label{tab:outage_dev}\n');
fprintf(fid, '\\begin{tabular}{rcccc}\n\\toprule\n');
fprintf(fid, 'SNR (dB) & Rayleigh & Rician & RayRand & RicRand \\\n\\midrule\n');
for idx = 1:length(SNR_dB)
    fprintf(fid, '%d & %.3f & %.3f & %.3f & %.3f \\\n', ...
        SNR_dB(idx), ...
        pct_err(Pout_Ray_sim(idx), Pout_Ray_ana(idx)), ...
        pct_err(Pout_Ric_sim(idx), Pout_Ric_ana(idx)), ...
        pct_err(Pout_RayRand_sim(idx), Pout_RayRand_ana(idx)), ...
        pct_err(Pout_RicRand_sim(idx), Pout_RicRand_ana(idx)));
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n\n');

% ---------- 表 3：样本数收敛演示（可复现版本） ----------
% 固定随机数种子，使该表格结果在每次运行 generate_report_data 时保持一致
demo_SNR = 10;
demo_gb  = 10^(demo_SNR/10);
demo_ana = 0.5 * erfc(sqrt(demo_gb));
rng(42);  % 设置随机种子
test_bits = 2*(randi([0 1], 5e6, 1)) - 1;
test_n    = sqrt(1/2) * (randn(5e6, 1) + 1j*randn(5e6, 1));
test_r    = sqrt(demo_gb) * test_bits + test_n;
Ns_test   = [1e5 5e5 1e6 5e6];

fprintf(fid, '\\begin{table}[htbp]\n\\centering\n');
fprintf(fid, '\\caption{AWGN BER 随样本数收敛（$\\gamma_b=%d$~dB，理论 BER = %.3e）}\n', demo_SNR, demo_ana);
fprintf(fid, '\\label{tab:sample_conv}\n');
fprintf(fid, '\\begin{tabular}{rcc}\n\\toprule\n');
fprintf(fid, '$N$ & 仿真 BER & 偏差 (\\%%) \\\n\\midrule\n');
for Nt = Ns_test
    ber_test = mean(real(test_r(1:Nt)) .* test_bits(1:Nt) < 0);
    err_test = abs(ber_test - demo_ana) / demo_ana * 100;
    fprintf(fid, '%.0f & %.3e & %.3f \\\n', Nt, ber_test, err_test);
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n');

% 关闭文件
fclose(fid);

fprintf('报告数据已生成在目录：%s\n', reportDir);
