%% ELEC9123: 设计任务 B（无线通信）
% 通过蒙特卡洛仿真验证 BPSK 系统在不同无线信道模型下的性能指标
%
% 作者: YU
% 学号: z5615351
% 文件名: z5615351_YU_DTB_2026.m
%
% 本脚本实现了六种无线信道模型下 BPSK 系统的蒙特卡洛仿真，
% 并将误码率（BER）和 outage 概率与理论解析式进行对比验证。
% 六种信道模型包括：
%   1) 仅 AWGN 信道
%   2) 仅拉普拉斯（Laplacian）脉冲噪声信道
%   3) 瑞利（Rayleigh）衰落 + AWGN
%   4) 莱斯（Rician）衰落 + AWGN
%   5) 随机部署用户 + 瑞利衰落 + AWGN
%   6) 随机部署用户 + 莱斯衰落 + AWGN

%% ===================== 环境清理 =====================
clc;        % 清空命令行窗口，保持输出整洁
clear;      % 清除工作区所有变量，避免历史变量干扰
close all;  % 关闭所有已打开的 Figure 窗口

%% ===================== 系统参数设置 =====================
% 下面所有参数均来自设计任务说明书（Task B）中的默认设定

C = 1.2;                         % 中断门限容量，单位 bps/Hz
                                 % 当瞬时信道容量 log2(1+SNR) < C 时，认为发生 outage
SNR_dB = 0:1:15;                 % 平均信噪比（每比特）的扫描范围，单位 dB，从 0 dB 到 15 dB，步长 1 dB
SNR_lin = 10.^(SNR_dB/10);       % 将 dB 信噪比转换为线性值 γ_b（平均 SNR per bit）
                                 % 公式：γ_b(linear) = 10^(γ_b(dB)/10)
R = 3;                           % 小区（圆形覆盖区域）半径，单位 m
nu = 2.2;                        % 路径损耗指数（path-loss exponent），典型取值 2~4
hd = 1;                          % 直射径（LOS）分量的幅度，取单位幅度
K = 5;                           % 莱斯 K 因子：直射径功率与散射径功率之比
N0 = 1;                          % 噪声功率/方差，这里取归一化值 1
N_samples = 1e6;                 % 每个 SNR 点的蒙特卡洛样本数（题目要求至少约 10^6）
gamma_th = 2^C - 1;              % outage 对应的 SNR 门限
                                 % 由 C = log2(1+γ_th) 可得 γ_th = 2^C - 1

% 打印任务信息
fprintf('ELEC9123 设计任务 B - 蒙特卡洛仿真\n');
fprintf('作者: YU (z5615351)\n');
fprintf('每个 SNR 点的样本数: %d\n\n', N_samples);

%% ===================== 任务 2：验证拉普拉斯随机变量 =====================
% 目标：使用逆 CDF 法生成服从拉普拉斯分布的噪声，并与理论 PDF/CDF 对比。
%
% 拉普拉斯分布 Laplace(0, b) 的方差为 2*b^2。
% 题目给定噪声方差为 N0，因此 2*b^2 = N0，得到 b = sqrt(N0/2)。
% 拉普拉斯分布的 CDF 为：
%   F(x) = 0.5 * exp(sqrt(2/N0)*x),            x <= 0
%   F(x) = 1 - 0.5 * exp(-sqrt(2/N0)*x),       x > 0
% 其逆 CDF（用于由均匀随机变量 u∈[0,1] 生成 x）为：
%   x = -b * sign(u-0.5) * log(1 - 2*|u-0.5|)

N_lap = 1e6;                     % 用于验证拉普拉斯分布的样本数
b_lap = sqrt(N0/2);              % 拉普拉斯分布的尺度参数 b
u_lap = rand(N_lap, 1);          % 生成 [0,1] 均匀分布随机变量
% 逆 CDF 变换：生成服从 Laplace(0, b_lap) 的随机变量
n_lap = -b_lap * sign(u_lap - 0.5) .* log(1 - 2*abs(u_lap - 0.5));

% 解析 PDF 和 CDF 的绘制点
x_lap = linspace(-6*b_lap, 6*b_lap, 1000);
% 理论 PDF：f(x) = (1/(2b)) * exp(-|x|/b)
f_lap_ana = (1/(2*b_lap)) * exp(-abs(x_lap)/b_lap);
% 理论 CDF：分段定义
F_lap_ana = zeros(size(x_lap));
F_lap_ana(x_lap <= 0) = 0.5 * exp(x_lap(x_lap <= 0)/b_lap);
F_lap_ana(x_lap > 0) = 1 - 0.5 * exp(-x_lap(x_lap > 0)/b_lap);

% 仿真经验 CDF：对每个 x 值，统计样本中小于等于 x 的比例
F_lap_sim = zeros(size(x_lap));
for i = 1:length(x_lap)
    F_lap_sim(i) = mean(n_lap <= x_lap(i));
end

% 绘制拉普拉斯分布验证图（PDF 和 CDF 各一个子图）
figure('Name', '拉普拉斯分布验证', 'Position', [100 100 900 400]);

subplot(1, 2, 1);
histogram(n_lap, 100, 'Normalization', 'pdf', 'FaceColor', [0.3 0.6 0.9]);
hold on;
plot(x_lap, f_lap_ana, 'r-', 'LineWidth', 2);
xlabel('x');
ylabel('PDF f_n(x)');
title('拉普拉斯噪声概率密度函数 (PDF)');
legend('仿真直方图', '理论 PDF', 'Location', 'best');
grid on;

subplot(1, 2, 2);
plot(x_lap, F_lap_sim, 'b-', 'LineWidth', 1.5);
hold on;
plot(x_lap, F_lap_ana, 'r--', 'LineWidth', 2);
xlabel('x');
ylabel('CDF F_n(x)');
title('拉普拉斯噪声累积分布函数 (CDF)');
legend('仿真 CDF', '理论 CDF', 'Location', 'best');
grid on;

% 计算仿真 CDF 与理论 CDF 之间的最大绝对误差（类似 Kolmogorov-Smirnov 统计量）
max_cdf_error = max(abs(F_lap_sim - F_lap_ana));
fprintf('任务 2：拉普拉斯 CDF 最大绝对误差 = %.4e\n\n', max_cdf_error);

%% ===================== 预先生成公共随机变量 =====================
% 下面这些随机变量与 SNR 无关，因此在进入 SNR 循环之前一次性生成，
% 可显著提高效率并保证不同 SNR 下使用同一组信道实现，便于对比。

% BPSK 调制符号：等概率取 +1 或 -1
bits = 2*(randi([0 1], N_samples, 1)) - 1;

% 瑞利衰落系数 h_s：复高斯 CN(0,1)，这里除以 sqrt(2) 使 E{|h_s|^2}=1
% |h_s| 服从瑞利分布，|h_s|^2 服从指数分布 Exp(1)
h_rayleigh = (randn(N_samples, 1) + 1j*randn(N_samples, 1)) / sqrt(2);

% 莱斯衰落系数 h：直射径 + 散射径
% h = sqrt(K/(K+1))*hd + sqrt(1/(K+1))*h_s
% 其中 hd 为确定性 LOS 分量，h_s 为 CN(0,1) 散射分量。
% 这样构造满足 E{|h|^2} = K/(K+1)*|hd|^2 + 1/(K+1)*1 = 1
h_los = sqrt(K/(K+1)) * hd;
h_scat = sqrt(1/(K+1)) * (randn(N_samples, 1) + 1j*randn(N_samples, 1)) / sqrt(2);
h_rician = h_los + h_scat;

% 随机用户部署距离 d：在半径为 R 的圆盘内均匀分布
% 圆盘内均匀分布时，距离 d 的 PDF 为 f_d(d) = 2d/R^2, 0 <= d <= R
% 生成方法：d = R * sqrt(rand)，因为 CDF F_d(d) = d^2/R^2，逆 CDF 即 d = R*sqrt(u)
d_user = R * sqrt(rand(N_samples, 1));

% 将衰落与随机距离合并，得到包含路径损耗 d^{-nu/2} 的等效信道
% 等效信道增益 = h / d^(nu/2)，因此接收功率额外乘以 d^{-nu}
h_rayleigh_rand = h_rayleigh ./ (d_user.^(nu/2));
h_rician_rand = h_rician ./ (d_user.^(nu/2));

%% ===================== 任务 3 & 4：BER 仿真 =====================
% 对六种信道模型分别进行蒙特卡洛 BER 仿真，并与理论解析值对比。

% 预分配存储 BER 结果的数组
BER_AWGN_sim = zeros(size(SNR_dB));   % AWGN 信道仿真 BER
BER_AWGN_ana = zeros(size(SNR_dB));   % AWGN 信道理论 BER
BER_Lap_sim  = zeros(size(SNR_dB));   % 拉普拉斯噪声信道仿真 BER
BER_Lap_ana  = zeros(size(SNR_dB));   % 拉普拉斯噪声信道理论 BER
BER_Ray_sim  = zeros(size(SNR_dB));   % 瑞利衰落+AWGN 仿真 BER
BER_Ray_ana  = zeros(size(SNR_dB));   % 瑞利衰落+AWGN 理论 BER
BER_Ric_sim  = zeros(size(SNR_dB));   % 莱斯衰落+AWGN 仿真 BER
BER_Ric_ana  = zeros(size(SNR_dB));   % 莱斯衰落+AWGN 理论 BER
BER_RayRand_sim = zeros(size(SNR_dB));% 随机用户+瑞利衰落 仿真 BER
BER_RayRand_ana = zeros(size(SNR_dB));% 随机用户+瑞利衰落 理论 BER
BER_RicRand_sim = zeros(size(SNR_dB));% 随机用户+莱斯衰落 仿真 BER
BER_RicRand_ana = zeros(size(SNR_dB));% 随机用户+莱斯衰落 理论 BER

% 遍历每个 SNR 点
for idx = 1:length(SNR_dB)
    gb = SNR_lin(idx);              % 当前平均 SNR（线性值）
    Eb = gb * N0;                   % 每比特能量：E_b = γ_b * N_0
    
    %% 1) AWGN 信道
    % 复高斯噪声：实部与虚部各为 N(0, N0/2)，总方差为 N0
    n_awgn = sqrt(N0/2) * (randn(N_samples, 1) + 1j*randn(N_samples, 1));
    % 接收信号：r = sqrt(Eb)*s + n
    r_awgn = sqrt(Eb) * bits + n_awgn;
    % BPSK 判决：real(r)*s < 0 表示判决错误（匹配滤波后取实部）
    BER_AWGN_sim(idx) = mean(real(r_awgn) .* bits < 0);
    % 理论 BER：P_b = 0.5 * erfc(sqrt(γ_b))
    BER_AWGN_ana(idx) = 0.5 * erfc(sqrt(gb));
    
    %% 2) 拉普拉斯脉冲噪声信道
    % 同样使用逆 CDF 法生成拉普拉斯噪声
    u = rand(N_samples, 1);
    n_lap = -b_lap * sign(u - 0.5) .* log(1 - 2*abs(u - 0.5));
    r_lap = sqrt(Eb) * bits + n_lap;
    BER_Lap_sim(idx) = mean(real(r_lap) .* bits < 0);
    % 理论 BER：P_b = 0.5 * exp(-sqrt(2*γ_b))
    BER_Lap_ana(idx) = 0.5 * exp(-sqrt(2*gb));
    
    %% 3) 瑞利衰落 + AWGN
    n_ray = sqrt(N0/2) * (randn(N_samples, 1) + 1j*randn(N_samples, 1));
    r_ray = h_rayleigh * sqrt(Eb) .* bits + n_ray;
    % 等效接收（最大比合并/匹配滤波）：将接收信号乘以信道共轭并取实部
    r_eq_ray = real(conj(h_rayleigh) .* r_ray);
    BER_Ray_sim(idx) = mean(r_eq_ray .* bits < 0);
    % 理论 BER：P_b = 0.5 * (1 - sqrt(γ_b/(1+γ_b)))
    BER_Ray_ana(idx) = 0.5 * (1 - sqrt(gb ./ (1 + gb)));
    
    %% 4) 莱斯衰落 + AWGN
    n_ric = sqrt(N0/2) * (randn(N_samples, 1) + 1j*randn(N_samples, 1));
    r_ric = h_rician * sqrt(Eb) .* bits + n_ric;
    r_eq_ric = real(conj(h_rician) .* r_ric);
    BER_Ric_sim(idx) = mean(r_eq_ric .* bits < 0);
    % 理论 BER 通过 MGF 积分表示计算（见底部 helper 函数 compute_rician_ber）
    BER_Ric_ana(idx) = compute_rician_ber(gb, K);
    
    %% 5) 随机部署用户 + 瑞利衰落 + AWGN
    n_rayrand = sqrt(N0/2) * (randn(N_samples, 1) + 1j*randn(N_samples, 1));
    r_rayrand = h_rayleigh_rand * sqrt(Eb) .* bits + n_rayrand;
    r_eq_rayrand = real(conj(h_rayleigh_rand) .* r_rayrand);
    BER_RayRand_sim(idx) = mean(r_eq_rayrand .* bits < 0);
    % 理论 BER：对距离 d 求条件期望
    % 条件 BER（给定距离 d）为 0.5*(1 - sqrt(γ_b/(d^ν + γ_b)))
    % 距离 d 的 PDF 为 2d/R^2，因此整体 BER = (2/R^2) * ∫ d * P_b(d) dd
    ber_cond_ray = @(d) 0.5 * (1 - sqrt(gb ./ (d.^nu + gb)));
    BER_RayRand_ana(idx) = (2/R^2) * integral(@(d) d .* ber_cond_ray(d), 0, R);
    
    %% 6) 随机部署用户 + 莱斯衰落 + AWGN
    n_ricrand = sqrt(N0/2) * (randn(N_samples, 1) + 1j*randn(N_samples, 1));
    r_ricrand = h_rician_rand * sqrt(Eb) .* bits + n_ricrand;
    r_eq_ricrand = real(conj(h_rician_rand) .* r_ricrand);
    BER_RicRand_sim(idx) = mean(r_eq_ricrand .* bits < 0);
    % 理论 BER：对每个距离 d，将有效 SNR 视为 γ_b/d^ν，调用 compute_rician_ber，再对 d 积分
    ber_cond_ric = @(d) arrayfun(@(dd) compute_rician_ber(gb/max(dd,eps)^nu, K), d);
    BER_RicRand_ana(idx) = (2/R^2) * integral(@(d) d .* ber_cond_ric(d), 0, R);
    
    fprintf('已完成 SNR = %d dB 的 BER 仿真\n', SNR_dB(idx));
end

%% ===================== 任务 3 & 4：Outage 概率仿真 =====================
% Outage 定义：瞬时信道容量 C_inst = log2(1 + SNR_inst) < C 的概率
% 等价于瞬时 SNR < γ_th = 2^C - 1 的概率

% 预分配 outage 结果数组
Pout_Ray_sim = zeros(size(SNR_dB));
Pout_Ray_ana = zeros(size(SNR_dB));
Pout_Ric_sim = zeros(size(SNR_dB));
Pout_Ric_ana = zeros(size(SNR_dB));
Pout_RayRand_sim = zeros(size(SNR_dB));
Pout_RayRand_ana = zeros(size(SNR_dB));
Pout_RicRand_sim = zeros(size(SNR_dB));
Pout_RicRand_ana = zeros(size(SNR_dB));

for idx = 1:length(SNR_dB)
    gb = SNR_lin(idx);
    T = gamma_th / gb;              % outage 门限折算到信道功率增益 |h|^2 上
                                    % 因为 SNR_inst = γ_b * |h|^2，所以 |h|^2 需超过 T = γ_th/γ_b
    
    %% 1) 瑞利衰落（d = 1 m）
    inst_snr_ray = gb * abs(h_rayleigh).^2;
    Pout_Ray_sim(idx) = mean(log2(1 + inst_snr_ray) < C);
    % 理论 outage：|h_s|^2 ~ Exp(1)，故 P_out = Pr(|h_s|^2 < T) = 1 - exp(-T)
    Pout_Ray_ana(idx) = 1 - exp(-T);
    
    %% 2) 莱斯衰落（d = 1 m）
    inst_snr_ric = gb * abs(h_rician).^2;
    Pout_Ric_sim(idx) = mean(log2(1 + inst_snr_ric) < C);
    % 理论 outage：|h|^2 服从非中心卡方分布，可用 Marcum Q 函数
    % P_out = 1 - Q_1(sqrt(2K), sqrt(2(K+1)T))
    % MATLAB 中 marcumq(a, b) 即 Q_1(a, b)
    Pout_Ric_ana(idx) = 1 - marcumq(sqrt(2*K), sqrt(2*(K+1)*T));
    
    %% 3) 随机部署用户 + 瑞利衰落
    inst_snr_rayrand = gb * abs(h_rayleigh).^2 ./ d_user.^nu;
    Pout_RayRand_sim(idx) = mean(log2(1 + inst_snr_rayrand) < C);
    % 理论 outage：利用圆盘均匀部署的 PDF 和指数衰落积分得到闭式
    % P_out = 1 - (2/(ν*R^2)) * T^(-2/ν) * γ(2/ν, T*R^ν)
    % 其中 γ(s, x) 为下不完全伽马函数，MATLAB 用 gammainc(x, s, 'lower')*gamma(s)
    Pout_RayRand_ana(idx) = 1 - (2/(nu*R^2)) * T^(-2/nu) * gamma(2/nu) * gammainc(T*R^nu, 2/nu, 'lower');
    
    %% 4) 随机部署用户 + 莱斯衰落
    inst_snr_ricrand = gb * abs(h_rician).^2 ./ d_user.^nu;
    Pout_RicRand_sim(idx) = mean(log2(1 + inst_snr_ricrand) < C);
    % 理论 outage：对每个距离 d 的条件 outage 用 Marcum Q 函数，再对 d 积分
    pout_cond_ric = @(d) 1 - marcumq(sqrt(2*K), sqrt(2*(K+1)*T*(d.^nu + eps)));
    Pout_RicRand_ana(idx) = (2/R^2) * integral(@(d) d .* pout_cond_ric(d), 0, R);
    
    fprintf('已完成 SNR = %d dB 的 Outage 仿真\n', SNR_dB(idx));
end

fprintf('\n');

%% ===================== 百分比偏差表格 =====================
% 计算仿真值与理论值之间的百分比偏差，用于量化蒙特卡洛估计的准确度。
% 为避免理论值为 0 时除零，分母取 max(ana, eps)，eps 为 MATLAB 最小正浮点数。
pct_err = @(sim, ana) abs(sim - ana) ./ max(ana, eps) * 100;

% 打印 BER 百分比偏差表
fprintf('--- BER 仿真值与理论值的百分比偏差 ---\n');
fprintf('SNR(dB) | AWGN    Laplace  Rayleigh Rician   RayRand  RicRand\n');
for idx = 1:length(SNR_dB)
    fprintf('%5d   | %6.3f%% %6.3f%% %6.3f%% %6.3f%% %6.3f%% %6.3f%%\n', ...
        SNR_dB(idx), ...
        pct_err(BER_AWGN_sim(idx), BER_AWGN_ana(idx)), ...
        pct_err(BER_Lap_sim(idx), BER_Lap_ana(idx)), ...
        pct_err(BER_Ray_sim(idx), BER_Ray_ana(idx)), ...
        pct_err(BER_Ric_sim(idx), BER_Ric_ana(idx)), ...
        pct_err(BER_RayRand_sim(idx), BER_RayRand_ana(idx)), ...
        pct_err(BER_RicRand_sim(idx), BER_RicRand_ana(idx)));
end

% 打印 Outage 百分比偏差表
fprintf('\n--- Outage 概率仿真值与理论值的百分比偏差 ---\n');
fprintf('SNR(dB) | Rayleigh Rician   RayRand  RicRand\n');
for idx = 1:length(SNR_dB)
    fprintf('%5d   | %6.3f%% %6.3f%% %6.3f%% %6.3f%%\n', ...
        SNR_dB(idx), ...
        pct_err(Pout_Ray_sim(idx), Pout_Ray_ana(idx)), ...
        pct_err(Pout_Ric_sim(idx), Pout_Ric_ana(idx)), ...
        pct_err(Pout_RayRand_sim(idx), Pout_RayRand_ana(idx)), ...
        pct_err(Pout_RicRand_sim(idx), Pout_RicRand_ana(idx)));
end
fprintf('\n');

%% ===================== 任务 5：图 1 - 仅噪声信道的 BER =====================
figure('Name', 'BER vs SNR - 仅噪声信道', 'Position', [100 550 800 600]);
% semilogy：Y 轴对数坐标，适合展示 BER 数量级变化
semilogy(SNR_dB, BER_AWGN_ana, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
% 仿真 BER 可能为 0，对数坐标会出错，因此用 max(BER, 1e-7) 截断下限
semilogy(SNR_dB, max(BER_AWGN_sim, 1e-7), 'b--s', 'LineWidth', 1.5, 'MarkerSize', 6);
semilogy(SNR_dB, BER_Lap_ana, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 6);
semilogy(SNR_dB, max(BER_Lap_sim, 1e-7), 'r--s', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('平均 SNR \gamma_b (dB)');
ylabel('误码率 (BER)');
title('AWGN 与拉普拉斯噪声信道下的 BER 对比');
legend('AWGN 理论', 'AWGN 仿真', ...
       'Laplacian 理论', 'Laplacian 仿真', ...
       'Location', 'southwest');
axis([0 15 1e-5 1]);

%% ===================== 任务 5：图 2 - 衰落信道的 BER =====================
figure('Name', 'BER vs SNR - 衰落信道', 'Position', [150 500 800 600]);
semilogy(SNR_dB, BER_Ray_ana, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5);
hold on;
semilogy(SNR_dB, max(BER_Ray_sim, 1e-7), 'b--s', 'LineWidth', 1.5, 'MarkerSize', 5);
semilogy(SNR_dB, BER_Ric_ana, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 5);
semilogy(SNR_dB, max(BER_Ric_sim, 1e-7), 'r--s', 'LineWidth', 1.5, 'MarkerSize', 5);
semilogy(SNR_dB, BER_RayRand_ana, 'g-^', 'LineWidth', 1.5, 'MarkerSize', 5);
semilogy(SNR_dB, max(BER_RayRand_sim, 1e-7), 'g--v', 'LineWidth', 1.5, 'MarkerSize', 5);
semilogy(SNR_dB, BER_RicRand_ana, 'm-d', 'LineWidth', 1.5, 'MarkerSize', 5);
semilogy(SNR_dB, max(BER_RicRand_sim, 1e-7), 'm--p', 'LineWidth', 1.5, 'MarkerSize', 5);
grid on;
xlabel('平均 SNR \gamma_b (dB)');
ylabel('误码率 (BER)');
title('衰落信道下 BPSK 的 BER 对比');
legend('Rayleigh 理论', 'Rayleigh 仿真', ...
       'Rician 理论', 'Rician 仿真', ...
       '随机 Rayleigh 理论', '随机 Rayleigh 仿真', ...
       '随机 Rician 理论', '随机 Rician 仿真', ...
       'Location', 'southwest');
axis([0 15 1e-5 1]);

%% ===================== 任务 5：图 3 - Outage 概率 =====================
figure('Name', 'Outage 概率 vs SNR', 'Position', [200 450 800 600]);
semilogy(SNR_dB, Pout_Ray_ana, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5);
hold on;
semilogy(SNR_dB, max(Pout_Ray_sim, 1e-7), 'b--s', 'LineWidth', 1.5, 'MarkerSize', 5);
semilogy(SNR_dB, Pout_Ric_ana, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 5);
semilogy(SNR_dB, max(Pout_Ric_sim, 1e-7), 'r--s', 'LineWidth', 1.5, 'MarkerSize', 5);
semilogy(SNR_dB, Pout_RayRand_ana, 'g-^', 'LineWidth', 1.5, 'MarkerSize', 5);
semilogy(SNR_dB, max(Pout_RayRand_sim, 1e-7), 'g--v', 'LineWidth', 1.5, 'MarkerSize', 5);
semilogy(SNR_dB, Pout_RicRand_ana, 'm-d', 'LineWidth', 1.5, 'MarkerSize', 5);
semilogy(SNR_dB, max(Pout_RicRand_sim, 1e-7), 'm--p', 'LineWidth', 1.5, 'MarkerSize', 5);
grid on;
xlabel('平均 SNR \gamma_b (dB)');
ylabel('Outage 概率 P_{out}');
title('Outage 概率 vs SNR (C = 1.2 bps/Hz)');
legend('Rayleigh 理论', 'Rayleigh 仿真', ...
       'Rician 理论', 'Rician 仿真', ...
       '随机 Rayleigh 理论', '随机 Rayleigh 仿真', ...
       '随机 Rician 理论', '随机 Rician 仿真', ...
       'Location', 'southwest');
axis([0 15 1e-5 1]);

%% ===================== 可选：样本数增加带来的误差减小演示 =====================
% 本部分通过逐步增加样本数，直观展示蒙特卡洛估计的收敛性。
% 以 AWGN 信道在 10 dB 处的 BER 为例，分别用 N = 1e5, 5e5, 1e6, 5e6 进行估计。
fprintf('--- 演示：样本数越多，BER 偏差越小 ---\n');
test_SNR_dB = 10;
test_gb = 10^(test_SNR_dB/10);
test_Eb = test_gb * N0;
test_bits = 2*(randi([0 1], 5e6, 1)) - 1;
test_n = sqrt(N0/2) * (randn(5e6, 1) + 1j*randn(5e6, 1));
test_r = sqrt(test_Eb) * test_bits + test_n;
test_ana = 0.5 * erfc(sqrt(test_gb));
Ns_test = [1e5 5e5 1e6 5e6];
for Nt = Ns_test
    ber_test = mean(real(test_r(1:Nt)) .* test_bits(1:Nt) < 0);
    err_test = abs(ber_test - test_ana) / test_ana * 100;
    fprintf('N = %7.0f : 仿真 BER = %.6e, 偏差 = %.4f%%\n', Nt, ber_test, err_test);
end

%% ===================== 辅助函数 =====================
function ber = compute_rician_ber(gb, K)
    % compute_rician_ber: 计算 BPSK 在莱斯衰落信道下的理论 BER
    % 输入:
    %   gb - 平均 SNR（线性值）
    %   K  - 莱斯 K 因子
    % 输出:
    %   ber - 理论误码率
    %
    % 采用 MGF（矩生成函数）积分表示：
    % P_b = (1/π) ∫_0^(π/2) M_γ(-1/sin²θ) dθ
    % 其中莱斯衰落 SNR 的 MGF 为：
    % M_γ(s) = (K+1)/(K+1 - s*gb) * exp(-K*gb*s/(K+1 - s*gb))
    % 代入 s = -1/sin²θ 后得到被积函数。
    % 对于极高 SNR，BER 数值上接近 0，直接返回 0 以避免数值积分不稳定。
    if gb > 1e4
        ber = 0;
        return;
    end
    fun = @(theta) (1/pi) * ((K+1) .* sin(theta).^2) ./ ((K+1).*sin(theta).^2 + gb) ...
                   .* exp(-K*gb ./ ((K+1).*sin(theta).^2 + gb));
    ber = integral(fun, 0, pi/2);
end
