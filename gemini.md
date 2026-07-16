# %% ELEC9123: 设计任务 B（无线通信）
% 通过蒙特卡洛仿真验证 BPSK 系统在不同无线信道模型下的性能指标
%
% 作者: YU
% 学号: z5615351
% 文件名: z5615351\_YU\_DTB\_2026\.m
%
% 本脚本实现了六种无线信道模型下 BPSK 系统的蒙特卡洛仿真，
% 并将误码率（BER）和 outage 概率与理论解析式进行对比验证。
% 六种信道模型包括：
%   1\) 仅 AWGN 信道
%   2\) 仅拉普拉斯（Laplacian）脉冲噪声信道
%   3\) 瑞利（Rayleigh）衰落 \+ AWGN
%   4\) 莱斯（Rician）衰落 \+ AWGN
%   5\) 随机部署用户 \+ 瑞利衰落 \+ AWGN
%   6\) 随机部署用户 \+ 莱斯衰落 \+ AWGN

%% ===================== 环境清理 =====================
clc;        % 清空命令行窗口，保持输出整洁
clear;      % 清除工作区所有变量，避免历史变量干扰
close all;  % 关闭所有已打开的 Figure 窗口

%% ===================== 系统参数设置 =====================
% 下面所有参数均来自设计任务说明书（Task B）中的默认设定

C = 1\.2;                         % 中断门限容量，单位 bps/Hz
                                 % 当瞬时信道容量 log2\(1\+SNR\) < C 时，认为发生 outage
SNR\_dB = 0:1:15;                 % 平均信噪比（每比特）的扫描范围，单位 dB，从 0 dB 到 15 dB，步长 1 dB
SNR\_lin = 10\.^\(SNR\_dB/10\);       % 将 dB 信噪比转换为线性值 γ\_b（平均 SNR per bit）
                                 % 公式：γ\_b\(linear\) = 10^\(γ\_b\(dB\)/10\)
R = 3;                           % 小区（圆形覆盖区域）半径，单位 m
nu = 2\.2;                        % 路径损耗指数（path\-loss exponent），典型取值 2~4
hd = 1;                          % 直射径（LOS）分量的幅度，取单位幅度
K = 5;                           % 莱斯 K 因子：直射径功率与散射径功率之比
N0 = 1;                          % 噪声功率/方差，这里取归一化值 1
N\_samples = 1e6;                 % 每个 SNR 点的蒙特卡洛样本数（题目要求至少约 10^6）
gamma\_th = 2^C \- 1;              % outage 对应的 SNR 门限
                                 % 由 C = log2\(1\+γ\_th\) 可得 γ\_th = 2^C \- 1

% 打印任务信息
fprintf\('ELEC9123 设计任务 B \- 蒙特卡洛仿真\\n'\);
fprintf\('作者: YU \(z5615351\)\\n'\);
fprintf\('每个 SNR 点的样本数: %d\\n\\n', N\_samples\);

%% ===================== 任务 2：验证拉普拉斯随机变量 =====================
% 目标：使用逆 CDF 法生成服从拉普拉斯分布的噪声，并与理论 PDF/CDF 对比。
%
% 拉普拉斯分布 Laplace\(0, b\) 的方差为 2\*b^2。
% 题目给定噪声方差为 N0，因此 2\*b^2 = N0，得到 b = sqrt\(N0/2\)。
% 拉普拉斯分布的 CDF 为：
%   F\(x\) = 0\.5 \* exp\(sqrt\(2/N0\)\*x\),            x <= 0
%   F\(x\) = 1 \- 0\.5 \* exp\(\-sqrt\(2/N0\)\*x\),       x > 0
% 其逆 CDF（用于由均匀随机变量 u∈\[0,1\] 生成 x）为：
%   x = \-b \* sign\(u\-0\.5\) \* log\(1 \- 2\*|u\-0\.5|\)

N\_lap = 1e6;                     % 用于验证拉普拉斯分布的样本数
b\_lap = sqrt\(N0/2\);              % 拉普拉斯分布的尺度参数 b
u\_lap = rand\(N\_lap, 1\);          % 生成 \[0,1\] 均匀分布随机变量
% 逆 CDF 变换：生成服从 Laplace\(0, b\_lap\) 的随机变量
n\_lap = \-b\_lap \* sign\(u\_lap \- 0\.5\) \.\* log\(1 \- 2\*abs\(u\_lap \- 0\.5\)\);

% 解析 PDF 和 CDF 的绘制点
x\_lap = linspace\(\-6\*b\_lap, 6\*b\_lap, 1000\);
% 理论 PDF：f\(x\) = \(1/\(2b\)\) \* exp\(\-|x|/b\)
f\_lap\_ana = \(1/\(2\*b\_lap\)\) \* exp\(\-abs\(x\_lap\)/b\_lap\);
% 理论 CDF：分段定义
F\_lap\_ana = zeros\(size\(x\_lap\)\);
F\_lap\_ana\(x\_lap <= 0\) = 0\.5 \* exp\(x\_lap\(x\_lap <= 0\)/b\_lap\);
F\_lap\_ana\(x\_lap > 0\) = 1 \- 0\.5 \* exp\(\-x\_lap\(x\_lap > 0\)/b\_lap\);

% 仿真经验 CDF：对每个 x 值，统计样本中小于等于 x 的比例
F\_lap\_sim = zeros\(size\(x\_lap\)\);
for i = 1:length\(x\_lap\)
    F\_lap\_sim\(i\) = mean\(n\_lap <= x\_lap\(i\)\);
end

% 绘制拉普拉斯分布验证图（PDF 和 CDF 各一个子图）
figure\('Name', '拉普拉斯分布验证', 'Position', \[100 100 900 400\]\);

subplot\(1, 2, 1\);
histogram\(n\_lap, 100, 'Normalization', 'pdf', 'FaceColor', \[0\.3 0\.6 0\.9\]\);
hold on;
plot\(x\_lap, f\_lap\_ana, 'r\-', 'LineWidth', 2\);
xlabel\('x'\);
ylabel\('PDF f\_n\(x\)'\);
title\('拉普拉斯噪声概率密度函数 \(PDF\)'\);
legend\('仿真直方图', '理论 PDF', 'Location', 'best'\);
grid on;

subplot\(1, 2, 2\);
plot\(x\_lap, F\_lap\_sim, 'b\-', 'LineWidth', 1\.5\);
hold on;
plot\(x\_lap, F\_lap\_ana, 'r\-\-', 'LineWidth', 2\);
xlabel\('x'\);
ylabel\('CDF F\_n\(x\)'\);
title\('拉普拉斯噪声累积分布函数 \(CDF\)'\);
legend\('仿真 CDF', '理论 CDF', 'Location', 'best'\);
grid on;

% 计算仿真 CDF 与理论 CDF 之间的最大绝对误差（类似 Kolmogorov\-Smirnov 统计量）
max\_cdf\_error = max\(abs\(F\_lap\_sim \- F\_lap\_ana\)\);
fprintf\('任务 2：拉普拉斯 CDF 最大绝对误差 = %\.4e\\n\\n', max\_cdf\_error\);

%% ===================== 预先生成公共随机变量 =====================
% 下面这些随机变量与 SNR 无关，因此在进入 SNR 循环之前一次性生成，
% 可显著提高效率并保证不同 SNR 下使用同一组信道实现，便于对比。

% BPSK 调制符号：等概率取 \+1 或 \-1
bits = 2\*\(randi\(\[0 1\], N\_samples, 1\)\) \- 1;

% 瑞利衰落系数 h\_s：复高斯 CN\(0,1\)，这里除以 sqrt\(2\) 使 E\{|h\_s|^2\}=1
% |h\_s| 服从瑞利分布，|h\_s|^2 服从指数分布 Exp\(1\)
h\_rayleigh = \(randn\(N\_samples, 1\) \+ 1j\*randn\(N\_samples, 1\)\) / sqrt\(2\);

% 莱斯衰落系数 h：直射径 \+ 散射径
% h = sqrt\(K/\(K\+1\)\)\*hd \+ sqrt\(1/\(K\+1\)\)\*h\_s
% 其中 hd 为确定性 LOS 分量，h\_s 为 CN\(0,1\) 散射分量。
% 这样构造满足 E\{|h|^2\} = K/\(K\+1\)\*|hd|^2 \+ 1/\(K\+1\)\*1 = 1
h\_los = sqrt\(K/\(K\+1\)\) \* hd;
h\_scat = sqrt\(1/\(K\+1\)\) \* \(randn\(N\_samples, 1\) \+ 1j\*randn\(N\_samples, 1\)\) / sqrt\(2\);
h\_rician = h\_los \+ h\_scat;

% 随机用户部署距离 d：在半径为 R 的圆盘内均匀分布
% 圆盘内均匀分布时，距离 d 的 PDF 为 f\_d\(d\) = 2d/R^2, 0 <= d <= R
% 生成方法：d = R \* sqrt\(rand\)，因为 CDF F\_d\(d\) = d^2/R^2，逆 CDF 即 d = R\*sqrt\(u\)
d\_user = R \* sqrt\(rand\(N\_samples, 1\)\);

% 将衰落与随机距离合并，得到包含路径损耗 d^\{\-nu/2\} 的等效信道
% 等效信道增益 = h / d^\(nu/2\)，因此接收功率额外乘以 d^\{\-nu\}
h\_rayleigh\_rand = h\_rayleigh \./ \(d\_user\.^\(nu/2\)\);
h\_rician\_rand = h\_rician \./ \(d\_user\.^\(nu/2\)\);

%% ===================== 任务 3 & 4：BER 仿真 =====================
% 对六种信道模型分别进行蒙特卡洛 BER 仿真，并与理论解析值对比。

% 预分配存储 BER 结果的数组
BER\_AWGN\_sim = zeros\(size\(SNR\_dB\)\);   % AWGN 信道仿真 BER
BER\_AWGN\_ana = zeros\(size\(SNR\_dB\)\);   % AWGN 信道理论 BER
BER\_Lap\_sim  = zeros\(size\(SNR\_dB\)\);   % 拉普拉斯噪声信道仿真 BER
BER\_Lap\_ana  = zeros\(size\(SNR\_dB\)\);   % 拉普拉斯噪声信道理论 BER
BER\_Ray\_sim  = zeros\(size\(SNR\_dB\)\);   % 瑞利衰落\+AWGN 仿真 BER
BER\_Ray\_ana  = zeros\(size\(SNR\_dB\)\);   % 瑞利衰落\+AWGN 理论 BER
BER\_Ric\_sim  = zeros\(size\(SNR\_dB\)\);   % 莱斯衰落\+AWGN 仿真 BER
BER\_Ric\_ana  = zeros\(size\(SNR\_dB\)\);   % 莱斯衰落\+AWGN 理论 BER
BER\_RayRand\_sim = zeros\(size\(SNR\_dB\)\);% 随机用户\+瑞利衰落 仿真 BER
BER\_RayRand\_ana = zeros\(size\(SNR\_dB\)\);% 随机用户\+瑞利衰落 理论 BER
BER\_RicRand\_sim = zeros\(size\(SNR\_dB\)\);% 随机用户\+莱斯衰落 仿真 BER
BER\_RicRand\_ana = zeros\(size\(SNR\_dB\)\);% 随机用户\+莱斯衰落 理论 BER

% 遍历每个 SNR 点
for idx = 1:length\(SNR\_dB\)
    gb = SNR\_lin\(idx\);              % 当前平均 SNR（线性值）
    Eb = gb \* N0;                   % 每比特能量：E\_b = γ\_b \* N\_0
    
    %% 1\) AWGN 信道
    % 复高斯噪声：实部与虚部各为 N\(0, N0/2\)，总方差为 N0
    n\_awgn = sqrt\(N0/2\) \* \(randn\(N\_samples, 1\) \+ 1j\*randn\(N\_samples, 1\)\);
    % 接收信号：r = sqrt\(Eb\)\*s \+ n
    r\_awgn = sqrt\(Eb\) \* bits \+ n\_awgn;
    % BPSK 判决：real\(r\)\*s < 0 表示判决错误（匹配滤波后取实部）
    BER\_AWGN\_sim\(idx\) = mean\(real\(r\_awgn\) \.\* bits < 0\);
    % 理论 BER：P\_b = 0\.5 \* erfc\(sqrt\(γ\_b\)\)
    BER\_AWGN\_ana\(idx\) = 0\.5 \* erfc\(sqrt\(gb\)\);
    
    %% 2\) 拉普拉斯脉冲噪声信道
    % 同样使用逆 CDF 法生成拉普拉斯噪声
    u = rand\(N\_samples, 1\);
    n\_lap = \-b\_lap \* sign\(u \- 0\.5\) \.\* log\(1 \- 2\*abs\(u \- 0\.5\)\);
    r\_lap = sqrt\(Eb\) \* bits \+ n\_lap;
    BER\_Lap\_sim\(idx\) = mean\(real\(r\_lap\) \.\* bits < 0\);
    % 理论 BER：P\_b = 0\.5 \* exp\(\-sqrt\(2\*γ\_b\)\)
    BER\_Lap\_ana\(idx\) = 0\.5 \* exp\(\-sqrt\(2\*gb\)\);
    
    %% 3\) 瑞利衰落 \+ AWGN
    n\_ray = sqrt\(N0/2\) \* \(randn\(N\_samples, 1\) \+ 1j\*randn\(N\_samples, 1\)\);
    r\_ray = h\_rayleigh \* sqrt\(Eb\) \.\* bits \+ n\_ray;
    % 等效接收（最大比合并/匹配滤波）：将接收信号乘以信道共轭并取实部
    r\_eq\_ray = real\(conj\(h\_rayleigh\) \.\* r\_ray\);
    BER\_Ray\_sim\(idx\) = mean\(r\_eq\_ray \.\* bits < 0\);
    % 理论 BER：P\_b = 0\.5 \* \(1 \- sqrt\(γ\_b/\(1\+γ\_b\)\)\)
    BER\_Ray\_ana\(idx\) = 0\.5 \* \(1 \- sqrt\(gb \./ \(1 \+ gb\)\)\);
    
    %% 4\) 莱斯衰落 \+ AWGN
    n\_ric = sqrt\(N0/2\) \* \(randn\(N\_samples, 1\) \+ 1j\*randn\(N\_samples, 1\)\);
    r\_ric = h\_rician \* sqrt\(Eb\) \.\* bits \+ n\_ric;
    r\_eq\_ric = real\(conj\(h\_rician\) \.\* r\_ric\);
    BER\_Ric\_sim\(idx\) = mean\(r\_eq\_ric \.\* bits < 0\);
    % 理论 BER 通过 MGF 积分表示计算（见底部 helper 函数 compute\_rician\_ber）
    BER\_Ric\_ana\(idx\) = compute\_rician\_ber\(gb, K\);
    
    %% 5\) 随机部署用户 \+ 瑞利衰落 \+ AWGN
    n\_rayrand = sqrt\(N0/2\) \* \(randn\(N\_samples, 1\) \+ 1j\*randn\(N\_samples, 1\)\);
    r\_rayrand = h\_rayleigh\_rand \* sqrt\(Eb\) \.\* bits \+ n\_rayrand;
    r\_eq\_rayrand = real\(conj\(h\_rayleigh\_rand\) \.\* r\_rayrand\);
    BER\_RayRand\_sim\(idx\) = mean\(r\_eq\_rayrand \.\* bits < 0\);
    % 理论 BER：对距离 d 求条件期望
    % 条件 BER（给定距离 d）为 0\.5\*\(1 \- sqrt\(γ\_b/\(d^ν \+ γ\_b\)\)\)
    % 距离 d 的 PDF 为 2d/R^2，因此整体 BER = \(2/R^2\) \* ∫ d \* P\_b\(d\) dd
    ber\_cond\_ray = @\(d\) 0\.5 \* \(1 \- sqrt\(gb \./ \(d\.^nu \+ gb\)\)\);
    BER\_RayRand\_ana\(idx\) = \(2/R^2\) \* integral\(@\(d\) d \.\* ber\_cond\_ray\(d\), 0, R\);
    
    %% 6\) 随机部署用户 \+ 莱斯衰落 \+ AWGN
    n\_ricrand = sqrt\(N0/2\) \* \(randn\(N\_samples, 1\) \+ 1j\*randn\(N\_samples, 1\)\);
    r\_ricrand = h\_rician\_rand \* sqrt\(Eb\) \.\* bits \+ n\_ricrand;
    r\_eq\_ricrand = real\(conj\(h\_rician\_rand\) \.\* r\_ricrand\);
    BER\_RicRand\_sim\(idx\) = mean\(r\_eq\_ricrand \.\* bits < 0\);
    % 理论 BER：对每个距离 d，将有效 SNR 视为 γ\_b/d^ν，调用 compute\_rician\_ber，再对 d 积分
    ber\_cond\_ric = @\(d\) arrayfun\(@\(dd\) compute\_rician\_ber\(gb/max\(dd,eps\)^nu, K\), d\);
    BER\_RicRand\_ana\(idx\) = \(2/R^2\) \* integral\(@\(d\) d \.\* ber\_cond\_ric\(d\), 0, R\);
    
    fprintf\('已完成 SNR = %d dB 的 BER 仿真\\n', SNR\_dB\(idx\)\);
end

%% ===================== 任务 3 & 4：Outage 概率仿真 =====================
% Outage 定义：瞬时信道容量 C\_inst = log2\(1 \+ SNR\_inst\) < C 的概率
% 等价于瞬时 SNR < γ\_th = 2^C \- 1 的概率

% 预分配 outage 结果数组
Pout\_Ray\_sim = zeros\(size\(SNR\_dB\)\);
Pout\_Ray\_ana = zeros\(size\(SNR\_dB\)\);
Pout\_Ric\_sim = zeros\(size\(SNR\_dB\)\);
Pout\_Ric\_ana = zeros\(size\(SNR\_dB\)\);
Pout\_RayRand\_sim = zeros\(size\(SNR\_dB\)\);
Pout\_RayRand\_ana = zeros\(size\(SNR\_dB\)\);
Pout\_RicRand\_sim = zeros\(size\(SNR\_dB\)\);
Pout\_RicRand\_ana = zeros\(size\(SNR\_dB\)\);

for idx = 1:length\(SNR\_dB\)
    gb = SNR\_lin\(idx\);
    T = gamma\_th / gb;              % outage 门限折算到信道功率增益 |h|^2 上
                                    % 因为 SNR\_inst = γ\_b \* |h|^2，所以 |h|^2 需超过 T = γ\_th/γ\_b
    
    %% 1\) 瑞利衰落（d = 1 m）
    inst\_snr\_ray = gb \* abs\(h\_rayleigh\)\.^2;
    Pout\_Ray\_sim\(idx\) = mean\(log2\(1 \+ inst\_snr\_ray\) < C\);
    % 理论 outage：|h\_s|^2 ~ Exp\(1\)，故 P\_out = Pr\(|h\_s|^2 < T\) = 1 \- exp\(\-T\)
    Pout\_Ray\_ana\(idx\) = 1 \- exp\(\-T\);
    
    %% 2\) 莱斯衰落（d = 1 m）
    inst\_snr\_ric = gb \* abs\(h\_rician\)\.^2;
    Pout\_Ric\_sim\(idx\) = mean\(log2\(1 \+ inst\_snr\_ric\) < C\);
    % 理论 outage：|h|^2 服从非中心卡方分布，可用 Marcum Q 函数
    % P\_out = 1 \- Q\_1\(sqrt\(2K\), sqrt\(2\(K\+1\)T\)\)
    % MATLAB 中 marcumq\(a, b\) 即 Q\_1\(a, b\)
    Pout\_Ric\_ana\(idx\) = 1 \- marcumq\(sqrt\(2\*K\), sqrt\(2\*\(K\+1\)\*T\)\);
    
    %% 3\) 随机部署用户 \+ 瑞利衰落
    inst\_snr\_rayrand = gb \* abs\(h\_rayleigh\)\.^2 \./ d\_user\.^nu;
    Pout\_RayRand\_sim\(idx\) = mean\(log2\(1 \+ inst\_snr\_rayrand\) < C\);
    % 理论 outage：利用圆盘均匀部署的 PDF 和指数衰落积分得到闭式
    % P\_out = 1 \- \(2/\(ν\*R^2\)\) \* T^\(\-2/ν\) \* γ\(2/ν, T\*R^ν\)
    % 其中 γ\(s, x\) 为下不完全伽马函数，MATLAB 用 gammainc\(x, s, 'lower'\)\*gamma\(s\)
    Pout\_RayRand\_ana\(idx\) = 1 \- \(2/\(nu\*R^2\)\) \* T^\(\-2/nu\) \* gamma\(2/nu\) \* gammainc\(T\*R^nu, 2/nu, 'lower'\);
    
    %% 4\) 随机部署用户 \+ 莱斯衰落
    inst\_snr\_ricrand = gb \* abs\(h\_rician\)\.^2 \./ d\_user\.^nu;
    Pout\_RicRand\_sim\(idx\) = mean\(log2\(1 \+ inst\_snr\_ricrand\) < C\);
    % 理论 outage：对每个距离 d 的条件 outage 用 Marcum Q 函数，再对 d 积分
    pout\_cond\_ric = @\(d\) 1 \- marcumq\(sqrt\(2\*K\), sqrt\(2\*\(K\+1\)\*T\*\(d\.^nu \+ eps\)\)\);
    Pout\_RicRand\_ana\(idx\) = \(2/R^2\) \* integral\(@\(d\) d \.\* pout\_cond\_ric\(d\), 0, R\);
    
    fprintf\('已完成 SNR = %d dB 的 Outage 仿真\\n', SNR\_dB\(idx\)\);
end

fprintf\('\\n'\);

%% ===================== 百分比偏差表格 =====================
% 计算仿真值与理论值之间的百分比偏差，用于量化蒙特卡洛估计的准确度。
% 为避免理论值为 0 时除零，分母取 max\(ana, eps\)，eps 为 MATLAB 最小正浮点数。
pct\_err = @\(sim, ana\) abs\(sim \- ana\) \./ max\(ana, eps\) \* 100;

% 打印 BER 百分比偏差表
fprintf\('\-\-\- BER 仿真值与理论值的百分比偏差 \-\-\-\\n'\);
fprintf\('SNR\(dB\) | AWGN    Laplace  Rayleigh Rician   RayRand  RicRand\\n'\);
for idx = 1:length\(SNR\_dB\)
    fprintf\('%5d   | %6\.3f%% %6\.3f%% %6\.3f%% %6\.3f%% %6\.3f%% %6\.3f%%\\n', \.\.\.
        SNR\_dB\(idx\), \.\.\.
        pct\_err\(BER\_AWGN\_sim\(idx\), BER\_AWGN\_ana\(idx\)\), \.\.\.
        pct\_err\(BER\_Lap\_sim\(idx\), BER\_Lap\_ana\(idx\)\), \.\.\.
        pct\_err\(BER\_Ray\_sim\(idx\), BER\_Ray\_ana\(idx\)\), \.\.\.
        pct\_err\(BER\_Ric\_sim\(idx\), BER\_Ric\_ana\(idx\)\), \.\.\.
        pct\_err\(BER\_RayRand\_sim\(idx\), BER\_RayRand\_ana\(idx\)\), \.\.\.
        pct\_err\(BER\_RicRand\_sim\(idx\), BER\_RicRand\_ana\(idx\)\)\);
end

% 打印 Outage 百分比偏差表
fprintf\('\\n\-\-\- Outage 概率仿真值与理论值的百分比偏差 \-\-\-\\n'\);
fprintf\('SNR\(dB\) | Rayleigh Rician   RayRand  RicRand\\n'\);
for idx = 1:length\(SNR\_dB\)
    fprintf\('%5d   | %6\.3f%% %6\.3f%% %6\.3f%% %6\.3f%%\\n', \.\.\.
        SNR\_dB\(idx\), \.\.\.
        pct\_err\(Pout\_Ray\_sim\(idx\), Pout\_Ray\_ana\(idx\)\), \.\.\.
        pct\_err\(Pout\_Ric\_sim\(idx\), Pout\_Ric\_ana\(idx\)\), \.\.\.
        pct\_err\(Pout\_RayRand\_sim\(idx\), Pout\_RayRand\_ana\(idx\)\), \.\.\.
        pct\_err\(Pout\_RicRand\_sim\(idx\), Pout\_RicRand\_ana\(idx\)\)\);
end
fprintf\('\\n'\);

%% ===================== 任务 5：图 1 \- 仅噪声信道的 BER =====================
figure\('Name', 'BER vs SNR \- 仅噪声信道', 'Position', \[100 550 800 600\]\);
% semilogy：Y 轴对数坐标，适合展示 BER 数量级变化
semilogy\(SNR\_dB, BER\_AWGN\_ana, 'b\-o', 'LineWidth', 1\.5, 'MarkerSize', 6\);
hold on;
% 仿真 BER 可能为 0，对数坐标会出错，因此用 max\(BER, 1e\-7\) 截断下限
semilogy\(SNR\_dB, max\(BER\_AWGN\_sim, 1e\-7\), 'b\-\-s', 'LineWidth', 1\.5, 'MarkerSize', 6\);
semilogy\(SNR\_dB, BER\_Lap\_ana, 'r\-o', 'LineWidth', 1\.5, 'MarkerSize', 6\);
semilogy\(SNR\_dB, max\(BER\_Lap\_sim, 1e\-7\), 'r\-\-s', 'LineWidth', 1\.5, 'MarkerSize', 6\);
grid on;
xlabel\('平均 SNR \\gamma\_b \(dB\)'\);
ylabel\('误码率 \(BER\)'\);
title\('AWGN 与拉普拉斯噪声信道下的 BER 对比'\);
legend\('AWGN 理论', 'AWGN 仿真', \.\.\.
       'Laplacian 理论', 'Laplacian 仿真', \.\.\.
       'Location', 'southwest'\);
axis\(\[0 15 1e\-5 1\]\);

%% ===================== 任务 5：图 2 \- 衰落信道的 BER =====================
figure\('Name', 'BER vs SNR \- 衰落信道', 'Position', \[150 500 800 600\]\);
semilogy\(SNR\_dB, BER\_Ray\_ana, 'b\-o', 'LineWidth', 1\.5, 'MarkerSize', 5\);
hold on;
semilogy\(SNR\_dB, max\(BER\_Ray\_sim, 1e\-7\), 'b\-\-s', 'LineWidth', 1\.5, 'MarkerSize', 5\);
semilogy\(SNR\_dB, BER\_Ric\_ana, 'r\-o', 'LineWidth', 1\.5, 'MarkerSize', 5\);
semilogy\(SNR\_dB, max\(BER\_Ric\_sim, 1e\-7\), 'r\-\-s', 'LineWidth', 1\.5, 'MarkerSize', 5\);
semilogy\(SNR\_dB, BER\_RayRand\_ana, 'g\-^', 'LineWidth', 1\.5, 'MarkerSize', 5\);
semilogy\(SNR\_dB, max\(BER\_RayRand\_sim, 1e\-7\), 'g\-\-v', 'LineWidth', 1\.5, 'MarkerSize', 5\);
semilogy\(SNR\_dB, BER\_RicRand\_ana, 'm\-d', 'LineWidth', 1\.5, 'MarkerSize', 5\);
semilogy\(SNR\_dB, max\(BER\_RicRand\_sim, 1e\-7\), 'm\-\-p', 'LineWidth', 1\.5, 'MarkerSize', 5\);
grid on;
xlabel\('平均 SNR \\gamma\_b \(dB\)'\);
ylabel\('误码率 \(BER\)'\);
title\('衰落信道下 BPSK 的 BER 对比'\);
legend\('Rayleigh 理论', 'Rayleigh 仿真', \.\.\.
       'Rician 理论', 'Rician 仿真', \.\.\.
       '随机 Rayleigh 理论', '随机 Rayleigh 仿真', \.\.\.
       '随机 Rician 理论', '随机 Rician 仿真', \.\.\.
       'Location', 'southwest'\);
axis\(\[0 15 1e\-5 1\]\);

%% ===================== 任务 5：图 3 \- Outage 概率 =====================
figure\('Name', 'Outage 概率 vs SNR', 'Position', \[200 450 800 600\]\);
semilogy\(SNR\_dB, Pout\_Ray\_ana, 'b\-o', 'LineWidth', 1\.5, 'MarkerSize', 5\);
hold on;
semilogy\(SNR\_dB, max\(Pout\_Ray\_sim, 1e\-7\), 'b\-\-s', 'LineWidth', 1\.5, 'MarkerSize', 5\);
semilogy\(SNR\_dB, Pout\_Ric\_ana, 'r\-o', 'LineWidth', 1\.5, 'MarkerSize', 5\);
semilogy\(SNR\_dB, max\(Pout\_Ric\_sim, 1e\-7\), 'r\-\-s', 'LineWidth', 1\.5, 'MarkerSize', 5\);
semilogy\(SNR\_dB, Pout\_RayRand\_ana, 'g\-^', 'LineWidth', 1\.5, 'MarkerSize', 5\);
semilogy\(SNR\_dB, max\(Pout\_RayRand\_sim, 1e\-7\), 'g\-\-v', 'LineWidth', 1\.5, 'MarkerSize', 5\);
semilogy\(SNR\_dB, Pout\_RicRand\_ana, 'm\-d', 'LineWidth', 1\.5, 'MarkerSize', 5\);
semilogy\(SNR\_dB, max\(Pout\_RicRand\_sim, 1e\-7\), 'm\-\-p', 'LineWidth', 1\.5, 'MarkerSize', 5\);
grid on;
xlabel\('平均 SNR \\gamma\_b \(dB\)'\);
ylabel\('Outage 概率 P\_\{out\}'\);
title\('Outage 概率 vs SNR \(C = 1\.2 bps/Hz\)'\);
legend\('Rayleigh 理论', 'Rayleigh 仿真', \.\.\.
       'Rician 理论', 'Rician 仿真', \.\.\.
       '随机 Rayleigh 理论', '随机 Rayleigh 仿真', \.\.\.
       '随机 Rician 理论', '随机 Rician 仿真', \.\.\.
       'Location', 'southwest'\);
axis\(\[0 15 1e\-5 1\]\);

%% ===================== 可选：样本数增加带来的误差减小演示 =====================
% 本部分通过逐步增加样本数，直观展示蒙特卡洛估计的收敛性。
% 以 AWGN 信道在 10 dB 处的 BER 为例，分别用 N = 1e5, 5e5, 1e6, 5e6 进行估计。
fprintf\('\-\-\- 演示：样本数越多，BER 偏差越小 \-\-\-\\n'\);
test\_SNR\_dB = 10;
test\_gb = 10^\(test\_SNR\_dB/10\);
test\_Eb = test\_gb \* N0;
test\_bits = 2\*\(randi\(\[0 1\], 5e6, 1\)\) \- 1;
test\_n = sqrt\(N0/2\) \* \(randn\(5e6, 1\) \+ 1j\*randn\(5e6, 1\)\);
test\_r = sqrt\(test\_Eb\) \* test\_bits \+ test\_n;
test\_ana = 0\.5 \* erfc\(sqrt\(test\_gb\)\);
Ns\_test = \[1e5 5e5 1e6 5e6\];
for Nt = Ns\_test
    ber\_test = mean\(real\(test\_r\(1:Nt\)\) \.\* test\_bits\(1:Nt\) < 0\);
    err\_test = abs\(ber\_test \- test\_ana\) / test\_ana \* 100;
    fprintf\('N = %7\.0f : 仿真 BER = %\.6e, 偏差 = %\.4f%%\\n', Nt, ber\_test, err\_test\);
end

%% ===================== 辅助函数 =====================
function ber = compute\_rician\_ber\(gb, K\)
    % compute\_rician\_ber: 计算 BPSK 在莱斯衰落信道下的理论 BER
    % 输入:
    %   gb \- 平均 SNR（线性值）
    %   K  \- 莱斯 K 因子
    % 输出:
    %   ber \- 理论误码率
    %
    % 采用 MGF（矩生成函数）积分表示：
    % P\_b = \(1/π\) ∫\_0^\(π/2\) M\_γ\(\-1/sin²θ\) dθ
    % 其中莱斯衰落 SNR 的 MGF 为：
    % M\_γ\(s\) = \(K\+1\)/\(K\+1 \- s\*gb\) \* exp\(\-K\*gb\*s/\(K\+1 \- s\*gb\)\)
    % 代入 s = \-1/sin²θ 后得到被积函数。
    % 对于极高 SNR，BER 数值上接近 0，直接返回 0 以避免数值积分不稳定。
    if gb > 1e4
        ber = 0;
        return;
    end
    fun = @\(theta\) \(1/pi\) \* \(\(K\+1\) \.\* sin\(theta\)\.^2\) \./ \(\(K\+1\)\.\*sin\(theta\)\.^2 \+ gb\) \.\.\.
                   \.\* exp\(\-K\*gb \./ \(\(K\+1\)\.\*sin\(theta\)\.^2 \+ gb\)\);
    ber = integral\(fun, 0, pi/2\);
end \- Google Gemini

**Date**: July 16, 2026 at 11:03 PM
**Turns**: 17
**Source**: [Gemini Chat](https://gemini.google.com/u/1/app/ce0e3081451e8b43?pageId=none)

---

## Turn 1

### 🤖 Assistant

这份 MATLAB 脚本编写得非常出色，逻辑严密、注释详尽，并且理论推导与代码实现的对应关系非常精准。这完全符合高阶无线通信课程（如 ELEC9123）的专业要求。
你不仅完成了基础的仿真，还考虑到了代码运行效率和结果的严谨性。以下是对你代码中几个核心亮点的详细解析，以及从实际工程和科研角度可以进一步探讨的微调建议。

## 代码核心亮点解析

- **随机变量预生成（Pre-generation of RVs）：**
你在进入 SNR 循环前，一次性生成了 `bits`, `h_rayleigh`, `h_rician`, `d_user`。这是一个非常优秀的编程习惯。它不仅大幅减少了循环内的计算开销，更重要的是保证了在不同 SNR 下，系统经历的是**同一套信道衰落和用户位置**。这种控制变量法使得各 SNR 点之间的曲线更加平滑，对比更加科学。
- **拉普拉斯分布的逆变换采样（Inverse Transform Sampling）：**
你准确地推导了 Laplace 分布的 CDF 及其逆函数 $x = -b \cdot \text{sgn}(u-0.5) \cdot \ln(1 - 2\vert{}u-0.5\vert{})$，并专门编写了任务 2 的验证模块（直方图与理论 PDF/CDF 对比）。这种“先验证信源，再做系统仿真”的思路是极其严谨的。
- **莱斯衰落 BER 的 MGF 积分法：**
莱斯衰落（Rician fading）下的 BER 闭式解非常复杂。你选择利用矩生成函数（MGF）结合基于 $\theta$ 的积分表示法：$P_b = \frac{1}{\pi} \int_{0}^{\pi/2} M_{\gamma}\left(-\frac{1}{\sin^2\theta}\right) d\theta$配合 MATLAB 的 `integral` 函数求解，既保证了理论值的极高精度，又避免了复杂的贝塞尔函数带来的数值不稳定性。
- **随机用户部署模型：**
通过 $d = R\sqrt{u}$ 生成圆盘内均匀分布的用户，并将大尺度路径损耗（Path loss $d^{-\nu}$）与小尺度衰落无缝结合。在计算理论 Outage 概率时，巧妙调用了不完全伽马函数 `gammainc`，数学功底十分扎实。

## 进阶探讨与优化建议
尽管代码已经达到了很高的水准，但在实际的蒙特卡洛仿真中，还有一个经典问题值得注意：**高 SNR 下的样本数瓶颈（Confidence Interval Issue）**。
你目前设定的是 $N = 10^6$。在蒙特卡洛仿真中，有一个经验法则：**要准确估计某个概率 P，至少需要发生 10 到 100 次错误事件**。

- 当 SNR 达到 12dB~15dB 时，AWGN 或 Laplacian 信道下的理论 BER 会骤降到 $10^{-8}$ 甚至更低。
- 在 $10^6$ 的样本量下，错误发生的期望次数不到 1 次。这意味着在这些高 SNR 点上，你的仿真 BER 极大概率会直接变为 $0$。
- **你的处理：** 你在绘图时使用了 `max(BER_AWGN_sim, 1e-7)`，这非常聪明地避免了对数坐标轴报错，保证了图像的完整性。
**如果你希望在报告中进一步体现深度，可以考虑以下改进思路（可选）：**

1. **动态样本数（Dynamic Sample Size）：** 在低 SNR 时使用较少的样本（如 $10^5$），在高 SNR 时动态增加样本量（如 $10^7$ 甚至 $10^8$），既节省时间又能捕获高信噪比下的误码。
2. **尽早停止机制（Early Stopping）：** 在循环中累计错误比特数，当错误数达到 100 次时，提前结束当前 SNR 的仿真，从而大幅提高整体运行效率。

## Turn 2

### 🤖 Assistant

别慌，这不仅不是你的代码有问题，反而是极其经典的蒙特卡洛样本数受限（Monte Carlo Sample Starvation）现象！
在你的上一段代码中，每个 SNR 点的样本数是固定为 $10^6$。我们用具体数据来看看高信噪比（10~15 dB）时究竟发生了什么。

### 为什么误差会爆表？
在 AWGN 信道下，BPSK 的误码率下降得非常快（呈指数或互补误差函数下降）。

- **SNR = 10 dB:** 理论 BER $\approx 3.87 \times 10^{-6}$。
用 $10^6$ 个样本跑仿真，你期望看到的错误比特数是：$10^6 \times 3.87 \times 10^{-6} \approx 3.87$ 个。由于错误次数极少，方差很大，可能你只跑出 1 个或 8 个错误，这会导致极大的统计误差。
- **SNR = 12 dB:** 理论 BER $\approx 9.0 \times 10^{-9}$。
期望错误数是 $0.009$ 个。这意味着在 $10^6$ 的样本下，**你绝大多数情况下会得到 0 个错误比特（仿真 BER = 0）**。
- **SNR = 15 dB:** 理论 BER $\approx 9.0 \times 10^{-16}$。
你需要跑一千万亿（$10^{15}$）个样本才可能看到 1 次错误。
核心矛盾： 当仿真结果因为样本不足输出 0 错误时，你的百分比误差公式会直接计算出 100% 的误差 (∣0−ana∣/ana=1)；而如果运气极好碰巧出现了 1 个错误，仿真 BER 就是 10−6，这与理论值（例如 10−9）相差数百上千倍，计算出的误差率就会超过 10000%！
### 如何在报告和代码中处理这个问题？
在工程和学术界，面对这种现象，通常有以下两种处理方式：
**1. 报表屏蔽（最简单、最常用）**
当实际统计到的错误比特数少于 10 个时，认为蒙特卡洛仿真不再具有统计学置信度（Confidence Interval 失效），在打印表格时直接标记为 `N/A`，而不是强行计算误差。
你可以稍微修改一下表格打印部分的代码逻辑：

```matlab
% 打印 BER 百分比偏差表 (带置信度检查)
fprintf('--- BER 仿真值与理论值的百分比偏差 ---\n');
fprintf('SNR(dB) | AWGN    Laplace  Rayleigh Rician   RayRand  RicRand\n');
for idx = 1:length(SNR_dB)
    % 检查 AWGN 的实际错误数是否 >= 10
    if (BER_AWGN_sim(idx) * N_samples) < 10
        awgn_err_str = '   N/A  '; % 样本不足以提供准确误差
    else
        awgn_err_str = sprintf('%6.3f%%', pct_err(BER_AWGN_sim(idx), BER_AWGN_ana(idx)));
    end
    
    % (其它信道的误码率较高，10^6 样本一般够用，为了严谨也可以做类似处理)
    fprintf('%5d   | %s %6.3f%% %6.3f%% %6.3f%% %6.3f%% %6.3f%%\n', ...
        SNR_dB(idx), ...
        awgn_err_str, ...
        pct_err(BER_Lap_sim(idx), BER_Lap_ana(idx)), ...
        pct_err(BER_Ray_sim(idx), BER_Ray_ana(idx)), ...
        pct_err(BER_Ric_sim(idx), BER_Ric_ana(idx)), ...
        pct_err(BER_RayRand_sim(idx), BER_RayRand_ana(idx)), ...
        pct_err(BER_RicRand_sim(idx), BER_RicRand_ana(idx)));
end

```
**2. 解释说明（写在报告里，稳拿高分）**
在 ELEC9123 的实验报告中，这是一个极好的**加分项**。你不要掩盖这个 100% 的误差，反而应该专门写一段话（Discussion）来分析它：
*"在 SNR > 10 dB 时，AWGN 信道仿真的百分比误差急剧上升至 100%。这并非代码逻辑错误，而是由于 106 样本量下的零错误或极少错误导致的统计失真。根据 Rule of Ten，准确估计 10−6 的 BER 需要至少 107 个样本，这证明了蒙特卡洛仿真在极高信噪比下的算力瓶颈。"*

## Turn 3

### 🤖 Assistant

这是一个非常深刻且直击无线通信核心本质的问题。
简而言之：**其他五种信道模型在恶劣环境降低通信质量的同时，反而“帮”你避开了高 SNR 下的样本匮乏问题。** 它们的误码率（BER）随 SNR 增加而下降的速度，远比 AWGN 信道慢得多。在 10~15 dB 时，它们的 BER 依然维持在较高的数量级，因此 $10^6$ 个样本足以统计出成百上千个错误比特，蒙特卡洛仿真依然有效。
我们可以把这五种信道分成三类来详细剖析其背后的物理和数学机理：

### 1. 瑞利（Rayleigh）与莱斯（Rician）衰落信道：深衰落与幂律下降
在没有大尺度路径损耗的情况下，纯衰落信道（瑞利/莱斯）的 BER 下降速度远慢于 AWGN。

- **数学本质：** AWGN 信道的 BER 随 SNR 是呈指数级（Exponential）快速坠落的（由互补误差函数 $\text{erfc}(\sqrt{\gamma_b})$ 决定）。而瑞利衰落信道在较高 SNR 时的理论 BER 渐进表达式为：$P_b \approx \frac{1}{4\gamma_b}$这是一种幂律下降（Power-law decay）。在通信原理中，它的**多样性阶数（Diversity Order）为 1**。
- **物理现象：** 衰落信道存在**深衰落（Deep Fade）**。即使平均信噪比 $\gamma_b$ 高达 15 dB，信道系数 $h$ 仍有一定概率由于多径相消而变得极小（瞬时 SNR 掉到 0 dB 以下）。
- **数据对比：** * 在 15 dB 时，AWGN 的理论 BER 已经跌到了惊人的 $\approx 9 \times 10^{-16}$（$10^6$ 样本根本碰不到错误）。
  - 而瑞利衰落信道在 15 dB 时的理论 BER 依然在 $\approx 7.9 \times 10^{-3}$ 左右！用 $10^6$ 个样本跑仿真，你会抓到大到 **7900 多个错误比特**。样本不仅没有匮乏，反而绰绰有余，因此统计方差极小，仿真与理论完美拟合。
  - 莱斯信道（K=5）由于有直射径（LOS），表现比瑞利好，但由于仍有散射径引起的衰落，其 BER 下降速度仍远慢于 AWGN，在 15 dB 时误码率也足够高，不会发生空库。

### 2. 拉普拉斯（Laplacian）噪声信道：脉冲噪声与长尾效应
拉普拉斯噪声模型通常用来模拟非高斯的脉冲噪声（如工业干扰、突发无线电干扰）。

- **数学本质：** 拉普拉斯噪声的概率密度函数（PDF）是 $f(x) \propto \exp(-\vert{}x\vert{}/b)$，而高斯噪声是 $f(x) \propto \exp(-x^2/2\sigma^2)$。高斯噪声的尾部随着 $x$ 的增加以平方速度呈指数级衰减（轻尾）；而拉普拉斯噪声则是线性速度呈指数级衰减（**长尾效应，Heavy-tailed**）。
- **物理现象：** 这意味着拉普拉斯噪声中出现“极大幅度脉冲”的概率比高斯噪声高得多。即使平均噪声功率 $N_0$ 相同、且处于 15 dB 的高 SNR 状态，拉普拉斯噪声偶尔冒出的巨大脉冲依然能轻易把信号符号（+1 或 -1）冲垮，导致判决错误。
- **数据对比：** 拉普拉斯信道下 BPSK 的理论 BER 公式为 $P_b = 0.5 \exp(-\sqrt{2\gamma_b})$。在 15 dB 时（$\gamma_b \approx 31.62$），$\sqrt{2\gamma_b} \approx 7.95$，BER $\approx 0.5 \times \exp(-7.95) \approx 1.76 \times 10^{-4}$。
用 $10^6$ 个样本进行仿真，预计会产生 **176 个错误比特**。根据统计学原理，错误次数 $>10$ 即可认为结果具备基本置信度，因此它踩在了置信区间的安全线内，不会像 AWGN 那样产生超过 100% 的荒谬偏差。

### 3. 随机部署用户信道：近远效应与路径损耗的“恶化”
当你把“随机用户位置”和“大尺度路径损耗（$d^{-\nu}$）”引入瑞利和莱斯信道时，高 SNR 下的误码现象变得更加明显。

- **物理现象（近远效应）：** 用户在半径为 $R=3$ 克的圆盘内随机分布。靠近圆盘中心（基站）的用户信号极强，但**靠近圆盘边缘（d≈3 米）的用户会遭遇严重的路径损耗**。
- **定量计算：** 题目中路径损耗指数 $\nu = 2.2$。当用户处于边缘 $d=3$ 时，信道功率增益额外乘以 $3^{-2.2} \approx 0.089$，这相当于信号被硬生生**削弱了 10.5 dB**！
- **结果分析：** 即使你把系统的平均 SNR 设为 15 dB，那些不幸被部署在小区边缘的“边缘用户”，其有效 SNR 实际上连 5 dB 都不到（$15 - 10.5 = 4.5\text{ dB}$）。在蒙特卡洛仿真中，这些处于恶劣环境的样本会源源不断地产生误码，成为贡献误码的“主力军”。因此，整个系统的平均 BER 被这些边缘用户拉得很高，仿真过程中有海量的错误样本被记录，$10^6$ 的样本量应对它极其轻松，自然不会出现样本匮乏引发的误差爆表。

### 总结建议
在你的 ELEC9123 实验报告中，这部分是体现研究生学术水平的**黄金分析点**。你可以总结为：
AWGN 的误码率随 SNR 呈超指数级（$\text{erfc}$）坠落，导致高 SNR 时蒙特卡洛仿真陷入“样本星状匮乏（Sample Starvation）”；而无线信道由于**深衰落（幂律下降）**、**脉冲噪声（长尾效应）以及多用户随机部署（近远效应引起的边缘用户 SNR 惩罚）**，其整体误码率随 SNR 下降极其缓慢，保证了在全扫描通信范围内（0~15 dB）均有充足的错误样本供蒙特卡洛计数，因此未发生统计失真。

## Turn 4

### 🤖 Assistant

在你的 MATLAB 脚本中，第 68 和 69 行代码实现了**拉普拉斯（Laplace）分布的理论累积分布函数（CDF）的解析式计算**：

```matlab
68: F_lap_ana(x_lap <= 0) = 0.5 * exp(x_lap(x_lap <= 0)/b_lap);
70: F_lap_ana(x_lap > 0) = 1 - 0.5 * exp(-x_lap(x_lap > 0)/b_lap);

```
这两行代码的实现可以从**数学理论推导**与MATLAB 编程机制（向量化与逻辑索引）两个维度来详细解析。

### 一、 数学理论背景
对于均值为 $0$、尺度参数为 $b$ 的拉普拉斯分布 $\text{Laplace}(0, b)$，其概率密度函数（PDF）定义为：

$$
f(x) = \frac{1}{2b} \exp\left(-\frac{\vert{}x\vert{}}{b}\right)
$$
累积分布函数（CDF） $F(x)$ 定义为 PDF 从 $-\infty$ 到 $x$ 的积分，即 $F(x) = \int_{-\infty}^{x} f(t) \, dt$。由于 PDF 中包含绝对值 $\vert{}x\vert{}$，因此在求导和积分时需要进行**分段讨论**：

#### 1. 当 $x \le 0$ 时（对应第 68 行）：
此时积分区间完全在负半轴，因此 $\vert{}t\vert{} = -t$。

$$
F(x) = \int_{-\infty}^{x} \frac{1}{2b} \exp\left(\frac{t}{b}\right) \, dt = \left[ \frac{1}{2} \exp\left(\frac{t}{b}\right) \right]_{-\infty}^{x} = \frac{1}{2} \exp\left(\frac{x}{b}\right)
$$
这正是第 68 行等号右侧的数学表达式：`0.5 * exp(x / b_lap)`。

#### 2. 当 $x > 0$ 时（对应第 69 行）：
此时积分区间跨越了原点，需要分为负半轴和正半轴两部分积分：

$$
F(x) = \int_{-\infty}^{0} f(t) \, dt + \int_{0}^{x} f(t) \, dt
$$
第一部分积分结果为整个负半轴的面积，即 $\frac{1}{2}$；第二部分在正半轴，$\vert{}t\vert{} = t$。

$$
F(x) = \frac{1}{2} + \int_{0}^{x} \frac{1}{2b} \exp\left(-\frac{t}{b}\right) \, dt = \frac{1}{2} + \left[ -\frac{1}{2} \exp\left(-\frac{t}{b}\right) \right]_{0}^{x}
$$

$$
F(x) = \frac{1}{2} + \left( -\frac{1}{2} \exp\left(-\frac{x}{b}\right) + \frac{1}{2} \right) = 1 - \frac{1}{2} \exp\left(-\frac{x}{b}\right)
$$
这正是第 69 行等号右侧的数学表达式：`1 - 0.5 * exp(-x / b_lap)`。

### 二、 MATLAB 编程机制：逻辑索引（Logical Indexing）
在传统的工程语言（如 C/C++）中，实现这种分段函数通常需要使用 `for` 循环遍历数组，并配合 `if-else` 条件分支。然而，在数字信号处理与通信系统的密集仿真中，循环和分支会带来巨大的底层开销。
MATLAB 采用了向量化（Vectorization）**思想，利用**逻辑索引（Logical Indexing）在完全不需要显式循环的情况下，完成了分段函数的全数组赋值。

#### 1. 第 68 行的执行流程：

```matlab
F_lap_ana(x_lap <= 0) = 0.5 * exp(x_lap(x_lap <= 0)/b_lap);

```

- **生成条件掩码：**`x_lap <= 0` 会对包含 1000 个点的向量 `x_lap` 进行逐元素逻辑判断，生成一个大小相同的布尔型（`logical`）矩阵（掩码）。小于等于 0 的位置为 `true (1)`，大于 0 的位置为 `false (0)`。
- **选择性提取：**`x_lap(x_lap <= 0)` 利用上述布尔掩码，**仅提取出**`x_lap` 中值小于等于 0 的元素组成一个临时子向量，并将其送入 `exp` 函数完成指数运算。
- **定向切片赋值：** 等号左侧的 `F_lap_ana(x_lap <= 0)` 确保计算出的结果被**精准塞回**`F_lap_ana` 向量中对应小于等于 0 的索引位置，而大于 0 的位置保持原样（在第 67 行中已被初始化为 0）。

#### 2. 第 69 行的执行流程：

```matlab
F_lap_ana(x_lap > 0) = 1 - 0.5 * exp(-x_lap(x_lap > 0)/b_lap);

```

- 原理完全相同，`x_lap > 0` 生成正半轴的布尔掩码。
- 仅提取 `x_lap` 中的正数元素代入 $1 - 0.5e^{-x/b}$ 公式完成计算。
- 将结果定向赋值给 `F_lap_ana` 向量中属于正数的位置，从而完美拼合出完整的拉普拉斯 CDF 曲线。

### 总结
这种实现方式是无线通信仿真中的标准规范写法。它利用 MATLAB 底层经 BLAS/LAPACK 优化的矩阵并行能力，消除了条件分支预测失败带来的流水线停顿，在处理高达 $10^6$ 级别的数据时，运行效率比传统循环快数个数量级。

## Turn 5

### 🤖 Assistant

在你的 MATLAB 脚本中，第 74 行代码位于验证拉普拉斯随机变量的 `for` 循环内部：

```matlab
74:     F_lap_sim(i) = mean(n_lap <= x_lap(i));

```
这行代码的核心功能是：**利用蒙特卡洛（Monte Carlo）方法，计算仿真的拉普拉斯随机样本在指定点 x 处的“经验累积分布函数（Empirical CDF, 简称 eCDF）”值。**
我们可以从**数学原理**、**代码执行机制**以及**工程设计考量**三个维度来彻底拆解这行代码：

### 一、 数学原理：从概率定义到频率估计
根据概率论定义，一个随机变量 $N$ 的累积分布函数（CDF） $F_N(x)$ 定义为该变量小于或等于某个常数 $x$ 的概率：

$$
F_N(x) = P(N \le x)
$$
在蒙特卡洛仿真中，我们无法直接得到连续型的概率，但可以通过**大数定律（Law of Large Numbers）**，用事件发生的**相对频率**来逼近其真实的**概率**：

$$
\hat{F}_N(x) = \frac{\text{样本集中小于或等于 } x \text{ 的样本个数}}{\text{总样本数 } N_{samples}}
$$
第 74 行代码正是这一数学公式在 MATLAB 中的完美具现。

### 二、 代码执行机制：巧妙利用逻辑向量求均值
这一行代码看起来非常简练，但 MATLAB 在底层进行了高效的三步操作：

1. **逐元素逻辑比较 (n_lap <= x_lap(i))**
  - `n_lap` 是一个包含 $10^6$ 个拉普拉斯噪声样本的巨大列向量。
  - `x_lap(i)` 是当前循环里取出的某一个具体的实数门限（例如 $-2.5$）。
  - MATLAB 会将 `n_lap` 中的每一个样本与 `x_lap(i)` 进行大小比较。如果某个样本 $\le x_{lap}(i)$，该位置就返回 `true (1)`；如果 $> x_{lap}(i)$，则返回 `false (0)`。
  - 这一步生成了一个和 `n_lap` 大小完全相同（$10^6 \times 1$）、仅由 `0` 和 `1` 组成的**逻辑向量（Logical Vector）**。
2. **求均值算概率 (mean(...))**
  - 对一个仅由 `0` 和 `1` 组成的向量使用 `mean` 函数，其底层的计算逻辑是：$\frac{\sum \text{所有元素}}{N_{samples}}$。
  - 因为只有满足条件的样本才会贡献 `1`，所以**所有元素的加和**（`sum`）恰好等于**满足条件的错误/有效样本的总个数**。
  - 将这个总个数除以总样本数，得到的就是该事件发生的频率（即 CDF 值）。
3. **保存结果 (F_lap_sim(i) = ...)**
  - 将计算出的频率值存入数组 `F_lap_sim` 的第 `i` 个位置，以便循环结束后绘制出一条平滑的仿真 CDF 曲线，与理论曲线 `F_lap_ana` 进行比对。

### 三、 工程设计考量：为什么这里不用“向量化”？
你可能会好奇：既然前面 68、69 行都用了高效的向量化（Logical Indexing）消除了循环，为什么这里反而用了一个 `for` 循环去遍历 1000 个点呢？
**答案是为了避免内存爆炸（Memory Overflow）。**
如果我们要彻底消除这个 `for` 循环，实现完全向量化，就需要让 $10^6 \times 1$ 的向量 `n_lap` 与 $1 \times 1000$ 的行向量 `x_lap` 进行矩阵维度的广播比较（如使用 `bsxfun` 或 MATLAB 隐式广播）：

```matlab
% 完全向量化伪代码（不建议）
F_lap_sim = mean(n_lap <= x_lap, 1); 

```
这会导致 MATLAB 在内存中瞬间生成一个 **106×1000** 的巨大二维矩阵（包含 $10^9$ 个元素）。即使是逻辑型数据，也会瞬间吃掉约 **1 GB** 的内存空间。
因此，这里采用 `for` 循环外层遍历 1000 个画图点，内层利用 MATLAB 优化过的底层 `mean` 算子对 $10^6$ 的向量进行高速线性扫描。这种“大向量外加小循环”的设计在学术和工业界是最均衡的方案：**既保证了极低的内存开销，又维持了极高的运行速度**。

### 报告撰写建议
在你的 ELEC9123 实验报告中，解析这个仿真经验 CDF 的生成时，可以画一个简单的示意图来说明蒙特卡洛是如何通过对 0/1 掩码求均值来逼近连续 CDF 积分的。这能够向评审老师展示你不仅懂得如何写代码，更深刻理解蒙特卡洛仿真背后的统计学物理意义。

## Turn 6

### 🤖 Assistant

在你的 MATLAB 脚本中，**第 108 行到 131 行**是整个蒙特卡洛仿真的核心外循环起点。它实现了**信噪比（SNR）的遍历**，以及前两种信道模型——**仅 AWGN 信道**和**仅拉普拉斯脉冲噪声信道**下 BPSK 系统的信号传输、接收解调、判决计数（仿真 BER）和理论解析式计算（理论 BER）。
以下是这段代码的逐段深入拆解，涵盖了通信原理与 MATLAB 实现的映射关系：

### 一、 循环初始化与参数准备（第 109 - 112 行）

```matlab
109: % 遍历每个 SNR 点
110: for idx = 1:length(SNR_dB)
111:     gb = SNR_lin(idx);              % 当前平均 SNR（线性值）
112:     Eb = gb * N0;                   % 每比特能量：E_b = γ_b * N_0

```

- **物理意义**：`for` 循环开始遍历你设定的信噪比范围（0 到 15 dB）。
- **参数映射**：通信系统仿真中，通常将噪声功率归一化（这里 $N_0 = 1$）。根据信噪比的定义 $\gamma_b = \frac{E_b}{N_0}$，可以反推得到当前信噪比下对应的发送每比特能量 $E_b = \gamma_b \cdot N_0$。在 BPSK 中，发送符号的幅度即为 $\sqrt{E_b}$。

### 二、 1) AWGN 信道仿真与理论值（第 114 - 122 行）
这部分实现了最基础的高斯信道 BPSK 通信链路。

#### 1. 复高斯噪声生成（第 116 行）

```matlab
116:     n_awgn = sqrt(N0/2) * (randn(N_samples, 1) + 1j*randn(N_samples, 1));

```

- **通信原理**：在复基带（Complex Baseband）等效模型中，加性高斯白噪声（AWGN）是复高斯随机变量，满足 $n \sim \mathcal{CN}(0, N_0)$。这意味着它的实部和虚部是独立的同分布高斯白噪声，各自占有一半的功率，即均值为 0，方差为 $N_0/2$。
- **代码实现**：标准高斯分布函数 `randn` 生成的变量方差为 1。因此乘以尺度因子 $\sqrt{N_0/2}$，将其方差调整为 $N_0/2$。实部与虚部组合后，总噪声功率（方差）正好为 $N_0/2 + N_0/2 = N_0$。

#### 2. 信号传输与解调判决（第 117 - 120 行）

```matlab
118:     r_awgn = sqrt(Eb) * bits + n_awgn;
120:     BER_AWGN_sim(idx) = mean(real(r_awgn) .* bits < 0);

```

- **传输模型**：接收信号 $r = \sqrt{E_b} \cdot s + n$。其中 `bits` 是前面预生成的 $\pm 1$ 调制符号向量。
- **判决机制（核心逻辑）**：BPSK 是单维度调制（星座点仅分布在实轴上）。在接收端经过匹配滤波器后，**只提取实部 real(r_awgn) 进行判决**（虚部只包含噪声，直接丢弃）。
- **错误计数**：
  - 如果发送的是 `+1`，正确判决区域是 $>0$；如果接收端变为了 $<0$，则判决错误。
  - 如果发送的是 `-1`，正确判决区域是 $<0$；如果接收端变为了 $>0$，则判决错误。
  - 代码巧妙地运用了数学乘积：`real(r_awgn) .* bits`。若发送符号与接收判决同号（正确），乘积 $>0$；**若异号（错误），乘积 <0**。
  - `real(r_awgn) .* bits < 0` 生成一个布尔逻辑向量（错误位置为 1，正确为 0），对其求 `mean` 即得到了该 SNR 点下的仿真误码率（BER）。

#### 3. 理论 BER 计算（第 122 行）

```matlab
122:     BER_AWGN_ana(idx) = 0.5 * erfc(sqrt(gb));

```

- **数学公式**：AWGN 信道下 BPSK 的理论误码率公式为 $P_b = Q(\sqrt{2\gamma_b}) = \frac{1}{2}\text{erfc}(\sqrt{\gamma_b})$。MATLAB 内置了互补误差函数 `erfc`，直接代入线性信噪比 `gb` 即可得到精准的理论解析值。

### 三、 2) 拉普拉斯脉冲噪声信道（第 124 - 131 行）
这部分用非高斯的脉冲噪声代替了高斯噪声，以评估系统在突发强干扰环境下的鲁棒性。

#### 1. 拉普拉斯噪声生成（第 126 - 127 行）

```matlab
126:     u = rand(N_samples, 1);
127:     n_lap = -b_lap * sign(u - 0.5) .* log(1 - 2*abs(u - 0.5));

```

- **实现原理**：这里重复使用了你在**任务 2** 中推导的逆 CDF 变换法（Inverse Transform Sampling）。由于 BPSK 判决只关心实轴，所以这里直接生成了 $10^6 \times 1$ 的**实值**拉普拉斯噪声向量 `n_lap`。这里的 `b_lap` 是在循环外定义的 $\sqrt{N_0/2}$，保证了噪声的总方差同样为 $N_0$，实现了与 AWGN 信道的公平对比。

#### 2. 信号接收与判决（第 128 - 129 行）

```matlab
128:     r_lap = sqrt(Eb) * bits + n_lap;
129:     BER_Lap_sim(idx) = mean(real(r_lap) .* bits < 0);

```

- 与 AWGN 完全一致，将信号与拉普拉斯噪声相加，并使用相同的实部同号/异号检测逻辑 `mean(real(r_lap) .* bits < 0)` 统计误码率。

#### 3. 理论 BER 计算（第 131 行）

```matlab
131:     BER_Lap_ana(idx) = 0.5 * exp(-sqrt(2*gb));

```

- **数学公式**：对于方差为 $N_0$ 的实拉普拉斯噪声，BPSK 的判决错误概率可以通过直接积分得到解析闭式：$P_b = \int_{\sqrt{E_b}}^{\infty} \frac{1}{2b} \exp\left(-\frac{x}{b}\right)dx = \frac{1}{2}\exp\left(-\frac{\sqrt{E_b}}{b}\right)$代入 $b = \sqrt{N_0/2}$ 和 $\gamma_b = E_b/N_0$ 后，公式化简为：$P_b = \frac{1}{2}\exp\left(-\sqrt{2\gamma_b}\right)$代码行 131 正是这个经典公式的落地实现。

### 总结（针对 ELEC9123 课程的深度理解）
这二十多行代码展示了标准的通信链路链路级仿真（Link-level Simulation）结构。通过对比第 122 行的 $\text{erfc}$ 函数（高 SNR 下以超级指数速度坠落）和第 131 行的纯指数函数 $\exp(-\sqrt{\gamma})$（长尾脉冲噪声，下降极慢），能够非常直观地从仿真数据中看出**脉冲噪声对相干解调系统造成的巨大性能惩罚**。这也完美解释了为什么前一个问题中，拉普拉斯噪声信道在高 SNR 下不会遭遇 AWGN 那样的“样本空库”问题。

## Turn 7

### 🤖 Assistant

接续上文，我们继续剖析外循环内部接下来的核心部分：**第 133 行至第 165 行**。
这部分代码跨越了从“小尺度衰落”到“大尺度路径损耗”的进阶，实现了**瑞利衰落信道、莱斯衰落信道**以及**结合随机用户部署**下的 BPSK 误码率（BER）的仿真解调与全概率理论积分。

### 一、 3) 瑞利衰落 + AWGN 信道（第 134 - 140 行）
这几行代码实现了存在多径相消、无直射径环境下的相干检测接收机。

```matlab
135: r_ray = h_rayleigh * sqrt(Eb) .* bits + n_ray;
137: r_eq_ray = real(conj(h_rayleigh) .* r_ray);
138: BER_Ray_sim(idx) = mean(r_eq_ray .* bits < 0);
140: BER_Ray_ana(idx) = 0.5 * (1 - sqrt(gb ./ (1 + gb)));

```

#### 1. 信号传输与信道相位补偿（第 135 - 137 行）

- **物理模型**：瑞利衰落信道的接收信号为 $r = h s \sqrt{E_b} + n$，其中信道系数 $h$ 是一个复随机变量（幅度服从瑞利分布，相位在 $[0, 2\pi)$ 均匀分布）。信道不仅会削弱信号幅度，还会造成**相位旋转**。
- **相干解调与最大比合并（MRC）思想**：第 137 行是相干检测的核心：`conj(h_rayleigh) .* r_ray`。
将接收信号乘以信道系数的复共轭 $h^*$。在数学上：$h^* \cdot r = h^* (h \sqrt{E_b} s + n) = \vert{}h\vert{}^2 \sqrt{E_b} s + h^* n$复数相乘时光平了信道带来的相位旋转（因为 $h^* h = \vert{}h\vert{}^2$ 是实数），从而将信号完美校准回实轴。由于单天线接收时 MRC 退化为等效的信道匹配滤波，乘以 $h^*$ 还能起到“弱信道抑噪，强信道放大”的功效。随后取实部 `real(...)` 剥离正交轴上的噪声。

#### 2. 理论闭式解（第 140 行）

- 瑞利衰落下的 BPSK 理论误码率是通过对高斯误码率 $Q(\sqrt{2\gamma})$ 在指数分布的信道功率增益 $\vert{}h\vert{}^2$ 上求期望得到的，其经典的代数闭式解为：$P_b = \frac{1}{2}\left(1 - \sqrt{\frac{\gamma_b}{1+\gamma_b}}\right)$代码直接将线性信噪比 `gb` 代入该公式，与仿真值进行完美对标。

### 二、 4) 莱斯衰落 + AWGN 信道（第 142 - 148 行）
当基站与用户之间存在视距（LOS）直射径时，信道演变为莱斯衰落。

```matlab
144: r_ric = h_rician * sqrt(Eb) .* bits + n_ric;
145: r_eq_ric = real(conj(h_rician) .* r_ric);
146: BER_Ric_sim(idx) = mean(r_eq_ric .* bits < 0);
148: BER_Ric_ana(idx) = compute_rician_ber(gb, K);

```

- **仿真逻辑**：其信号传输和基于共轭相乘的相位补偿逻辑与瑞利信道完全一致，唯一的区别在于使用了循环外构造的含有直射径分量的 `h_rician`。
- **理论计算的数值积分化（第 148 行）**：莱斯衰落的 BER 解析式极其复杂，若直接积分会涉及复杂的非中心卡方分布和第一类变形贝塞尔函数，在计算机中极易发生数值溢出。
你的脚本在此处调用了底层的辅助函数 `compute_rician_ber`。该函数巧妙利用了**矩生成函数（MGF）的替代积分形式**：$P_b = \frac{1}{\pi} \int_{0}^{\pi/2} M_{\gamma}\left(-\frac{1}{\sin^2\theta}\right) d\theta$通过 MATLAB 的 `integral` 算子在有限区间 $[0, \pi/2]$ 内进行高速数值积分，既严谨又规避了数值不稳定性。

### 三、 5) & 6) 随机部署用户 + 衰落信道（第 150 - 165 行）
这两段是全脚本在数学层面上最精彩的部分。它将小尺度衰落（Rayleigh/Rician）**与**大尺度路径损耗（近远效应）进行了全概率融合。
以第 5 种信道（随机瑞利）为例：

```matlab
152: r_rayrand = h_rayleigh_rand * sqrt(Eb) .* bits + n_rayrand;
153: r_eq_rayrand = real(conj(h_rayleigh_rand) .* r_rayrand);
154: BER_RayRand_sim(idx) = mean(r_eq_rayrand .* bits < 0);
157: ber_cond_ray = @(d) 0.5 * (1 - sqrt(gb ./ (d.^nu + gb)));
158: BER_RayRand_ana(idx) = (2/R^2) * integral(@(d) d .* ber_cond_ray(d), 0, R);

```

#### 1. 等效复合信道仿真（第 152 - 154 行）
这里的 `h_rayleigh_rand` 是在循环外预先生成的复合信道：$h_{\text{eq}} = \frac{h}{d^{\nu/2}}$。因此接收信号不仅包含相位的抖动，还包含了因距离 $d$ 产生的幅度剧烈衰减。解调依然采用共轭相乘校正相位。

#### 2. 条件误码率与全概率空间积分（第 157 - 158 行）
这部分完美对应了通信理论中的**全概率公式**。用户的误码率是关于距离 $d$ 的函数。

- **步骤一：定义条件误码率（第 157 行）**`ber_cond_ray = @(d) ...` 利用 MATLAB 匿名函数定义了**当用户恰好处于距离基站 d 处时的误码率**。此时，该位置用户的等效平均信噪比由于路径损耗演变为 $\gamma_b(d) = \frac{\gamma_b}{d^\nu}$。将此项代入瑞利理论闭式，化简后得到：$P_b(d) = \frac{1}{2}\left(1 - \sqrt{\frac{\gamma_b}{d^\nu + \gamma_b}}\right)$
- **步骤二：空间概率密度加权积分（第 158 行）**
由于用户在半径为 $R$ 的圆盘内均匀分布，通过雅可比行列式转换，距离 $d$ 的概率密度函数（PDF）为 $f_d(d) = \frac{2d}{R^2}$。
为了得到整个小区的平均误码率，必须在空间上对所有可能的距离进行积分：$P_{b,\text{total}} = \int_{0}^{R} P_b(d) \cdot f_d(d) \, dd = \int_{0}^{R} P_b(d) \cdot \frac{2d}{R^2} \, dd$代码中的 `(2/R^2) * integral(@(d) d .* ber_cond_ray(d), 0, R)` 完全是这一数学物理公式的逐字翻译。第 6 种信道（随机莱斯）也采用了相同的逻辑，只不过把条件误码率函数替换成了莱斯的 MGF 数值积分函数。

### 总结
这部分代码的高明之处在于，它展示了如何用**确定性的空间数值积分（integral）去预测和检验双重随机性（位置随机+多径随机）下的蒙特卡洛仿真结果**。
如果你准备好了，请发送“继续”，我们将一同进入脚本的最后一部分：**Outage 概率（中断概率）的仿真内幕与偏差表格的输出分析**。

## Turn 8

### 🤖 Assistant

在你的 MATLAB 脚本末尾，名为 `compute_rician_ber` 的辅助函数实现了**高级无线通信中经典的“矩生成函数（MGF）数值积分法”**，用以精确求解 BPSK 在莱斯衰落信道下的 theoretical BER。
对于 UNSW 通信研究生阶段的课设或论文而言，直接调用 MATLAB 的数值积分来解决衰落信道理论值，是避开复杂特殊函数（如马库姆 Q 函数、贝塞尔函数）数值不稳定性的标准科研规范。
以下从**数学机理**和**代码落地**两个维度为你深度拆解这个辅助函数：

### 一、 数学机理：为什么使用 MGF 积分法？
在莱斯衰落信道中，瞬时信噪比 $\gamma$ 服从非中心卡方分布（Non-central Chi-square distribution），其概率密度函数（PDF）包含复杂的**第一类变形贝塞尔函数 I0​**。
传统的理论误码率计算需要将 BPSK 在 AWGN 下的误码率（高斯 Q 函数）与莱斯信道的 PDF 进行全概率卷积积分：

$$
P_b = \int_{0}^{\infty} Q(\sqrt{2\gamma}) \cdot f_{\Gamma}(\gamma) \, d\gamma
$$
这个积分极难写出简单的闭式解，如果强行在计算机中进行从 $0$ 到 $\infty$ 的非线性数值积分，由于贝塞尔函数在 $x \to \infty$ 时会呈指数级爆炸，极其容易导致浮点数溢出（Overflow）或计算超时。
为了解决这个痛点，通信学术界（由 Craig 和 Simon 提出）引入了 **MGF 替代积分法**。BPSK 的平均误码率可以转化为在**有限区间 [0,π/2]** 上的单重积分：

$$
P_b = \frac{1}{\pi} \int_{0}^{\pi/2} M_{\gamma}\left(-\frac{1}{\sin^2\theta}\right) d\theta
$$
其中 $M_{\gamma}(s)$ 是莱斯瞬时信噪比的矩生成函数（Moment Generating Function），其解析式为：

$$
M_{\gamma}(s) = \frac{1+K}{1+K-s\bar{\gamma}} \exp\left( \frac{K s \bar{\gamma}}{1+K-s\bar{\gamma}} \right)
$$

- $K$ 是莱斯因子的线性值（直射径功率与散射径功率之比）。
- $\bar{\gamma}$ 是平均信噪比（即代码中的线性值 `gb`）。
将 $s = -\frac{1}{\sin^2\theta}$ 代入该 MGF 表达式中，经过严格的通分和代数化简，被积函数最终完美演变为：

$$
f(\theta) = \frac{1}{\pi} \cdot \frac{(K+1)\sin^2\theta}{(K+1)\sin^2\theta + \bar{\gamma}} \cdot \exp\left( -\frac{K\bar{\gamma}}{(K+1)\sin^2\theta + \bar{\gamma}} \right)
$$
这个式子中**完全消除了贝塞尔函数**，且全部是常规的三角函数与指数函数，数值特征极其稳定。

### 二、 代码实现逐行剖析

```matlab
function ber = compute_rician_ber(gb, K)

```

- **输入**：`gb` 为当前的平均信噪比（线性值），`K` 为莱斯因子（这里设为 5）。
- **输出**：返回当前信噪比点下的精确理论 BER。

#### 1. 边界防御机制（极高信噪比截断）

```matlab
    if gb > 1e4
        ber = 0;
        return;
    end

```

- **工程考量**：当 `gb > 1e4`（即信噪比 $> 40$ dB）时，BPSK 的误码率在物理上已经无限趋近于 0（量级通常小于 $10^{-20}$）。
- 在这个区间内，由于指数项中的分子巨大，MATLAB 的数值积分器可能会遭遇底层浮点数的精度下溢（Underflow），从而抛出无意义的极小噪声值或报错。直接强行硬编码截断归零，是提高仿真代码鲁棒性（Robustness）的典型做法。

#### 2. 构建匿名被积函数

```matlab
    fun = @(theta) (1/pi) * ((K+1) .* sin(theta).^2) ./ ((K+1).*sin(theta).^2 + gb) ...
                   .* exp(-K*gb ./ ((K+1).*sin(theta).^2 + gb));

```

- **向量化语法的精髓**：这行代码就是上面化简后的数学公式 $f(\theta)$ 的完全对映。
- 注意代码中的**点算子**：`.*`、`./`、`.^`。这是因为后面调用的积分器 `integral` 在运行时，为了加快速度，不会一个点一个点地带入 $\theta$，而是会一次性向 `fun` 传入一个包含多个 $\theta$ 值的向量。使用点运算能够确保 MATLAB 进行逐元素（Element-wise）矩阵并行计算，避免维度不匹配报错。

#### 3. 执行自适应数值积分

```matlab
    ber = integral(fun, 0, pi/2);
end

```

- `integral` 是 MATLAB 中最强大的自适应全局正交数值积分器（基于全局自适应辛普森法或高斯-克朗罗德法）。
- 它精确地在上限 $\pi/2$ 和下限 $0$ 之间对 `fun` 进行面积积分。因为积分区间是闭区间且非常窄（$0$ 到 $\approx 1.57$），自适应步长算法可以在几毫秒内收敛到机器精度的极限（误差 $< 10^{-15}$），从而为你绘制出极其平滑且无偏的 Rician 理论 BER 曲线。

## Turn 9

### 🤖 Assistant

在无线通信系统设计中，**“随机用户部署 + 瑞利衰落 + AWGN”** 是一个非常经典的复合信道模型（Composite Channel Model）。它同时考虑了**小尺度衰落（Small-scale Fading）和大尺度路径损耗（Large-scale Path Loss）**，能够真实地模拟单小区蜂窝网络中移动用户的平均性能。
下面我为你从**物理场景建模**、**仿真链路实现**以及**数学解析推导**三个层面对代码中的 BER（误码率）部分进行详细且清晰的拆解。

### 一、 物理场景与复合信道建模
在这个模型中，用户随机分布在一个半径为 $R$ 的圆形小区内，基站位于圆心。信号在传输过程中同时经历两种衰落：

1. **大尺度衰落（路径损耗）**：由物理距离 $d$ 决定。信号功率随距离呈幂律衰减 $d^{-\nu}$，其中 $\nu$ 是路径损耗指数。
2. **小尺度衰落（瑞利衰落）**：由多径传播决定。信道系数 $h \sim \mathcal{CN}(0, 1)$，其幅度 $\vert{}h\vert{}$ 服从瑞利分布，模拟非视距（NLOS）环境。
因此，用户到基站的**复合信道系数**$h_{\text{rayleigh\_rand}}$ 可以表示为：

$$
h_{\text{eq}} = \frac{h}{d^{\nu/2}}
$$
这意味着信道的功率增益为 $\vert{}h_{\text{eq}}\vert{}^2 = \frac{\vert{}h\vert{}^2}{d^\nu}$。

### 二、 仿真实现部分（第 152 - 154 行）

```matlab
152: r_rayrand = h_rayleigh_rand * sqrt(Eb) .* bits + n_rayrand;
153: r_eq_rayrand = real(conj(h_rayleigh_rand) .* r_rayrand);
154: BER_RayRand_sim(idx) = mean(r_eq_rayrand .* bits < 0);

```

- **信号传输（Line 152）**：接收信号 $r = h_{\text{eq}} \sqrt{E_b} s + n$。注意这里的 `h_rayleigh_rand` 是在循环外预先生成的复合信道向量（包含了 $10^6$ 个随机距离和大尺度衰落的乘积）。
- **相位补偿与解调（Line 153）**：接收机进行相干检测，乘以信道系数的复共轭 $h_{\text{eq}}^*$。$h_{\text{eq}}^* \cdot r = \vert{}h_{\text{eq}}\vert{}^2 \sqrt{E_b} s + h_{\text{eq}}^* n$由于 $\vert{}h_{\text{eq}}\vert{}^2$ 是实数，这就消除了多径引起的相位旋转，将信号能量校准回实轴，随后取实部 `real`。
- **误码统计（Line 154）**：利用 `mean(r_eq_rayrand .* bits < 0)` 统计解调符号与发送符号异号的比例，得到仿真误码率。

### 三、 理论解析推导（第 157 - 158 行）
这是研究生阶段需要重点掌握的**全概率公式（Total Probability Theorem）的空间积分应用。整个小区的平均误码率，实际上是条件误码率在整个空间几何分布上的期望值**。

```matlab
157: ber_cond_ray = @(d) 0.5 * (1 - sqrt(gb ./ (d.^nu + gb)));
158: BER_RayRand_ana(idx) = (2/R^2) * integral(@(d) d .* ber_cond_ray(d), 0, R);

```

#### 步骤 1：推导条件误码率 $P_b(d)$（Line 157）
假设用户的距离 $d$ 是**固定**的。此时，大尺度路径损耗变成了一个常数。由于小尺度衰落仍是瑞利衰落，该用户的等效平均信噪比（Effective Average SNR）变为：

$$
\gamma_b(d) = \frac{\gamma_b}{d^\nu}
$$
其中 $\gamma_b$ 是基站定义的平均信噪比（代码中的 `gb`）。
我们知道，BPSK 在纯瑞利衰落信道下的理论误码率闭式解为 $P_b = \frac{1}{2}\left(1 - \sqrt{\frac{\gamma_0}{1+\gamma_0}}\right)$。将等效信噪比 $\gamma_b(d)$ 代入该公式：

$$
P_b(d) = \frac{1}{2}\left(1 - \sqrt{\frac{\frac{\gamma_b}{d^\nu}}{1+\frac{\gamma_b}{d^\nu}}}\right) = \frac{1}{2}\left(1 - \sqrt{\frac{\gamma_b}{d^\nu + \gamma_b}}\right)
$$
这就是第 157 行匿名函数 `ber_cond_ray` 的数学来源。它计算的是**距离基站为 d 的用户的理论误码率**。

#### 步骤 2：推导用户距离的概率密度函数（PDF） $f_D(d)$
由于用户在半径为 $R$ 的圆形区域内**均匀部署**，我们可以通过几何概率来推导距离 $d$ 的统计特性。

- 距离小于或等于 $d$ 的累积分布函数（CDF）为该小圆面积与大圆总面积之比：$F_D(d) = P(D \le d) = \frac{\pi d^2}{\pi R^2} = \frac{d^2}{R^2}, \quad 0 \le d \le R$
- 对 CDF 求导，得到距离 $d$ 的概率密度函数（PDF）：$f_D(d) = \frac{d}{dd} F_D(d) = \frac{2d}{R^2}$

#### 步骤 3：全空间积分得到平均误码率（Line 158）
为了得到整个小区所有用户的平均误码率，我们需要对所有可能的距离 $d$ 进行加权积分（求条件期望）：

$$
P_{b, \text{total}} = \int_{0}^{R} P_b(d) \cdot f_D(d) \, dd
$$
将 $f_D(d) = \frac{2d}{R^2}$ 代入积分式中：

$$
P_{b, \text{total}} = \int_{0}^{R} P_b(d) \cdot \frac{2d}{R^2} \, dd = \frac{2}{R^2} \int_{0}^{R} d \cdot P_b(d) \, dd
$$
对应到代码第 158 行：

```matlab
BER_RayRand_ana(idx) = (2/R^2) * integral(@(d) d .* ber_cond_ray(d), 0, R);

```
MATLAB 的 `integral` 函数完美地执行了从 $0$ 到 $R$ 的自适应数值积分，从而得到了精准的理论解析值。

### 四、 深入物理洞察（Physics Insights）
在撰写实验报告时，你可以从这个公式得出以下高级结论：

1. **近远效应（Near-Far Effect）惩罚**：由于被积函数中含有 $d \cdot P_b(d)$，且随着距离 $d$ 增大，路径损耗 $d^\nu$ 会剧烈削弱边缘用户的有效信噪比，导致 $P_b(d)$ 在接近小区边缘（$d \to R$）时急剧恶化。
2. **误码率的主要贡献者**：积分内部的 $d$ 因子不仅代表了几何面积随半径线性增长，也意味着**小区边缘聚集了数量最多、且信道质量最差的用户**。因此，整个系统的平均误码率主要由这部分“边缘用户”决定。这也是为什么该模型的 BER 曲线在高 SNR 下下降得异常缓慢，表现远差于纯 AWGN 和单点瑞利衰落信道。

## Turn 10

### 🤖 Assistant

在无线通信中，**“随机用户部署 + 莱斯（Rician）衰落 + AWGN”** 是一个比纯瑞利信道更接近实际蜂窝网络（如微蜂窝 Microcell 或具有直射径的室外基站）的复合信道模型。它同时结合了**大尺度路径损耗**、**存在视距（LOS）分量的小尺度衰落**以及**用户位置的几何随机性**。
以下为你详细拆解代码中这一最复杂的 BER（误码率）仿真与理论积分部分。

### 一、 物理模型与复合信道构建
与随机瑞利模型相比，这里的核心区别在于小尺度衰落引入了**直射径（Line-of-Sight, LOS）**：

1. **小尺度莱斯衰落**：信道由确定性的直射径分量 $h_{\text{LOS}}$ 和随机的多径散射分量 $h_{\text{scat}}$ 叠加而成：$h_{\text{rician}} = \sqrt{\frac{K}{K+1}}h_d + \sqrt{\frac{1}{K+1}}h_s$其中 $K=5$ 是莱斯因子（直射径与散射径的功率比）。当 $K \to 0$ 时，退化为瑞利衰落；当 $K \to \infty$ 时，退化为 AWGN 信道。
2. **复合信道系数**：结合路径损耗指数 $\nu=2.2$ 和随机距离 $d$，等效复合信道增益为：$h_{\text{eq}} = \frac{h_{\text{rician}}}{d^{\nu/2}}$

### 二、 仿真链路实现部分（第 160 - 163 行）

```matlab
160: n_ricrand = sqrt(N0/2) * (randn(N_samples, 1) + 1j*randn(N_samples, 1));
161: r_ricrand = h_rician_rand * sqrt(Eb) .* bits + n_ricrand;
162: r_eq_ricrand = real(conj(h_rician_rand) .* r_ricrand);
163: BER_RicRand_sim(idx) = mean(r_eq_ricrand .* bits < 0);

```

- **信号传输（Line 160-161）**：生成复高斯白噪声 `n_ricrand`（总方差为 $N_0$）。发送信号经过复合信道 `h_rician_rand` 并叠加噪声，得到接收信号 $r = h_{\text{eq}}\sqrt{E_b}s + n$。
- **相干解调与相位校正（Line 162）**：由于莱斯信道系数是复数（包含相位旋转），接收端必须乘以信道系数的复共轭 $\text{conj}(h_{\text{eq}})$ 进行**同相相干解调**，将信号能量拉回实轴，并取实部 `real` 滤除正交轴噪声。
- **错误计数（Line 163）**：通过判决矩阵乘积是否小于 0 (`.* bits < 0`) 判断是否误码，并求均值得到该 SNR 点的仿真 BER。

### 三、 理论解析与空间积分部分（第 164 - 165 行）
这是整个脚本中数学逻辑和 MATLAB 编程技巧最密集的两行代码：

```matlab
164: ber_cond_ric = @(d) arrayfun(@(dd) compute_rician_ber(gb/max(dd,eps)^nu, K), d);
165: BER_RicRand_ana(idx) = (2/R^2) * integral(@(d) d .* ber_cond_ric(d), 0, R);

```

#### 1. 条件误码率的定义与 `arrayfun` 的编程妙用（Line 164）
根据全概率公式，我们首先需要知道**当用户固定在距离 d 处时的条件误码率 Pb​(d)**。

- 当距离 $d$ 固定，大尺度路径损耗确定，用户的有效平均信噪比变为 $\gamma_b(d) = \frac{\gamma_b}{d^\nu}$（即代码中的 `gb/dd^nu`）。
- 我们需要调用前文提到的辅助函数 `compute_rician_ber` 来计算这个 SNR 下的莱斯理论误码率。
**★ 为什么这里要嵌套一个 arrayfun？**
MATLAB 的数值积分器 `integral` 在运行时，为了提高效率，会向被积函数传入一个**包含多个 d 值的向量**。然而，我们的辅助函数 `compute_rician_ber` 内部本身也使用了 `integral`（对 $\theta$ 积分），它**只能接收标量（单个数值）的信噪比输入**。
如果直接把向量 `d` 传给 `compute_rician_ber`，代码会因为维度不匹配而崩溃。

- `arrayfun(@(dd) ..., d)` 的作用是**将向量 d 拆解开来**，将其中的每一个元素（单个距离 `dd`）逐个带入后面的莱斯 BER 函数进行计算，最后再拼回成一个输出向量。这在 MATLAB 中属于高级的**向量化兼容处理**。
- `max(dd, eps)` 是一个**边界防御机制**：当积分上限接近原点 $d \to 0$ 时，避免分母为 0 导致浮点数溢出（`NaN` 或 `Inf`）。

#### 2. 全小区几何概率空间积分（Line 165）
用户在半径为 $R=3$ 的圆盘内均匀分布，其距离的概率密度函数（PDF）为 $f_D(d) = \frac{2d}{R^2}$。
利用全概率公式，对整个小区的空间进行积分以求得平均误码率：

$$
P_{b, \text{total}} = \int_{0}^{R} P_b(d) \cdot f_D(d) \, dd = \int_{0}^{R} P_b(d) \cdot \frac{2d}{R^2} \, dd = \frac{2}{R^2} \int_{0}^{R} d \cdot P_b(d) \, dd
$$
对应代码：`(2/R^2) * integral(@(d) d .* ber_cond_ric(d), 0, R)`。
通过这行自适应正交积分，MATLAB 完美算出了包含双重随机性（位置空间随机 + 莱斯小尺度衰落随机）下的最终理论误码率。

### 四、 深度物理洞察（UNSW 研究生报告加分要点）
在撰写 ELEC9123 实验报告时，你可以将**随机瑞利**与**随机莱斯**的结果进行深度对比分析：

1. **直射径（LOS）的红利**：在相同的平均信噪比下，**随机莱斯模型的 BER 表现明显优于随机瑞利模型**。这是因为莱斯信道中的直射径分量（$K=5$）提供了稳定的信号能量，极大地减小了信道发生“深衰落（Deep Fade）”的概率。
2. **近远效应的不可抗拒性**：尽管小尺度表现变好了，但如果你观察高 SNR 区域（12~15 dB）的曲线斜率，你会发现随机莱斯的下降速度依然非常缓慢。这说明**大尺度路径损耗（近远效应）仍然是决定蜂窝网络小区平均性能的绝对主导因素**。聚集在小区边缘（$d \to R$）的大量边缘用户，其有效 SNR 被路径损耗严重剥削，即使有直射径，也无法完全弥补大尺度衰减带来的性能红利衰退。

## Turn 11

### 🤖 Assistant

在 MATLAB 中，`dd` 并不是什么特殊的系统关键字，它只是一个**自定义的临时循环变量（标量）**。起名叫 `dd` 纯粹是为了和外层的距离向量 `d` 做出区分，寓意 “`d` 向量中的某一个具体元素（Element of d）”。
为了让你彻底明白为什么要这样写，我们需要拆解 MATLAB 的向量化（Vectorization）机制和 `arrayfun` 的底层逻辑。

### 1. 核心矛盾：向量（d）与 标量函数（compute_rician_ber）

- **d 是一个向量（Vector）：** MATLAB 的数值积分器 `integral` 在执行自适应积分时，为了追求运算速度，**绝对不会**一个点一个点地带入距离。它会一次性生成一组距离采样点，组合成一个**向量 d**（例如 `[0.1, 0.5, 1.2, 2.0]`）整体传给你的匿名函数。
- **compute_rician_ber 只能接收标量（Scalar）：**
看一下你的辅助函数内部，它自己也调用了一个 `integral(fun, 0, pi/2)` 对 $\theta$ 进行积分。MATLAB 的积分器内部是无法直接嵌套另一个向量化积分器的。这意味着，`compute_rician_ber` 在同一时间**只能计算某一个固定信噪比下的误码率**，它的输入必须是一个单独的数（标量）。
如果直接写 `compute_rician_ber(gb/d^nu, K)`，那就是把一堆距离（向量 `d`）直接怼进了只能吃单个数字的函数里，MATLAB 瞬间就会报 **"维度不匹配 (Matrix dimensions must agree)"** 的灾难性错误。

### 2. `arrayfun` 的解耦角色与 `dd` 的诞生
为了解决这个矛盾，你用到了 `arrayfun`。它的本质是**在底层用 C 语言优化过的一个隐式 for 循环**。

```matlab
ber_cond_ric = @(d) arrayfun(@(dd) compute_rician_ber(gb/max(dd,eps)^nu, K), d);

```
这行代码的内部执行逻辑等价于下面这段传统的 `for` 循环：

```matlab
% 假设 integral 传进来的向量 d 长度为 N
N = length(d);
result = zeros(1, N);

for i = 1:N
    dd = d(i);  % <--- 这就是 dd！它每次只从向量 d 里取出一个具体的标量值
    
    % 现在输入是标量 dd 了，可以安全地调用辅助函数
    result(i) = compute_rician_ber(gb / max(dd, eps)^nu, K); 
end

```
在 `arrayfun(@(dd) ..., d)` 语法中：

1. 最后的 `, d)` 告诉 MATLAB：“去把向量 `d` 给我拆了。”
2. 前面的 `@(dd)` 定义了一个内部的匿名函数，**dd 就是这个隐式循环的迭代变量**。每一次循环，`arrayfun` 会从 `d` 中提取一个元素赋值给 `dd`，然后代入后面的公式进行计算。

### 3. 如果改成别的名字可以吗？
**完全可以。**`dd` 只是作者个人的命名习惯。为了验证这一点，你可以把这行代码改成：

```matlab
ber_cond_ric = @(d) arrayfun(@(x) compute_rician_ber(gb/max(x,eps)^nu, K), d);

```
这里的 `x` 代替了 `dd`，代码的运行结果和效率完全一模一样。

### 总结

- **d**：是外层 `integral` 传进来的**整个距离向量**。
- **dd**（或 `x`）：是 `arrayfun` 内部每次抽出来的**单个距离标量**，专门用来伺候只能吃标量的 `compute_rician_ber` 函数。
这种在匿名函数内部嵌套 `arrayfun` 并在内部匿名函数中使用复合变量名（如 `dd`）的写法，是 MATLAB 高级矩阵编程中处理“标量函数无法处理向量输入”这一经典冲突的标准教科书式解法。

## Turn 12

### 🤖 Assistant

在 MATLAB 中，`eps` 是一个内置的常数，称为**机器 epsilon（Machine Epsilon）**，代表**浮点数的相对精度（Floating-point relative accuracy）**。
对于标准的研究生级数值计算与通信系统仿真，深入理解 `eps` 的数学本质和工程用途至关重要。以下为你详细拆解它的含义以及它在你的无线通信脚本中的核心作用：

### 一、 `eps` 的数学与计算机本质
在计算机中，实数是以二进制浮点数（IEEE 754 标准）的形式存储的。因为存储位数有限（MATLAB 默认的双精度 `double` 占用 64 位），计算机无法表示连续的实数轴，只能表示一系列离散的数字点。

- **定义**：`eps` 定义为 **1.0 与下一个比 1.0 大的浮点数之间的距离**。
- **具体数值**：在双精度格式下，`eps` 的精确值是 $2^{-52}$，大约等于：$\text{eps} \approx 2.2204 \times 10^{-16}$
- **特性**：如果你在 MATLAB 命令行输入 `1.0 + eps`，计算机会承认它大于 1；但如果你输入 `1.0 + (eps/2)`，由于超出了浮点数的存储精度，计算机在做四舍五入后，结果会直接退化为 `1.0`。

### 二、 在无线通信仿真代码中的核心作用：数值防御
在你的脚本第 164 行中，`eps` 出现在这句代码里：

```matlab
ber_cond_ric = @(d) arrayfun(@(dd) compute_rician_ber(gb/max(dd,eps)^nu, K), d);

```
这里的 `max(dd, eps)` 是一道非常经典的**数值边界防御机制（Boundary Defense Mechanism）**，其目的是**防止分母为 0 导致的程序崩溃（奇点保护）**。

#### 1. 致命的“分母归零”冲突
在计算单小区全概率空间积分时，MATLAB 的数值积分器 `integral` 会对距离 $d$ 在区间 $[0, R]$ 上进行自适应采样积分。

- 当积分器采样到圆心位置，即 **d=0** 时，大尺度路径损耗项 $d^\nu$ 就会变成 $0^\nu = 0$。
- 此时，条件信噪比公式为 $\gamma_b(d) = \frac{\gamma_b}{d^\nu}$。如果不做任何处理，代码就会执行 `gb / 0`。
- 在数学上，`gb / 0` 会直接触发计算机输出 **Inf（无穷大）**。
- 随后，这个 `Inf` 被传入后面的莱斯误码率或中断概率函数中，在涉及指数运算（如 $\exp(-\text{Inf})$）或进一步的矩阵运算时，极易诱发底层浮点数产生 **NaN（Not a Number，非数）**。一旦积分器收到 `NaN`，整个自适应步长算法就会瞬间报错并中断执行。

#### 2. `max(dd, eps)` 的解题思路
通过引入 `max(dd, eps)`，你为系统设定了一个极其微小的正下限：

- 当用户距离基站 $d > 0$ 时，`max(dd, eps)` 依然等于 $d$，不影响任何大尺度路径损耗的物理真实性。
- 当积分器无限逼近圆心，甚至直接取 $d = 0$ 时，`max(0, eps)` 会强行将距离替换为 $\text{eps} \approx 2.22 \times 10^{-16}$。
- 此时，分母变成了 $\text{eps}^\nu$，计算出的条件信噪比虽然极其巨大（如 $10^{30}$ 数量级），但在计算机内存中它依然是一个**合法且有限的浮点数**。莱斯误码率函数能够平滑地输出一个趋近于 $0$ 的合法误码率，从而保护了数值积分的连续性，避免了崩溃。

### 三、 工程总结
在科学计算与通信工程中，`eps` 常被用作一种**正则化扰动项（Regularization Factor）或数值保护伞**。当你写出的数学公式中存在潜在的积分原点奇点（Singularity，如 $\frac{1}{x}$，$\ln(x)$，$x^{-\nu}$ 等），且积分下限包含 $0$ 时，利用 `eps` 进行截断是国际学术界和工程界最标准的规范写法。

## Turn 13

### 🤖 Assistant

这段代码实现了**BPSK系统在纯瑞利衰落信道（当距离固定为 d=1 米，即忽略大尺度路径损耗）下的瞬时信噪比计算，以及仿真和理论“中断概率（Outage Probability）”的对比验证**。
在移动通信中，中断（Outage）是指当信道质量变得太差，导致瞬时信道容量无法满足系统设定的最低传输速率要求时的状态。
以下为你从**物理量映射**、**仿真逻辑**和**理论推导**三个维度进行详细且清晰的拆解：

### 一、 瞬时信噪比计算（第一行）

```matlab
inst_snr_ray = gb * abs(h_rayleigh).^2;

```

- **通信物理量**：这行代码计算的是**瞬时信噪比 γinst​**（Instantaneous SNR）。
- **数学公式**：$\gamma_{\text{inst}} = \gamma_b \cdot \vert{}h\vert{}^2$其中：
  - `gb`（$\gamma_b$）是系统设定的**平均信噪比**（线性值）。
  - `h_rayleigh`（$h$）是小尺度瑞利衰落系数，其幅度 $\vert{}h\vert{}$ 服从瑞利分布。
  - `abs(h_rayleigh).^2`（$\vert{}h\vert{}^2$）是**信道功率增益**（Channel Power Gain）。由于瑞利衰落的功率归一化为 1（$E[\vert{}h\vert{}^2]=1$），因此 $\vert{}h\vert{}^2$ 服从均值为 1 的**指数分布**（Exponential Distribution）。

### 二、 蒙特卡洛仿真中断概率（第二行）

```matlab
Pout_Ray_sim(idx) = mean(log2(1 + inst_snr_ray) < C);

```

- **通信原理（中断的定义）**：
根据香农公式（Shannon's Formula），瞬时信道的无差错最高传输速率（即瞬时信道容量 $C_{\text{inst}}$）为：$C_{\text{inst}} = \log_2(1 + \gamma_{\text{inst}}) \quad (\text{bps/Hz})$系统设定了一个最低的目标容量门限 `C`（例如代码开头的 $1.2\text{ bps/Hz}$）。如果瞬时容量跌破了这个门限，系统就会发生**中断**。
- **代码执行机制**：
  - `log2(1 + inst_snr_ray)`：计算这 $10^6$ 个随机信道状态下，每个样本点对应的瞬时信道容量。
  - `< C`：逐元素进行逻辑判断。如果瞬时容量小于 `C`（发生中断），该位置返回 `1 (true)`；否则返回 `0 (false)`。这会生成一个由 0 和 1 组成的逻辑向量。
  - `mean(...)`：对该逻辑向量求均值。本质上就是统计 $\frac{\text{中断发生次数}}{\text{总仿真样本数}}$，从而利用蒙特卡洛频率法估计出当前的**仿真中断概率**。

### 三、 理论解析式计算（第三、四行）

```matlab
% 理论 outage：|h_s|^2 ~ Exp(1)，故 P_out = Pr(|h_s|^2 < T) = 1 - exp(-T)
Pout_Ray_ana(idx) = 1 - exp(-T);

```
这是研究生阶段需要熟练推导的经典平坦衰落信道指标，其严密的数学推导过程如下：

1. **门限转化**：
根据定义，中断概率 $P_{\text{out}}$ 是瞬时容量小于门限 $C$ 的概率：$P_{\text{out}} = P(C_{\text{inst}} < C) = P\left(\log_2(1 + \gamma_b \vert{}h\vert{}^2) < C\right)$对不等式两边同时做指数变换和化简：$1 + \gamma_b \vert{}h\vert{}^2 < 2^C \implies \vert{}h\vert{}^2 < \frac{2^C - 1}{\gamma_b}$在你的循环开头，定义了 $T = \frac{\gamma_{\text{th}}}{\gamma_b} = \frac{2^C - 1}{\gamma_b}$（即代码中的 `T = gamma_th / gb;`）。
因此，中断的条件完美等价于：**信道功率增益 ∣h∣2 小于折算门限 T**。$P_{\text{out}} = P(\vert{}h\vert{}^2 < T)$
2. **利用 CDF 求解**：
正如前文所述，对于纯瑞利衰落，信道功率增益 $X = \vert{}h\vert{}^2$ 服从参数为 $\lambda = 1$ 的标准指数分布（$\text{Exp}(1)$）。
指数分布的概率密度函数（PDF）为 $f_X(x) = e^{-x}$（当 $x \ge 0$ 时）。
那么，其累积分布函数（CDF）即为：$F_X(x) = P(X \le x) = \int_{0}^{x} e^{-t} dt = 1 - e^{-x}$将折算门限 $T$ 代入该 CDF 公式，直接得出理论中断概率的闭式解（Closed-form Expression）：$P_{\text{out}} = F_X(T) = 1 - e^{-T}$代码第 174 行的 `1 - exp(-T)` 正是这个最终数学解析式的直接落地。

### 四、 UNSW 研究生学术报告分析建议
在撰写 ELEC9123 的设计任务报告时，针对这几行代码对应的图 3（Outage 概率曲线），你可以写入以下深度分析：

1. **高信噪比下的中断表现（Outage Floor 现象）**：
当平均信噪比 $\gamma_b \to \infty$ 时，门限 $T = \frac{2^C - 1}{\gamma_b} \to 0$。此时由泰勒展开可得：$P_{\text{out}} = 1 - e^{-T} \approx T = \frac{2^C - 1}{\gamma_b}$这说明在对数坐标轴（`semilogy`）上，瑞利信道的中断概率曲线随着 SNR 的增大，会表现为一条**斜率为 -1 的直线**（即信噪比每增加 10 dB，中断概率下降一个数量级）。这能够非常直观地通过你的仿真曲线得到验证。
2. **与随机部署模型的对比引申**：
此处 $d = 1\text{ m}$ 作为基准线，其性能是最优的。当后续引入随机用户部署（大尺度衰落）后，由于边缘用户的有效信噪比被路径损耗剧烈剥削，折算门限 $T$ 会因距离而变大，从而导致小区的平均中断概率（Outage）大幅恶化。这部分纯瑞利理论值是你后续复合信道分析的重要基准（Baseline）。

## Turn 14

### 🤖 Assistant

这段代码实现了**BPSK系统在纯莱斯（Rician）衰落信道下（固定距离 d=1 米，忽略大尺度路径损耗）的瞬时信噪比计算，以及仿真和理论“中断概率（Outage Probability）”的对比验证**。
与纯瑞利信道不同，莱斯信道模拟的是存在强直射径（LOS, Line-of-Sight）的无线通信环境。这种直射径的存在，使得信道的统计特性从指数分布演变为非中心卡方分布。
以下为您从**瞬时信噪比物理模型**、**蒙特卡洛仿真逻辑**以及**基于 Marcum Q 函数的理论推导**三个维度进行深入剖析：

### 一、 瞬时信噪比计算（第一行）

```matlab
inst_snr_ric = gb * abs(h_rician).^2;

```

- **物理模型**：计算莱斯衰落环境下的**瞬时信噪比 γinst​**。
- **数学公式**：$\gamma_{\text{inst}} = \gamma_b \cdot \vert{}h_{\text{rician}}\vert{}^2$其中：
  - `gb`（$\gamma_b$）为系统的**平均信噪比**。
  - `h_rician`（$h$）是包含了直射径和小尺度散射径的复合复信道系数，其幅度 $\vert{}h\vert{}$ 服从莱斯分布。
  - `abs(h_rician).^2`（$\vert{}h\vert{}^2$）是信道功率增益。由于小尺度衰落的功率已归一化（$E[\vert{}h\vert{}^2]=1$），这里的 $\vert{}h\vert{}^2$ 服从自由度为 2 的**非中心卡方分布（Non-central Chi-square Distribution）**。

### 二、 蒙特卡洛仿真中断概率（第二行）

```matlab
Pout_Ric_sim(idx) = mean(log2(1 + inst_snr_ric) < C);

```

- **通信原理**：同样基于香农公式，瞬时信道容量 $C_{\text{inst}} = \log_2(1 + \gamma_{\text{inst}})$。如果该瞬时容量小于系统预设的目标传输速率门限 `C`（如 1.2 bps/Hz），则触发系统“中断”。
- **机制拆解**：
  1. `log2(1 + inst_snr_ric)` 算出一维数组中 $10^6$ 个独立信道状态分别对应的瞬时容量。
  2. `< C` 进行逐元素布尔判断，发生中断的位置标记为 `1`，未中断标记为 `0`。
  3. `mean(...)` 通过对 0/1 逻辑向量求均值，统计出中断事件发生的相对频率，即为蒙特卡洛估计出的**仿真中断概率**。

### 三、 理论解析式与 Marcum Q 函数（第三至五行）

```matlab
% 理论 outage：|h|^2 服从非中心卡方分布，可用 Marcum Q 函数
% P_out = 1 - Q_1(sqrt(2K), sqrt(2(K+1)T))
% MATLAB 中 marcumq(a, b) 即 Q_1(a, b)
Pout_Ric_ana(idx) = 1 - marcumq(sqrt(2*K), sqrt(2*(K+1)*T));

```
这里的数学推导在无线通信理论中具有极高的理论权重，其推导核心如下：

#### 1. 门限折算
与瑞利信道一致，中断的本质是瞬时容量跌破门限，可完美等价化简为**信道功率增益 ∣h∣2 小于折算门限 T**（其中 $T = \frac{2^C - 1}{\gamma_b}$）：

$$
P_{\text{out}} = P(C_{\text{inst}} < C) = P(\vert{}h\vert{}^2 < T)
$$

#### 2. 非中心卡方分布与 Marcum Q 函数的引入
由于信道存在直射径，信道功率增益 $X = \vert{}h\vert{}^2$ 满足非中心卡方分布。其累积分布函数（CDF）无法表示为简单的指数函数或初等函数，必须借助特殊的数学工具——**一阶广义 Marcum Q 函数（First-order Marcum Q-function, Q1​(a,b)）**。
Marcum Q 函数的标准数学定义为：

$$
Q_1(a, b) = \int_{b}^{\infty} x \exp\left( -\frac{x^2 + a^2}{2} \right) I_0(ax) \, dx
$$
在通信中，它物理上代表了**非中心卡方变量超过某个门限的尾部概率（即未中断的概率，或成功传输的概率 Psucc​）**：

$$
P_{\text{succ}} = P(\vert{}h\vert{}^2 \ge T) = Q_1\left( \sqrt{2K}, \sqrt{2(K+1)T} \right)
$$
其中参数的物理含义为：

- 参数 $a = \sqrt{2K}$：与直射径（LOS）分量的确定性功率直接相关。
- 参数 $b = \sqrt{2(K+1)T}$：将折算门限 $T$ 映射到标准化随机变量的空间中。

#### 3. 中断概率的最终求解
由于中断概率是成功传输概率的对立事件（补集）：

$$
P_{\text{out}} = 1 - P_{\text{succ}} = 1 - Q_1\left( \sqrt{2K}, \sqrt{2(K+1)T} \right)
$$
MATLAB 内置了高效计算该积分的特殊函数 `marcumq(a, b)`。代码第 185 行的 `1 - marcumq(sqrt(2*K), sqrt(2*(K+1)*T))` 正是这一高阶通信理论公式的完全落地。

### 四、 UNSW 研究生学术报告分析引申（加分项）
在您的 ELEC9123 实验报告中，对于莱斯信道中断概率与瑞利信道中断概率的对比，可以提出以下深度的学术见解：

1. **直射径带来的“陡峭度”红利**：
如果您观察生成的图 3，会发现**纯莱斯信道（d=1 m）的中断概率曲线随 SNR 增加而下降的速度（斜率）明显快于纯瑞利信道**。在通信原理中，这意味着虽然两者的多样性阶数（Diversity Order）在渐进意义上都为 1，但莱斯信道拥有极大的**编码增益（Coding Gain）**。强直射径（$K=5$）像一根“定海神针”，将瞬时信道增益牢牢锁定在均值附近，极大地抑制了信道功率向原点（$0$）塌陷引发的深衰落。
2. **莱斯因子 K 的物理边界演变**：
可以在报告的 Discussion 中补充：当 $K=0$ 时，$\sqrt{2K}=0$，Marcum Q 函数退化为标准指数分布积分，公式将完美还原为瑞利衰落的 $1 - e^{-T}$；而当 $K \to \infty$ 时，信道没有任何随机多径，系统将退化为纯 AWGN 信道，中断概率曲线将表现为类似阶跃函数的陡峭陡降。这体现了 Marcum Q 解析式在边界条件下的数学完备性。

## Turn 15

### 🤖 Assistant

这段代码实现了**BPSK系统在“随机用户部署 + 瑞利衰落 + AWGN”复合信道模型下的瞬时信噪比计算，以及仿真和理论“中断概率（Outage Probability）”的对比验证**。
这是整个仿真脚本中在数学推导上非常漂亮的一段。它同时考虑了**大尺度的几何路径损耗（Path Loss）与小尺度的多径瑞利衰落（Rayleigh Fading）**。
以下为你从**物理量映射**、**蒙特卡洛仿真逻辑**以及**解析闭式解的严格数学推导**三个维度进行深入剖析：

### 一、 瞬时信噪比计算（第一行）

```matlab
inst_snr_rayrand = gb * abs(h_rayleigh).^2 ./ d_user.^nu;

```

- **通信物理模型**：这行代码计算的是包含双重随机性的**瞬时信噪比 γinst​**。
- **数学公式**：$\gamma_{\text{inst}}(h, d) = \gamma_b \cdot \frac{\vert{}h\vert{}^2}{d^\nu}$其中：
  - `gb` ($\gamma_b$) 是基站处的**平均信噪比**。
  - `abs(h_rayleigh).^2` ($\vert{}h\vert{}^2$) 是小尺度瑞利衰落功率增益，服从均值为 1 的**指数分布**。
  - `d_user.^nu` ($d^\nu$) 是大尺度路径损耗，其中 `d_user` ($d$) 是用户在圆形小区内随机分布的距离，`nu` ($\nu$) 是路径损耗指数。

### 二、 蒙特卡洛仿真中断概率（第二行）

```matlab
Pout_RayRand_sim(idx) = mean(log2(1 + inst_snr_rayrand) < C);

```

- **仿真机制**：瞬时信道容量 $C_{\text{inst}} = \log_2(1 + \gamma_{\text{inst}})$。如果这 $10^6$ 个随机样本（每个样本点都包含一个随机位置和独立的信道衰落）计算出的瞬时容量小于最低传输门限 `C`，系统即判定发生中断。通过 `mean(...)` 统计这 $10^6$ 个样本中发生中断的频率，即为**仿真中断概率**。

### 三、 理论闭式解与不完全伽马函数（第三至五行）

```matlab
% 理论 outage：利用圆盘均匀部署的 PDF 和指数衰落积分得到闭式
% P_out = 1 - (2/(ν*R^2)) * T^(-2/ν) * γ(2/ν, T*R^ν)
% 其中 γ(s, x) 为下不完全伽马函数，MATLAB 用 gammainc(x, s, 'lower')*gamma(s)
Pout_RayRand_ana(idx) = 1 - (2/(nu*R^2)) * T^(-2/nu) * gamma(2/nu) * gammainc(T*R^nu, 2/nu, 'lower');

```
这里的理论解析式并非数值积分，而是一个**纯代数闭式解（Closed-form Expression）**。其严格的数学推导步骤如下：

#### 步骤 1：全概率积分公式构建
与前面一致，中断条件可以转化为瞬时信噪比小于门限 $\gamma_{\text{th}} = 2^C - 1$：

$$
P_{\text{out}} = P\left( \gamma_b \frac{\vert{}h\vert{}^2}{d^\nu} < \gamma_{\text{th}} \right) = P\left( \vert{}h\vert{}^2 < \frac{\gamma_{\text{th}}}{\gamma_b} d^\nu \right) = P(\vert{}h\vert{}^2 < T \cdot d^\nu)
$$
其中 $T = \gamma_{\text{th}} / \gamma_b$。由于距离 $d$ 是随机变量，根据全概率公式，我们需要对用户在半径为 $R$ 的圆盘空间内的几何分布 PDF $f_D(d) = \frac{2d}{R^2}$ 进行积分：

$$
P_{\text{out}} = \int_{0}^{R} P\left( \vert{}h\vert{}^2 < T \cdot d^\nu \mid d \right) \cdot f_D(d) \, dd
$$

#### 步骤 2：对小尺度衰落求条件概率
由于固定距离 $d$ 后，小尺度功率增益 $\vert{}h\vert{}^2$ 服从标准指数分布（CDF 为 $1 - e^{-x}$），因此条件中断概率为：

$$
P\left( \vert{}h\vert{}^2 < T \cdot d^\nu \mid d \right) = 1 - \exp(-T \cdot d^\nu)
$$

#### 步骤 3：带入大尺度空间积分并进行元变换
将条件概率带回全空间积分式中：

$$
P_{\text{out}} = \int_{0}^{R} \left( 1 - \exp(-T \cdot d^\nu) \right) \cdot \frac{2d}{R^2} \, dd = 1 - \frac{2}{R^2} \int_{0}^{R} d \cdot \exp(-T \cdot d^\nu) \, dd
$$
这是一个非常经典的积分。为了求解右边的积分项 $I = \int_{0}^{R} d \cdot \exp(-T \cdot d^\nu) \, dd$，我们引入**换元法**：
令 $x = T \cdot d^\nu \implies d = \left(\frac{x}{T}\right)^{1/\nu}$，对两边求导可得 $dd = \frac{1}{\nu T} \left(\frac{x}{T}\right)^{\frac{1}{\nu}-1} dx$。
代入积分项，并将积分上下限由 $d \in [0, R]$ 转换为 $x \in [0, T \cdot R^\nu]$：

$$
I = \int_{0}^{T R^\nu} \left(\frac{x}{T}\right)^{1/\nu} \cdot e^{-x} \cdot \frac{1}{\nu T} \left(\frac{x}{T}\right)^{\frac{1}{\nu}-1} \, dx = \frac{1}{\nu T^{2/\nu}} \int_{0}^{T R^\nu} x^{\frac{2}{\nu}-1} e^{-x} \, dx
$$

#### 步骤 4：引入下不完全伽马函数（Lower Incomplete Gamma Function）
数学上，**下不完全伽马函数**定义为：

$$
\gamma(s, z) = \int_{0}^{z} x^{s-1} e^{-x} \, dx
$$
对比我们的积分项，正好满足 $s = \frac{2}{\nu}$，$z = T \cdot R^\nu$。因此：

$$
I = \frac{1}{\nu T^{2/\nu}} \cdot \gamma\left(\frac{2}{\nu}, T R^\nu\right)
$$
将积分结果 $I$ 带回 $P_{\text{out}}$ 的总公式中，并利用恒等式 $\frac{2}{R^2} \cdot \frac{1}{\nu} = \frac{2}{\nu R^2}$ 化简，最终得到：

$$
P_{\text{out}} = 1 - \frac{2}{\nu R^2} T^{-2/\nu} \cdot \gamma\left(\frac{2}{\nu}, T R^\nu\right)
$$

#### 步骤 5：MATLAB 函数映射
在 MATLAB 中，内置的 `gammainc(z, s, 'lower')` 计算的是正规化（Normalized）的下不完全伽马函数，即：

$$
\text{gammainc}(z, s, \text{'lower'}) = \frac{\gamma(s, z)}{\Gamma(s)}
$$
其中 $\Gamma(s)$ 是标准的伽马函数（MATLAB 中用 `gamma(s)` 表示）。
因此，真实的下不完全伽马函数在 MATLAB 中必须表示为：

$$
\gamma(s, z) = \text{gamma}(s) \cdot \text{gammainc}(z, s, \text{'lower'})
$$
把 $s = 2/\nu$ 和 $z = T \cdot R^\nu$ 代入，代码第 201 行的代数式与上述数学推导完全一致。

### 四、 UNSW 研究生学术报告分析要点（高分加分项）
在撰写 ELEC9123 实验报告时，针对该复合信道的中断概率曲线，可以提出以下高阶学术洞察：

1. **路径损耗指数 ν 的敏感性分析**：
闭式解中包含 $T^{-2/\nu}$ 项。当平均信噪比很高时（$\gamma_b \to \infty$），折算门限 $T = \gamma_{\text{th}}/\gamma_b \to 0$。此时根据不完全伽马函数的级数展开性质，整个小区的平均中断概率渐进表现为：$P_{\text{out}} \propto \left(\frac{1}{\gamma_b}\right)^{2/\nu}$这意味着在对数坐标轴上，曲线的斜率变为了 **−2/ν**。在纯瑞利信道中斜率是 $-1$，而当 $\nu = 2.2$ 时，斜率变为了 $-2/2.2 \approx -0.909$。**曲线变缓了！** 这在物理上证明了：**大尺度路径损耗的引入，直接损害了小区的平均中断多样性增益（Diversity Gain）**。
2. **边缘用户对中断性能的统治力**：
不完全伽马函数的积分上限是 $T \cdot R^\nu$。小区半径 $R=3\text{ m}$ 经过 $\nu=2.2$ 次方放大后变为 $3^{2.2} \approx 11.26$。这表明**小区边缘的用户由于距离远、路径损耗大，贡献了绝大多数的中断事件**。这种利用特殊函数解析式直接预测复合随机系统性能的方法，比单纯的蒙特卡洛仿真更能体现出蜂窝系统规划中近远效应的物理本质。

## Turn 16

### 🤖 Assistant

这段代码实现了**BPSK系统在“随机用户部署 + 莱斯（Rician）衰落 + AWGN”复合信道模型下的瞬时信噪比计算，以及仿真和理论“中断概率（Outage Probability）”的对比验证**。
这是整个脚本中在数学架构上最健壮、最贴近实际微蜂窝基站规划的指标计算部分。由于小尺度衰落引入了包含视距（LOS）分量的莱斯信道，使得空间积分无法像瑞利信道那样化简为纯代数的闭式解，因此这里采用了极其严谨的半解析数值积分（Semi-analytical Numerical Integration）方法。
以下为您从**物理模型构建**、**蒙特卡洛仿真**以及**全概率空间数值积分**三个维度进行深度拆解：

### 一、 瞬时信噪比与复合信道（第一行）

```matlab
inst_snr_ricrand = gb * abs(h_rician).^2 ./ d_user.^nu;

```

- **物理模型**：计算同时兼顾大尺度路径损耗和小尺度莱斯衰落的**复合瞬时信噪比 γinst​**。
- **数学公式**：$\gamma_{\text{inst}}(h, d) = \gamma_b \cdot \frac{\vert{}h_{\text{rician}}\vert{}^2}{d^\nu}$其中：
  - `gb` ($\gamma_b$) 为基站发射侧定义的**平均信噪比**。
  - `abs(h_rician).^2` ($\vert{}h\vert{}^2$) 为莱斯衰落功率增益，由于存在直射径（$K=5$），它服从**非中心卡方分布**。
  - `d_user.^nu` ($d^\nu$) 为大尺度路径损耗，其中距离 $d$ 在圆形小区内均匀分布。

### 二、 蒙特卡洛仿真中断概率（第二行）

```matlab
Pout_RicRand_sim(idx) = mean(log2(1 + inst_snr_ricrand) < C);

```

- **仿真机制**：瞬时容量 $C_{\text{inst}} = \log_2(1 + \gamma_{\text{inst}})$。通过矩阵化运算，一次性判断 $10^6$ 个包含随机位置和独立莱斯衰落的样本是否跌破目标容量门限 `C`。利用 `mean(...)` 统计这 $10^6$ 个点中发生中断（返回 1）的比例，得到**仿真中断概率**。

### 三、 理论全空间数值积分与匿名函数嵌套（第三至五行）

```matlab
% 理论 outage：对每个距离 d 的条件 outage 用 Marcum Q 函数，再对 d 积分
pout_cond_ric = @(d) 1 - marcumq(sqrt(2*K), sqrt(2*(K+1)*T*(d.^nu + eps)));
Pout_RicRand_ana(idx) = (2/R^2) * integral(@(d) d .* pout_cond_ric(d), 0, R);

```
这两行代码展示了如何运用全概率公式（Total Probability Theorem）在几何空间上对高阶无线通信指标进行求解。

#### 1. 构建条件中断概率函数 `pout_cond_ric`（Line 206）
首先，假设用户的空间位置被固定在距离基站 $d$ 处。此时，大尺度路径损耗固定为 $d^\nu$。
由于小尺度衰落是纯莱斯衰落，根据中断的定义，瞬时信噪比跌破门限 $\gamma_{\text{th}} = 2^C - 1$ 的条件为：

$$
P_{\text{out}}(d) = P\left( \gamma_b \frac{\vert{}h\vert{}^2}{d^\nu} < \gamma_{\text{th}} \right) = P\left( \vert{}h\vert{}^2 < \frac{\gamma_{\text{th}}}{\gamma_b} d^\nu \right) = P(\vert{}h\vert{}^2 < T \cdot d^\nu)
$$
其中 $T = \gamma_{\text{th}} / \gamma_b$（即循环外的 `T = gamma_th / gb`）。
我们知道，纯莱斯衰落在功率门限为 $T_{\text{eff}}$ 时的理论中断概率由一阶广义 Marcum Q 函数给出：$1 - Q_1(\sqrt{2K}, \sqrt{2(K+1)T_{\text{eff}}})$。
在这里，**等效功率门限变成了 Teff​=T⋅dν**。将此项代入，即可得到特定距离 $d$ 下的条件中断概率解析式：

$$
P_{\text{out}}(d) = 1 - Q_1\left( \sqrt{2K}, \sqrt{2(K+1)T \cdot d^\nu} \right)
$$
代码第 206 行的匿名函数正是这个公式的完美对映。`d.^nu + eps` 同样是一道数值防御，防止在原点 $d=0$ 时可能引发的开方或矩阵维度计算异常。

#### 2. 全小区几何空间加权积分（Line 207）
由于用户在半径为 $R$ 的圆盘内均匀分布，距离 $d$ 的概率密度函数（PDF）为 $f_D(d) = \frac{2d}{R^2}$。
为了得到全小区的平均中断概率，必须对整个小区的空间几何分布进行全概率积分：

$$
P_{\text{out}, \text{total}} = \int_{0}^{R} P_{\text{out}}(d) \cdot f_D(d) \, dd = \int_{0}^{R} P_{\text{out}}(d) \cdot \frac{2d}{R^2} \, dd = \frac{2}{R^2} \int_{0}^{R} d \cdot P_{\text{out}}(d) \, dd
$$
对应代码：

```matlab
Pout_RicRand_ana(idx) = (2/R^2) * integral(@(d) d .* pout_cond_ric(d), 0, R);

```
**★ 为什么这里不用闭式解，而用 integral 数值积分？**
在之前的“随机用户+瑞利衰落”中，条件中断概率是纯指数函数 $1 - \exp(-T d^\nu)$，形式简单，可以通过换元法积出不完全伽马函数的闭式解。而在此处，莱斯的条件中断概率包含了 **Marcum Q 函数（本质上是一个无穷积分）**。将 Marcum Q 函数乘以 $d$ 再对 $d$ 进行空间积分，在数学上**无法推导出由初等函数或常见特殊函数表示的代数闭式解**。
因此，脚本在这里采用了学术界最标准的做法：利用 MATLAB 的自适应高斯-克朗罗德复合积分器 `integral`，在空间尺度 $[0, R]$ 上进行高精度的数值求解。

### 四 bias、UNSW 研究生学术报告分析引申（高分加分项）
在您的 ELEC9123 实验报告中，对于该复合模型的中断性能表现，可以展开如下深度物理洞察：

1. **“直射径红利”与“近远效应”的博弈（Trade-off）**：
将**随机莱斯中断曲线**与**随机瑞利中断曲线**进行横向对比。您会发现，在相同的平均信噪比下，随机莱斯的中断概率整体更低。这说明**直射径（K=5）的存在不仅能降低误码率，还能有效降低系统发生中断的概率**。然而，随着平均 SNR 增大，随机莱斯曲线的斜率变缓趋势与随机瑞利高度一致。这从数学上证明了：**大尺度路径损耗引起的小区边缘“功率剥削”，对中断多样性阶数（Diversity Order）的负面物理影响是全局性的，小尺度上的直射径只能改善系统规划的基准线（编码增益），无法扭转近远效应带来的斜率惩罚**。
2. **半解析验证法的学术优越性**：
在报告中可以指出：本任务通过蒙特卡洛链路级仿真（Simulation）与基于 Marcum Q 函数的空间自适应数值积分（Theoretical Semianalytical）两条完全独立的路径得到了高度拟合的结果。这不仅验证了仿真器中物理信道生成的正确性，也证明了利用现代通信理论（Marcum Q 函数 + 全概率空间积分）预测复杂复合随机网络业务中断性能的精确度与可行性。

## Turn 17

### 🤖 Assistant

这份 MATLAB 脚本编写得非常出色，逻辑严密、注释详尽，并且理论推导与代码实现的对应关系非常精准。这完全符合高阶无线通信课程（如 ELEC9123）的专业要求。
你不仅完成了基础的仿真，还考虑到了代码运行效率和结果的严谨性。以下是对你代码中几个核心亮点的详细解析，以及从实际工程和科研角度可以进一步探讨的微调建议。

## 代码核心亮点解析

- **随机变量预生成（Pre-generation of RVs）：**
你在进入 SNR 循环前，一次性生成了 `bits`, `h_rayleigh`, `h_rician`, `d_user`。这是一个非常优秀的编程习惯。它不仅大幅减少了循环内的计算开销，更重要的是保证了在不同 SNR 下，系统经历的是**同一套信道衰落和用户位置**。这种控制变量法使得各 SNR 点之间的曲线更加平滑，对比更加科学。
- **拉普拉斯分布的逆变换采样（Inverse Transform Sampling）：**
你准确地推导了 Laplace 分布的 CDF 及其逆函数 $x = -b \cdot \text{sgn}(u-0.5) \cdot \ln(1 - 2\vert{}u-0.5\vert{})$，并专门编写了任务 2 的验证模块（直方图与理论 PDF/CDF 对比）。这种“先验证信源，再做系统仿真”的思路是极其严谨的。
- **莱斯衰落 BER 的 MGF 积分法：**
莱斯衰落（Rician fading）下的 BER 闭式解非常复杂。你选择利用矩生成函数（MGF）结合基于 $\theta$ 的积分表示法：$P_b = \frac{1}{\pi} \int_{0}^{\pi/2} M_{\gamma}\left(-\frac{1}{\sin^2\theta}\right) d\theta$配合 MATLAB 的 `integral` 函数求解，既保证了理论值的极高精度，又避免了复杂的贝塞尔函数带来的数值不稳定性。
- **随机用户部署模型：**
通过 $d = R\sqrt{u}$ 生成圆盘内均匀分布的用户，并将大尺度路径损耗（Path loss $d^{-\nu}$）与小尺度衰落无缝结合。在计算理论 Outage 概率时，巧妙调用了不完全伽马函数 `gammainc`，数学功底十分扎实。

## 进阶探讨与优化建议
尽管代码已经达到了很高的水准，但在实际的蒙特卡洛仿真中，还有一个经典问题值得注意：**高 SNR 下的样本数瓶颈（Confidence Interval Issue）**。
你目前设定的是 $N = 10^6$。在蒙特卡洛仿真中，有一个经验法则：**要准确估计某个概率 P，至少需要发生 10 到 100 次错误事件**。

- 当 SNR 达到 12dB~15dB 时，AWGN 或 Laplacian 信道下的理论 BER 会骤降到 $10^{-8}$ 甚至更低。
- 在 $10^6$ 的样本量下，错误发生的期望次数不到 1 次。这意味着在这些高 SNR 点上，你的仿真 BER 极大概率会直接变为 $0$。
- **你的处理：** 你在绘图时使用了 `max(BER_AWGN_sim, 1e-7)`，这非常聪明地避免了对数坐标轴报错，保证了图像的完整性。
**如果你希望在报告中进一步体现深度，可以考虑以下改进思路（可选）：**

---

*Exported from [Voyager](https://github.com/Nagi-ovo/gemini-voyager)*  
*Generated on July 16, 2026 at 11:03 PM*