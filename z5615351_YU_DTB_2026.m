%% ELEC9123: Design Task B (Wireless Communications)
% Performance Analysis Verification via Monte Carlo Simulations
%
% Author: YU
% zID: z5615351
% Filename: z5615351_YU_DTB_2025.m
%
% This script implements Monte Carlo simulations for a BPSK system under
% six wireless channel models and verifies BER and outage probability
% against analytical expressions.

clc;
clear;
close all;

%% ========================== SYSTEM PARAMETERS ==========================
C = 1.2;                         % Outage threshold (bps/Hz)
SNR_dB = 0:1:15;                 % Average SNR per bit in dB
SNR_lin = 10.^(SNR_dB/10);       % Linear average SNR per bit (gamma_b)
R = 3;                           % Cell radius (m)
nu = 2.2;                        % Path loss exponent
hd = 1;                          % Magnitude of LOS component
K = 5;                           % Rician K-factor
N0 = 1;                          % Noise variance
N_samples = 1e6;                 % Monte Carlo samples (>= 1e6)
gamma_th = 2^C - 1;              % SNR threshold for outage

fprintf('ELEC9123 Design Task B - Monte Carlo Simulation\n');
fprintf('Author: YU (z5615351)\n');
fprintf('Samples per SNR point: %d\n\n', N_samples);

%% ============== TASK 2: VERIFY LAPLACIAN RANDOM VARIABLES ==============
% Generate Laplacian noise with variance N0 using inverse CDF method.
% For Laplace(0, b) with variance 2b^2 = N0, b = sqrt(N0/2).
% CDF: F(x) = 0.5*exp(sqrt(2/N0)*x) for x<=0,
%            = 1 - 0.5*exp(-sqrt(2/N0)*x) for x>0.

N_lap = 1e6;
b_lap = sqrt(N0/2);
u_lap = rand(N_lap, 1);
% Inverse CDF: x = -b_lap * sign(u-0.5) * log(1 - 2|u-0.5|)
n_lap = -b_lap * sign(u_lap - 0.5) .* log(1 - 2*abs(u_lap - 0.5));

% Analytical PDF and CDF
x_lap = linspace(-6*b_lap, 6*b_lap, 1000);
f_lap_ana = (1/(2*b_lap)) * exp(-abs(x_lap)/b_lap);
F_lap_ana = zeros(size(x_lap));
F_lap_ana(x_lap <= 0) = 0.5 * exp(x_lap(x_lap <= 0)/b_lap);
F_lap_ana(x_lap > 0) = 1 - 0.5 * exp(-x_lap(x_lap > 0)/b_lap);

% Empirical CDF
F_lap_sim = zeros(size(x_lap));
for i = 1:length(x_lap)
    F_lap_sim(i) = mean(n_lap <= x_lap(i));
end

% Plot verification
figure('Name', 'Laplacian Verification', 'Position', [100 100 900 400]);
subplot(1, 2, 1);
histogram(n_lap, 100, 'Normalization', 'pdf', 'FaceColor', [0.3 0.6 0.9]);
hold on;
plot(x_lap, f_lap_ana, 'r-', 'LineWidth', 2);
xlabel('x');
ylabel('PDF f_n(x)');
title('Laplacian Noise PDF');
legend('Simulated', 'Analytical', 'Location', 'best');
grid on;

subplot(1, 2, 2);
plot(x_lap, F_lap_sim, 'b-', 'LineWidth', 1.5);
hold on;
plot(x_lap, F_lap_ana, 'r--', 'LineWidth', 2);
xlabel('x');
ylabel('CDF F_n(x)');
title('Laplacian Noise CDF');
legend('Simulated', 'Analytical', 'Location', 'best');
grid on;

% Compute and display KS-like maximum absolute CDF error
max_cdf_error = max(abs(F_lap_sim - F_lap_ana));
fprintf('TASK 2: Laplacian CDF max absolute error = %.4e\n\n', max_cdf_error);

%% ===================== PRE-COMPUTE COMMON RANDOM VARIABLES =====================
% BPSK symbols
bits = 2*(randi([0 1], N_samples, 1)) - 1;   % +1 or -1 with equal probability

% Fading coefficients (independent of SNR, reused for all SNR points)
% Rayleigh fading: h_s ~ CN(0,1) => |h_s|^2 ~ Exp(1)
h_rayleigh = (randn(N_samples, 1) + 1j*randn(N_samples, 1)) / sqrt(2);

% Rician fading: h = sqrt(K/(K+1))*hd + sqrt(1/(K+1))*h_s
% E{|h|^2} = K/(K+1)*|hd|^2 + 1/(K+1)*1 = 1
h_los = sqrt(K/(K+1)) * hd;
h_scat = sqrt(1/(K+1)) * (randn(N_samples, 1) + 1j*randn(N_samples, 1)) / sqrt(2);
h_rician = h_los + h_scat;

% Random user deployment: distance d ~ f_d(d) = 2d/R^2, 0 <= d <= R
% Generate by d = R*sqrt(rand)
d_user = R * sqrt(rand(N_samples, 1));

% Combined fading + random deployment
h_rayleigh_rand = h_rayleigh ./ (d_user.^(nu/2));   % includes path loss d^{-nu/2}
h_rician_rand = h_rician ./ (d_user.^(nu/2));

%% ===================== TASK 3 & 4: BER SIMULATIONS =====================
% Preallocate BER arrays
BER_AWGN_sim = zeros(size(SNR_dB));
BER_AWGN_ana = zeros(size(SNR_dB));
BER_Lap_sim  = zeros(size(SNR_dB));
BER_Lap_ana  = zeros(size(SNR_dB));
BER_Ray_sim  = zeros(size(SNR_dB));
BER_Ray_ana  = zeros(size(SNR_dB));
BER_Ric_sim  = zeros(size(SNR_dB));
BER_Ric_ana  = zeros(size(SNR_dB));
BER_RayRand_sim = zeros(size(SNR_dB));
BER_RayRand_ana = zeros(size(SNR_dB));
BER_RicRand_sim = zeros(size(SNR_dB));
BER_RicRand_ana = zeros(size(SNR_dB));

for idx = 1:length(SNR_dB)
    gb = SNR_lin(idx);
    Eb = gb * N0;   % Energy per bit
    
    %% AWGN only
    n_awgn = sqrt(N0/2) * (randn(N_samples, 1) + 1j*randn(N_samples, 1));
    r_awgn = sqrt(Eb) * bits + n_awgn;
    BER_AWGN_sim(idx) = mean(real(r_awgn) .* bits < 0);
    BER_AWGN_ana(idx) = 0.5 * erfc(sqrt(gb));
    
    %% Laplacian impulse noise only
    u = rand(N_samples, 1);
    n_lap = -b_lap * sign(u - 0.5) .* log(1 - 2*abs(u - 0.5));
    r_lap = sqrt(Eb) * bits + n_lap;
    BER_Lap_sim(idx) = mean(real(r_lap) .* bits < 0);
    BER_Lap_ana(idx) = 0.5 * exp(-sqrt(2*gb));
    
    %% Rayleigh fading + AWGN
    n_ray = sqrt(N0/2) * (randn(N_samples, 1) + 1j*randn(N_samples, 1));
    r_ray = h_rayleigh * sqrt(Eb) .* bits + n_ray;
    r_eq_ray = real(conj(h_rayleigh) .* r_ray);
    BER_Ray_sim(idx) = mean(r_eq_ray .* bits < 0);
    BER_Ray_ana(idx) = 0.5 * (1 - sqrt(gb ./ (1 + gb)));
    
    %% Rician fading + AWGN
    n_ric = sqrt(N0/2) * (randn(N_samples, 1) + 1j*randn(N_samples, 1));
    r_ric = h_rician * sqrt(Eb) .* bits + n_ric;
    r_eq_ric = real(conj(h_rician) .* r_ric);
    BER_Ric_sim(idx) = mean(r_eq_ric .* bits < 0);
    % Analytical: using MGF-based integral representation
    BER_Ric_ana(idx) = compute_rician_ber(gb, K);
    
    %% Randomly deployed users over Rayleigh fading + AWGN
    n_rayrand = sqrt(N0/2) * (randn(N_samples, 1) + 1j*randn(N_samples, 1));
    r_rayrand = h_rayleigh_rand * sqrt(Eb) .* bits + n_rayrand;
    r_eq_rayrand = real(conj(h_rayleigh_rand) .* r_rayrand);
    BER_RayRand_sim(idx) = mean(r_eq_rayrand .* bits < 0);
    % Analytical via conditional expectation over distance (equivalent to PDF in Eq. 8)
    ber_cond_ray = @(d) 0.5 * (1 - sqrt(gb ./ (d.^nu + gb)));
    BER_RayRand_ana(idx) = (2/R^2) * integral(@(d) d .* ber_cond_ray(d), 0, R);
    
    %% Randomly deployed users over Rician fading + AWGN
    n_ricrand = sqrt(N0/2) * (randn(N_samples, 1) + 1j*randn(N_samples, 1));
    r_ricrand = h_rician_rand * sqrt(Eb) .* bits + n_ricrand;
    r_eq_ricrand = real(conj(h_rician_rand) .* r_ricrand);
    BER_RicRand_sim(idx) = mean(r_eq_ricrand .* bits < 0);
    % Analytical via conditional expectation: average Rician BER over d
    ber_cond_ric = @(d) arrayfun(@(dd) compute_rician_ber(gb/max(dd,eps)^nu, K), d);
    BER_RicRand_ana(idx) = (2/R^2) * integral(@(d) d .* ber_cond_ric(d), 0, R);
    
    fprintf('BER simulation completed for SNR = %d dB\n', SNR_dB(idx));
end

%% ===================== TASK 3 & 4: OUTAGE PROBABILITY =====================
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
    T = gamma_th / gb;
    
    %% Rayleigh fading (d = 1 m)
    inst_snr_ray = gb * abs(h_rayleigh).^2;
    Pout_Ray_sim(idx) = mean(log2(1 + inst_snr_ray) < C);
    Pout_Ray_ana(idx) = 1 - exp(-T);
    
    %% Rician fading (d = 1 m)
    inst_snr_ric = gb * abs(h_rician).^2;
    Pout_Ric_sim(idx) = mean(log2(1 + inst_snr_ric) < C);
    Pout_Ric_ana(idx) = 1 - marcumq(sqrt(2*K), sqrt(2*(K+1)*T));
    
    %% Randomly deployed users over Rayleigh fading
    inst_snr_rayrand = gb * abs(h_rayleigh).^2 ./ d_user.^nu;
    Pout_RayRand_sim(idx) = mean(log2(1 + inst_snr_rayrand) < C);
    % Analytical: P_out = 1 - (2/(nu*R^2)) * T^(-2/nu) * gamma(2/nu, T*R^nu)
    Pout_RayRand_ana(idx) = 1 - (2/(nu*R^2)) * T^(-2/nu) * gamma(2/nu) * gammainc(T*R^nu, 2/nu, 'lower');
    
    %% Randomly deployed users over Rician fading
    inst_snr_ricrand = gb * abs(h_rician).^2 ./ d_user.^nu;
    Pout_RicRand_sim(idx) = mean(log2(1 + inst_snr_ricrand) < C);
    % Analytical via conditional expectation over distance
    pout_cond_ric = @(d) 1 - marcumq(sqrt(2*K), sqrt(2*(K+1)*T*(d.^nu + eps)));
    Pout_RicRand_ana(idx) = (2/R^2) * integral(@(d) d .* pout_cond_ric(d), 0, R);
    
    fprintf('Outage simulation completed for SNR = %d dB\n', SNR_dB(idx));
end

fprintf('\n');

%% ===================== PERCENTAGE DEVIATION TABLES =====================
% Avoid division by zero in percentage error
pct_err = @(sim, ana) abs(sim - ana) ./ max(ana, eps) * 100;

fprintf('--- BER Percentage Deviation (Simulated vs Analytical) ---\n');
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

fprintf('\n--- Outage Probability Percentage Deviation (Simulated vs Analytical) ---\n');
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

%% ===================== TASK 5: PLOT 1 - BER NOISE ONLY =====================
figure('Name', 'BER vs SNR - Noise Only Channels', 'Position', [100 550 800 600]);
semilogy(SNR_dB, BER_AWGN_ana, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
semilogy(SNR_dB, max(BER_AWGN_sim, 1e-7), 'b--s', 'LineWidth', 1.5, 'MarkerSize', 6);
semilogy(SNR_dB, BER_Lap_ana, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 6);
semilogy(SNR_dB, max(BER_Lap_sim, 1e-7), 'r--s', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('Average SNR \gamma_b (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR for AWGN and Laplacian Noise Channels');
legend('AWGN Analytical', 'AWGN Simulated', ...
       'Laplacian Analytical', 'Laplacian Simulated', ...
       'Location', 'southwest');
axis([0 15 1e-5 1]);

%% ===================== TASK 5: PLOT 2 - BER FADING CHANNELS =====================
figure('Name', 'BER vs SNR - Fading Channels', 'Position', [150 500 800 600]);
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
xlabel('Average SNR \gamma_b (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR for Fading Channels (BPSK)');
legend('Rayleigh Analytical', 'Rayleigh Simulated', ...
       'Rician Analytical', 'Rician Simulated', ...
       'Random Rayleigh Analytical', 'Random Rayleigh Simulated', ...
       'Random Rician Analytical', 'Random Rician Simulated', ...
       'Location', 'southwest');
axis([0 15 1e-5 1]);

%% ===================== TASK 5: PLOT 3 - OUTAGE PROBABILITY =====================
figure('Name', 'Outage Probability vs SNR', 'Position', [200 450 800 600]);
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
xlabel('Average SNR \gamma_b (dB)');
ylabel('Outage Probability P_{out}');
title('Outage Probability vs SNR (C = 1.2 bps/Hz)');
legend('Rayleigh Analytical', 'Rayleigh Simulated', ...
       'Rician Analytical', 'Rician Simulated', ...
       'Random Rayleigh Analytical', 'Random Rayleigh Simulated', ...
       'Random Rician Analytical', 'Random Rician Simulated', ...
       'Location', 'southwest');
axis([0 15 1e-5 1]);

%% ===================== OPTIONAL: DEMONSTRATE ERROR REDUCTION =====================
% This section shows that increasing N reduces the percentage deviation.
% It compares AWGN BER at 10 dB for N = 1e5, 5e5, 1e6, 5e6.
fprintf('--- Demonstration: BER deviation decreases with more samples ---\n');
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
    fprintf('N = %7.0f : BER_sim = %.6e, Deviation = %.4f%%\n', Nt, ber_test, err_test);
end

%% ===================== HELPER FUNCTIONS =====================
function ber = compute_rician_ber(gb, K)
    % Compute analytical BER for BPSK over Rician fading with average SNR gb.
    % Uses the MGF-based integral representation:
    % P_b = (1/pi) * integral_0^(pi/2) M_gamma(-1/sin^2(theta)) d(theta)
    % where M_gamma(s) = (K+1)/(K+1 - s*gb) * exp(-K*gb*s/(K+1 - s*gb)).
    % For very high SNR, BER is numerically zero.
    if gb > 1e4
        ber = 0;
        return;
    end
    fun = @(theta) (1/pi) * ((K+1) .* sin(theta).^2) ./ ((K+1).*sin(theta).^2 + gb) ...
                   .* exp(-K*gb ./ ((K+1).*sin(theta).^2 + gb));
    ber = integral(fun, 0, pi/2);
end
