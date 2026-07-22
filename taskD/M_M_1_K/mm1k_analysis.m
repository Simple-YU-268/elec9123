function [rho, N_theory, T_theory, P_block] = mm1k_analysis(lambda_hourly, service_minutes, K)
% mm1k_analysis: theoretical metrics for M/M/1/K queue.
% Inputs:
%   lambda_hourly   - arrival rate in customers per hour
%   service_minutes - mean service time in minutes
%   K               - system capacity (max customers in system)
% Outputs:
%   rho     - traffic intensity (lambda/mu)
%   N_theory- mean number of customers in system
%   T_theory- mean time in system per accepted customer
%   P_block - probability that an arriving customer is blocked

mu = 60 / service_minutes;  % customers per hour
rho = lambda_hourly / mu;

if rho == 1
    p0 = 1 / (K + 1);
    pn = p0 * ones(1, K + 1);
else
    % Steady-state probability of n customers in system
    n = 0:K;
    pn_unnorm = rho.^n;
    p0 = 1 / sum(pn_unnorm);
    pn = p0 * pn_unnorm;
end

N_theory = sum(n .* pn);
P_block = pn(end);          % blocking probability = P(K customers in system)

% Effective arrival rate
lambda_eff = lambda_hourly * (1 - P_block);

% Throughput/conservation: lambda_eff = mu * (1 - p0) = mu * server utilisation
% Mean time in system via Little's law applied to accepted customers
T_theory = N_theory / lambda_eff;  % hours
end
