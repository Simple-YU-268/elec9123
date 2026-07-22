function [rho, N_theory, T_theory] = mm1_analysis(lambda_hourly, service_minutes)
% mm1_analysis: computes theoretical metrics for an M/M/1 queue.
% Inputs:
%   lambda_hourly   - arrival rate in customers per hour
%   service_minutes - mean service time in minutes
% Outputs:
%   rho      - server utilisation
%   N_theory - mean number of customers in the system (N)
%   T_theory - mean time in system per customer (Little's law)

% Convert service time to hours so units match lambda
mu = 60 / service_minutes;          % service rate (customers per hour)
rho = lambda_hourly / mu;

if rho >= 1
    error('rho = %.4f >= 1; M/M/1 system unstable.', rho);
end

N_theory = rho / (1 - rho);
T_theory = N_theory / lambda_hourly; % Little's law: N = lambda*T
end
