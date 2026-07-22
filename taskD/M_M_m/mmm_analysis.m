function [N_theory, T_theory, P_wait_theory] = mmm_analysis(lambda_hourly, service_minutes, m)
% mmm_analysis: Erlang-C style M/M/m metrics for a single class of servers.
% Inputs:
%   lambda_hourly   - arrivals per hour
%   service_minutes - mean service time in minutes
%   m               - number of identical servers
% Outputs:
%   N_theory     - mean number of customers in system
%   T_theory     - mean time in system (hours)
%   P_wait_theory- probability an arriving customer must wait

lambda = lambda_hourly / 60;  % per minute
mu = 1 / service_minutes;     % per minute
rho = lambda / (m * mu);      % per-server utilisation

if rho >= 1
    error('rho = %.4f >= 1; M/M/%d unstable.', rho, m);
end

% Erlang-C formula: P(wait) = ( (m*rho)^m / m! ) / ( (1-rho) * sum_{k=0}^{m-1} ... + ... )
% Compute a0 = p0^{-1}
% a0 = sum_{n=0}^{m-1} (m*rho)^n / n! + (m*rho)^m / m! * 1/(1-rho)
a = lambda / mu;  % offered load in erlangs
term1 = 0;
for n = 0:(m-1)
    term1 = term1 + (a^n) / factorial(n);
end
term2 = (a^m) / factorial(m) * (1 / (1 - rho));
inv_p0 = term1 + term2;
p0 = 1 / inv_p0;

P_wait_theory = p0 * (a^m) / factorial(m) * (1 / (1 - rho));

% Mean queue length Lq = P_wait * rho / (1 - rho)
Lq = P_wait_theory * rho / (1 - rho);
% Mean in service = a = lambda/mu
N_theory = Lq + a;
% Mean time in system via Little's law
T_theory = N_theory / lambda_hourly;  % hours
end
