function [N_theory, T_theory] = mg1_analysis(lambda_hourly, service_minutes, dist_type)
% mg1_analysis: theoretical M/G/1 mean customers and time via P-K formula.
% Inputs:
%   lambda_hourly   - arrivals per hour
%   service_minutes - mean service time in minutes
%   dist_type       - 'exponential', 'uniform', or 'deterministic'
% Outputs:
%   N_theory - mean number in system
%   T_theory - mean time in system (hours)

lambda = lambda_hourly / 60;  % per minute
E_X = service_minutes;        % minutes
rho = lambda * E_X;

if rho >= 1
    error('rho = %.4f >= 1; M/G/1 unstable.', rho);
end

switch lower(dist_type)
    case 'exponential'
        E_X2 = 2 * (E_X^2);               % E[X^2] for exponential
    case 'uniform'
        a = 2.5;
        b = 7.5;
        % mean = (a+b)/2 = 5; variance = (b-a)^2/12
        E_X2 = (a^2 + a*b + b^2) / 3;     % E[X^2] for uniform(a,b)
    case 'deterministic'
        E_X2 = E_X^2;
    otherwise
        error('Unknown distribution type: %s', dist_type);
end

% Pollaczek-Khinchin mean number in system
N_theory = rho + (lambda^2 * E_X2) / (2 * (1 - rho));
T_theory = N_theory / lambda_hourly;  % hours
end
