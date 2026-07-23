function [N_sim, T_sim, P_block_sim, num_events] = mm1k_simulation(lambda_hourly, service_minutes, K, target_events, seed)
% mm1k_simulation: discrete-event simulation of M/M/1/K.
% Inputs:
%   lambda_hourly, service_minutes, K as above
%   target_events - minimum number of arrival/departure events
%   seed          - optional RNG seed

if nargin < 5 || isempty(seed)
    rng('shuffle');
else
    rng(seed);
end

lambda = lambda_hourly / 60;
mu = 1 / service_minutes;

max_customers = 2 * target_events;
arr_times = zeros(max_customers, 1);
dep_times = zeros(max_customers, 1);
blocked_count = 0;

server_busy = false;
queue = [];
next_customer_id = 0;
last_customer = 0;
num_in_system = 0;

t = 0;
next_arrival = exprnd(1 / lambda);
next_departure = inf;

N_integral = 0;
num_events = 0;
num_served = 0;

while num_events < target_events
    if next_arrival <= next_departure
        % Arrival event
        dt = next_arrival - t;
        N_integral = N_integral + num_in_system * dt;
        t = next_arrival;
        num_events = num_events + 1;

        if num_in_system < K
            % Accept customer
            last_customer = last_customer + 1;
            arr_times(last_customer) = t;

            if ~server_busy
                server_busy = true;
                next_customer_id = last_customer;
                dep_times(next_customer_id) = t + exprnd(1 / mu);
                next_departure = dep_times(next_customer_id);
                num_in_system = 1;
            else
                queue(end + 1) = last_customer; %#ok<AGROW>
                num_in_system = num_in_system + 1;
            end
        else
            % System full: block
            blocked_count = blocked_count + 1;
        end

        next_arrival = t + exprnd(1 / lambda);
    else
        % Departure event
        dt = next_departure - t;
        N_integral = N_integral + num_in_system * dt;
        t = next_departure;
        num_events = num_events + 1;
        num_served = num_served + 1;
        num_in_system = num_in_system - 1;

        if ~isempty(queue)
            next_customer_id = queue(1);
            queue(1) = [];
            dep_times(next_customer_id) = t + exprnd(1 / mu);
            next_departure = dep_times(next_customer_id);
            % queued customer stays in the system; count unchanged
        else
            server_busy = false;
            next_departure = inf;
        end
    end
end

N_sim = N_integral / t;

% Warm-up exclusion for sojourn time statistics
T = dep_times(1:num_served) - arr_times(1:num_served);
warm_up = max(1, floor(0.1 * num_served));
T_sim = mean(T(warm_up:num_served)); % simulation clock is already in minutes

% Blocking probability: total arrivals attempted = last_customer + blocked_count
% Arrival attempts include all arrival events (accepted + blocked)
P_block_sim = blocked_count / (last_customer + blocked_count);
end
