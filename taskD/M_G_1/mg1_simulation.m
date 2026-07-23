function [N_sim, T_sim, num_events] = mg1_simulation(lambda_hourly, service_minutes, dist_type, target_events, seed)
% mg1_simulation: M/G/1 discrete-event simulation with arbitrary service dist.

if nargin < 5 || isempty(seed)
    rng('shuffle');
else
    rng(seed);
end

lambda = lambda_hourly / 60;  % per minute
mu = 1 / service_minutes;      % only used for exponential case

max_customers = 2 * target_events;
arr_times = zeros(max_customers, 1);
dep_times = zeros(max_customers, 1);

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
        dt = next_arrival - t;
        N_integral = N_integral + num_in_system * dt;
        t = next_arrival;
        num_events = num_events + 1;

        last_customer = last_customer + 1;
        arr_times(last_customer) = t;

        if ~server_busy
            server_busy = true;
            next_customer_id = last_customer;
            S = generate_service(dist_type, service_minutes, mu);
            dep_times(next_customer_id) = t + S;
            next_departure = dep_times(next_customer_id);
            num_in_system = 1;
        else
            queue(end + 1) = last_customer; %#ok<AGROW>
            num_in_system = num_in_system + 1;
        end

        next_arrival = t + exprnd(1 / lambda);
    else
        dt = next_departure - t;
        N_integral = N_integral + num_in_system * dt;
        t = next_departure;
        num_events = num_events + 1;
        num_served = num_served + 1;
        num_in_system = num_in_system - 1;

        if ~isempty(queue)
            next_customer_id = queue(1);
            queue(1) = [];
            S = generate_service(dist_type, service_minutes, mu);
            dep_times(next_customer_id) = t + S;
            next_departure = dep_times(next_customer_id);
            % queued customer stays in the system; count unchanged
        else
            server_busy = false;
            next_departure = inf;
        end
    end
end

N_sim = N_integral / t;
T = dep_times(1:num_served) - arr_times(1:num_served);
warm_up = max(1, floor(0.1 * num_served));
T_sim = mean(T(warm_up:num_served)); % simulation clock is already in minutes
end

function S = generate_service(dist_type, mean_minutes, mu)
    switch lower(dist_type)
        case 'exponential'
            S = exprnd(1 / mu);
        case 'uniform'
            a = 2.5;
            b = 7.5;
            S = a + (b - a) * rand;
        case 'deterministic'
            S = mean_minutes;
        otherwise
            error('Unknown dist_type');
    end
end
