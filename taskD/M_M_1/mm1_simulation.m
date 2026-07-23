function [N_sim, T_sim, num_events] = mm1_simulation(lambda_hourly, service_minutes, target_events, seed)
% mm1_simulation: discrete-event simulation of a stable M/M/1 queue.
% Inputs:
%   lambda_hourly   - arrival rate in customers per hour
%   service_minutes - mean service time in minutes
%   target_events   - minimum number of events to simulate
%   seed            - optional RNG seed (default: rng('shuffle'))
% Outputs:
%   N_sim     - time-averaged number of customers in system
%   T_sim     - averaged total system time per served customer (minutes)
%   num_events- number of events actually simulated

if nargin < 4 || isempty(seed)
    rng('shuffle');
else
    rng(seed);
end

lambda = lambda_hourly / 60;  % per minute
mu = 1 / service_minutes;      % per minute

% Pre-allocate storage
max_customers = 2 * target_events;
arr_times = zeros(max_customers, 1);
dep_times = zeros(max_customers, 1);

% Event state
t = 0;
next_arrival = exprnd(1 / lambda);
next_departure = inf;
server_busy = false;
queue = [];             % FIFO list of customer indices waiting
next_customer_id = 0;   % customer currently in service (if any)
last_customer = 0;      % total arrivals so far

% Time-average N
num_in_system = 0;
last_event_time = 0;
N_integral = 0;

num_events = 0;
num_served = 0;

while num_events < target_events
    if next_arrival <= next_departure
        % Arrival event
        dt = next_arrival - t;
        N_integral = N_integral + num_in_system * dt;
        t = next_arrival;
        last_event_time = t;
        num_events = num_events + 1;

        last_customer = last_customer + 1;
        arr_times(last_customer) = t;

        if ~server_busy
            % Start service immediately
            server_busy = true;
            next_customer_id = last_customer;
            dep_times(next_customer_id) = t + exprnd(1 / mu);
            next_departure = dep_times(next_customer_id);
            num_in_system = 1;
        else
            % Join queue
            queue(end + 1) = last_customer; %#ok<AGROW>
            num_in_system = num_in_system + 1;
        end

        next_arrival = t + exprnd(1 / lambda);

    else
        % Departure event
        dt = next_departure - t;
        N_integral = N_integral + num_in_system * dt;
        t = next_departure;
        last_event_time = t;
        num_events = num_events + 1;
        num_served = num_served + 1;
        num_in_system = num_in_system - 1;

        if ~isempty(queue)
            % Start service for next customer in FIFO order
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

% Time-average over the simulated interval [0, t]
N_sim = N_integral / t;

% Sojourn times for served customers (exclude first 10% as warm-up)
valid_served = (1:num_served) + floor(0.1 * num_served);
valid_served = valid_served(valid_served <= num_served);
T = dep_times(1:num_served) - arr_times(1:num_served);
T_sim = mean(T(valid_served)); % simulation clock is already in minutes
end
