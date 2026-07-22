function [N_sim, T_sim, P_wait_sim, num_events] = mmm_simulation(lambda_hourly, service_minutes, m, target_events, seed)
% mmm_simulation: discrete-event simulation for M/M/m with m parallel servers.

if nargin < 5 || isempty(seed)
    rng('shuffle');
else
    rng(seed);
end

lambda = lambda_hourly / 60;  % per minute
mu = 1 / service_minutes;     % per minute

max_customers = 2 * target_events;
arr_times = zeros(max_customers, 1);
dep_times = zeros(max_customers, 1);

% Event list: one departure time per server (inf if idle)
next_departures = inf(m, 1);
server_busy = false(m, 1);
server_customer = zeros(m, 1); % which customer each server is serving

queue = [];       % FIFO waiting customers
last_customer = 0;
num_in_system = 0;

t = 0;
next_arrival = exprnd(1 / lambda);

N_integral = 0;
num_events = 0;
num_served = 0;
num_waited = 0;  % customers that had to wait before service

while num_events < target_events
    [min_departure, dep_server] = min(next_departures);

    if next_arrival <= min_departure
        % Arrival event
        dt = next_arrival - t;
        N_integral = N_integral + num_in_system * dt;
        t = next_arrival;
        num_events = num_events + 1;

        last_customer = last_customer + 1;
        arr_times(last_customer) = t;

        free_servers = find(~server_busy, 1);
        if ~isempty(free_servers)
            % Go directly to service
            s = free_servers(1);
            server_busy(s) = true;
            server_customer(s) = last_customer;
            dep_times(last_customer) = t + exprnd(1 / mu);
            next_departures(s) = dep_times(last_customer);
        else
            % All servers busy, join queue
            num_waited = num_waited + 1;
            queue(end + 1) = last_customer; %#ok<AGROW>
        end
        num_in_system = num_in_system + 1;

        next_arrival = t + exprnd(1 / lambda);
    else
        % Departure event from dep_server
        dt = min_departure - t;
        N_integral = N_integral + num_in_system * dt;
        t = min_departure;
        num_events = num_events + 1;
        num_served = num_served + 1;
        num_in_system = num_in_system - 1;

        cid = server_customer(dep_server);
        dep_times(cid) = t;

        if ~isempty(queue)
            % Start service for the next queued customer
            next_cid = queue(1);
            queue(1) = [];
            server_customer(dep_server) = next_cid;
            dep_times(next_cid) = t + exprnd(1 / mu);
            next_departures(dep_server) = dep_times(next_cid);
            num_in_system = num_in_system + 1; % customer moved from queue to service
        else
            server_busy(dep_server) = false;
            server_customer(dep_server) = 0;
            next_departures(dep_server) = inf;
        end
    end
end

N_sim = N_integral / t;
T = dep_times(1:num_served) - arr_times(1:num_served);
warm_up = max(1, floor(0.1 * num_served));
T_sim = mean(T(warm_up:num_served)) * 60;
P_wait_sim = num_waited / last_customer;
end
