function mm1_main()
% mm1_main: M/M/1 - compares theoretical and simulation results.

lambda_hourly = 10;        % customers per hour
service_minutes = 5;       % mean service time in minutes

[rho, N_theory, T_theory] = mm1_analysis(lambda_hourly, service_minutes);

% Convert T_theory to minutes for plotting
T_theory_minutes = T_theory * 60;

fprintf('=== M/M/1 Queue ===\n');
fprintf('Arrival rate  = %.2f cust/hour\n', lambda_hourly);
fprintf('Mean service  = %.2f minutes\n', service_minutes);
fprintf('Utilisation rho = %.4f\n', rho);
fprintf('Theoretical N   = %.4f\n', N_theory);
fprintf('Theoretical T   = %.4f minutes\n', T_theory_minutes);

% Run simulation with increasing fidelity
event_counts = [1e4, 1e5, 1e6];
num_runs = length(event_counts);
N_sim = zeros(num_runs, 1);
T_sim = zeros(num_runs, 1);

for i = 1:num_runs
    [N_sim(i), T_sim(i), ~] = mm1_simulation(lambda_hourly, service_minutes, event_counts(i), 42);
    fprintf('Sim (%6.0e events): N=%.4f, T=%.4f min\n', event_counts(i), N_sim(i), T_sim(i));
end

% Bar chart comparison
figure('Name','M/M/1 Comparison');

subplot(1,2,1);
bar([N_theory, N_sim']);
hold on;
yline(N_theory, 'r--', 'LineWidth', 1.5);
hold off;
title('Mean number in system N');
ylabel('N');
set(gca, 'XTickLabel', {'Theory','1e4','1e5','1e6'});

subplot(1,2,2);
bar([T_theory_minutes, T_sim']);
hold on;
yline(T_theory_minutes, 'r--', 'LineWidth', 1.5);
hold off;
title('Mean time in system T (minutes)');
ylabel('T (min)');
set(gca, 'XTickLabel', {'Theory','1e4','1e5','1e6'});

end
