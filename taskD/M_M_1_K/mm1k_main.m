function mm1k_main()
% mm1k_main: M/M/1/K comparison for service times 5 and 6 minutes.

lambda_hourly = 10;
K = 10;
service_times = [5, 6];  % minutes
num_cases = length(service_times);

event_counts = [1e4, 1e5, 1e6];
num_runs = length(event_counts);

fprintf('=== M/M/1/K Queue (K = %d) ===\n', K);

for c = 1:num_cases
    s = service_times(c);
    [rho, N_theory, T_theory, P_block] = mm1k_analysis(lambda_hourly, s, K);
    T_theory_min = T_theory * 60;

    fprintf('\n--- Service time = %.1f min ---\n', s);
    fprintf('rho = %.4f, P_block = %.4f\n', rho, P_block);
    fprintf('Theory: N = %.4f, T = %.4f min\n', N_theory, T_theory_min);

    N_sim = zeros(num_runs, 1);
    T_sim = zeros(num_runs, 1);
    P_block_sim = zeros(num_runs, 1);

    for r = 1:num_runs
        [N_sim(r), T_sim(r), P_block_sim(r), ~] = mm1k_simulation(lambda_hourly, s, K, event_counts(r), 42);
        fprintf('Sim %6.0e: N=%.4f, T=%.4f min, P_block=%.4f\n', event_counts(r), N_sim(r), T_sim(r), P_block_sim(r));
    end

    figure('Name', sprintf('M/M/1/K S=%.1f min', s));

    subplot(1, 3, 1);
    bar([N_theory, N_sim']);
    title('Mean N');
    ylabel('N');
    set(gca, 'XTickLabel', {'Theory','1e4','1e5','1e6'});

    subplot(1, 3, 2);
    bar([T_theory_min, T_sim']);
    title('Mean T (min)');
    ylabel('T (min)');
    set(gca, 'XTickLabel', {'Theory','1e4','1e5','1e6'});

    subplot(1, 3, 3);
    bar([P_block, P_block_sim']);
    title('Blocking probability');
    ylabel('P_{block}');
    set(gca, 'XTickLabel', {'Theory','1e4','1e5','1e6'});
end
end
