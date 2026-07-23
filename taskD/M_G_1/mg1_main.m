function mg1_main()
% mg1_main: M/G/1 comparison for exponential, uniform, deterministic services.

lambda_hourly = 10;
mean_service = 5;  % minutes
dist_types = {'exponential', 'uniform', 'deterministic'};
event_counts = [1e4, 1e5, 1e6];
num_dists = length(dist_types);
num_runs = length(event_counts);

fprintf('=== M/G/1 Queue ===\n');

N_theory = zeros(num_dists, 1);
T_theory = zeros(num_dists, 1);

for d = 1:num_dists
    dt = dist_types{d};
    [N_theory(d), T_theory(d)] = mg1_analysis(lambda_hourly, mean_service, dt);
    T_theory_min = T_theory(d) * 60;
    fprintf('%s: N_theory=%.4f, T_theory=%.4f min\n', dt, N_theory(d), T_theory_min);

    N_sim = zeros(num_runs, 1);
    T_sim = zeros(num_runs, 1);
    for r = 1:num_runs
        [N_sim(r), T_sim(r), ~] = mg1_simulation(lambda_hourly, mean_service, dt, event_counts(r), 42);
        fprintf('  Sim %6.0e: N=%.4f, T=%.4f min\n', event_counts(r), N_sim(r), T_sim(r));
    end

    figure('Name', sprintf('M/G/1 %s', dt));

    subplot(1,2,1);
    bar([N_theory(d), N_sim']);
    title(sprintf('Mean N (%s)', dt));
    ylabel('N');
    set(gca, 'XTickLabel', {'Theory','1e4','1e5','1e6'});

    subplot(1,2,2);
    bar([T_theory_min, T_sim']);
    title(sprintf('Mean T (min) (%s)', dt));
    ylabel('T (min)');
    set(gca, 'XTickLabel', {'Theory','1e4','1e5','1e6'});
end
end
