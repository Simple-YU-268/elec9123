function mmm_main()
% mmm_main: Two-class M/M/m airport check-in scenario.

% Business class
lambda_b = 20;   % per hour
m_b = 2;
s_b = 2.5;     % minutes

% Economy class
lambda_e = 100; % per hour
m_e = 5;
s_e = 2.5;      % minutes

event_counts = [1e4, 1e5, 1e6];
num_runs = length(event_counts);

fprintf('=== Two-class M/M/m Airport Check-in ===\n');

% Business class
[Nth_b, Tth_b, Pwth_b] = mmm_analysis(lambda_b, s_b, m_b);
Tth_b_min = Tth_b * 60;
fprintf('\nBusiness (lambda=%d/h, m=%d, S=%.1f min):\n', lambda_b, m_b, s_b);
fprintf('Theory: N=%.4f, T=%.4f min, P_wait=%.4f\n', Nth_b, Tth_b_min, Pwth_b);

N_sim_b = zeros(num_runs, 1); T_sim_b = zeros(num_runs, 1); Pw_sim_b = zeros(num_runs, 1);
for r = 1:num_runs
    [N_sim_b(r), T_sim_b(r), Pw_sim_b(r), ~] = mmm_simulation(lambda_b, s_b, m_b, event_counts(r), 42);
    fprintf('Sim %6.0e: N=%.4f, T=%.4f min, P_wait=%.4f\n', event_counts(r), N_sim_b(r), T_sim_b(r), Pw_sim_b(r));
end

% Economy class
[Nth_e, Tth_e, Pwth_e] = mmm_analysis(lambda_e, s_e, m_e);
Tth_e_min = Tth_e * 60;
fprintf('\nEconomy (lambda=%d/h, m=%d, S=%.1f min):\n', lambda_e, m_e, s_e);
fprintf('Theory: N=%.4f, T=%.4f min, P_wait=%.4f\n', Nth_e, Tth_e_min, Pwth_e);

N_sim_e = zeros(num_runs, 1); T_sim_e = zeros(num_runs, 1); Pw_sim_e = zeros(num_runs, 1);
for r = 1:num_runs
    [N_sim_e(r), T_sim_e(r), Pw_sim_e(r), ~] = mmm_simulation(lambda_e, s_e, m_e, event_counts(r), 42);
    fprintf('Sim %6.0e: N=%.4f, T=%.4f min, P_wait=%.4f\n', event_counts(r), N_sim_e(r), T_sim_e(r), Pw_sim_e(r));
end

% Combined bar chart
figure('Name','Two-class M/M/m Comparison');

subplot(1, 3, 1);
bar([Nth_b, mean(N_sim_b); Nth_e, mean(N_sim_e)]);
set(gca,'XTickLabel', {'Business','Economy'});
legend('Theory','Simulation','Location','best');
title('Mean N');
ylabel('N');

subplot(1, 3, 2);
bar([Tth_b_min, mean(T_sim_b); Tth_e_min, mean(T_sim_e)]);
set(gca,'XTickLabel', {'Business','Economy'});
legend('Theory','Simulation','Location','best');
title('Mean T (min)');
ylabel('T (min)');

subplot(1, 3, 3);
bar([Pwth_b, mean(Pw_sim_b); Pwth_e, mean(Pw_sim_e)]);
set(gca,'XTickLabel', {'Business','Economy'});
legend('Theory','Simulation','Location','best');
title('P_{wait}');
ylabel('P_{wait}');
end
