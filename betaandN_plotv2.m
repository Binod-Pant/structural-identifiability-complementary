% Plot SIR Model Results from .mat File
clear; clc; close all;

%% Load Results
filename = 'SIR_estimation_results.mat';
if ~exist(filename, 'file')
    error('File %s not found! Please update the filename.', filename);
end

load(filename);
fprintf('Loaded results from: %s\n', filename);

%% Extract Data
% True parameters and curves
beta_true = results.true_params.beta;
gamma_true = results.true_params.gamma;
N_true = results.true_params.N;
S0_true = results.true_params.S0;
I0_true = results.true_params.I0;
R0_true = results.true_params.R0;

t = results.time;
S_true = results.true_curves.S;
I_true = results.true_curves.I;
R_true = results.true_curves.R;
incidence_true = results.true_curves.incidence;

% Accepted results
accepted_count = results.accepted_count;
total_iterations = results.total_iterations;
accepted_params = results.accepted_params;
accepted_S = results.accepted_curves.S;
accepted_I = results.accepted_curves.I;
accepted_R = results.accepted_curves.R;
accepted_incidence = results.accepted_curves.incidence;
accepted_errors = results.accepted_errors;
sse_threshold = results.sse_threshold;

% R_naught data
accepted_R_naught = results.accepted_R_naught;
R_naught_true = beta_true * S0_true / N_true * (1/gamma_true);

fprintf('Accepted solutions: %d out of %d iterations\n', accepted_count, total_iterations);

%% Plotting



% First row: Incidence and SIR curves
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
nexttile;
h1 = plot(t, accepted_incidence{1}, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1); hold on;
for i = 2:accepted_count
   plot(t, accepted_incidence{i}, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
end
h2 = plot(t, incidence_true, 'r--', 'LineWidth', 3);
xlabel('Time');
title(' (a) Incidence ($\beta \times S \times I/N$)', 'Interpreter', 'latex', 'FontSize', 14);
legend([h2, h1], {'Observation','Fitted'}, 'Location', 'northeast','FontSize', 11);
grid on;

%% plot beta vs N
nexttile
if accepted_count > 0
    beta_estimates = accepted_params(:,1);
    N_estimates = sum(accepted_params(:,2:4), 2);
h1 = plot(N_estimates, beta_estimates, '*', 'MarkerSize', 8, 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5]);
hold on
h2 = plot(N_true, beta_true, 'co', 'MarkerSize', 10, 'MarkerFaceColor', 'green');
xlabel('Population ($N$)', 'Interpreter', 'latex');
    ylabel('Transmission rate ($\beta$)', 'Interpreter', 'latex', 'FontSize', 14);
    title('(b) $\beta$ vs N', 'Interpreter', 'latex', 'FontSize', 14);
end

  % Add theoretical line: beta/N = constant
   N_range = linspace(min([N_estimates; N_true]) * 0.9, max([N_estimates; N_true]) * 1.1, 100);
   theoretical_beta = (beta_true/N_true) * N_range;
   h3 = plot(N_range, theoretical_beta, 'b--', 'LineWidth', 2);
    legend([h2, h1,h3], {'True value','Estimated value', 'Theoretical relationship'}, 'Location', 'northwest','FontSize', 11);
grid on;
xlim([0.7, 2])


%% Plot R_e(t) - Effective Reproduction Number
nexttile

% Calculate R_e(t) for all accepted solutions
R_e_curves = cell(accepted_count, 1);
for i = 1:accepted_count
    beta_i = accepted_params(i, 1);
    N_i = sum(accepted_params(i, 2:4));
    S_i = accepted_S{i};
    R_e_curves{i} = beta_i * S_i / (N_i * gamma_true);
end

% Calculate true R_e(t)
R_e_true = beta_true * S_true / (N_true * gamma_true);

% Plot all accepted R_e(t) curves
h1 = plot(t, R_e_curves{1}, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1); hold on;
for i = 2:accepted_count
    plot(t, R_e_curves{i}, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
end

% Plot true R_e(t)
h2 = plot(t, R_e_true, 'g--', 'LineWidth', 3);

% Add threshold line at R_e = 1
h3 = plot([t(1), t(end)], [1, 1], 'k--', 'LineWidth', 2);

xlabel('Time');
ylabel('$\mathcal{R}_e(t)$', 'Interpreter', 'latex', 'FontSize', 14);
title('(c) Effective reproduction number ($\mathcal{R}_e(t)$)', 'Interpreter', 'latex', 'FontSize', 14);
legend([h2, h1, h3], {'True', 'Estimated', '$\mathcal{R}_e(t)=1$'},'Interpreter', 'latex', 'Location', 'northeast','FontSize', 11);
grid on;
%%

%% plot S, I, R

nexttile
h1 = plot(t, accepted_S{1}, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1); hold on;
for i = 2:accepted_count
   plot(t, accepted_S{i}, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
end
h2 = plot(t, S_true, 'g--', 'LineWidth', 3);
xlabel('Time');
title('(d) Susceptible (S)', 'Interpreter', 'latex', 'FontSize', 14);
legend([h2, h1], {'True','Estimated'}, 'Location', 'northeast','FontSize', 11);
grid on;

nexttile
h1 = plot(t, accepted_I{1}, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1); hold on;
for i = 2:accepted_count
   plot(t, accepted_I{i}, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
end
h2 = plot(t, I_true, 'g--', 'LineWidth', 3);
xlabel('Time');
title('(e) Infectious (I)', 'Interpreter', 'latex', 'FontSize', 14);
legend([h2, h1], {'True','Estimated'}, 'Location', 'northeast','FontSize', 11);
grid on;

nexttile
h1 = plot(t, accepted_R{1}, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1); hold on;
for i = 2:accepted_count
    plot(t, accepted_R{i}, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
end
h2 = plot(t, R_true, 'g--', 'LineWidth', 3);
xlabel('Time');
title('(f) Recovered (R)', 'Interpreter', 'latex', 'FontSize', 14);
legend([h2, h1], {'True','Estimated'}, 'Location', 'southeast','FontSize', 11);
grid on;
ylim([0 1])

% % Second row: Parameter estimation results
% subplot(3,4,5)
% if accepted_count > 0
%     beta_estimates = accepted_params(:,1);
%     plot(1:accepted_count, beta_estimates, 'bo-', 'LineWidth', 2); hold on;
%     %plot([1, accepted_count], [beta_true, beta_true], 'ro', 'LineWidth', 2);
%     plot(1:accepted_count, beta_true * ones(1, accepted_count), 'ro', 'LineWidth', 2);
%     xlabel('Accepted Solution #'); 
%     %ylabel('$\beta$ estimate', 'Interpreter', 'latex');
%     title('$\beta$ estimate', 'Interpreter', 'latex');
%     %title('Beta Estimation');
% end
% grid on;

% subplot(3,4,6)
% if accepted_count > 0
%     N_estimates = sum(accepted_params(:,2:4), 2);
%     plot(1:accepted_count, N_estimates, 'bo-', 'LineWidth', 2); hold on;
%     plot([1, accepted_count], [N_true, N_true], 'r-', 'LineWidth', 2);
%     %plot(1:accepted_count, beta_true * ones(1, accepted_count), 'ro', 'LineWidth', 2);
%     xlabel('Accepted Solution #'); 
%     %ylabel('N estimate');
%     title('Population ($N$) estimate', 'Interpreter', 'latex');
% end
% grid on;



% subplot(3,4,8)
% if accepted_count > 0
%     semilogy(1:accepted_count, accepted_errors, 'bo-', 'LineWidth', 2); hold on;
%     plot([1, accepted_count], [sse_threshold, sse_threshold], 'g--', 'LineWidth', 2);
%     xlabel('Accepted Solution #'); 
%     %ylabel('SSE');
%     title('SSE Values');
% end
% grid on;

% % Third row: Initial conditions
% subplot(3,4,9)
% if accepted_count > 0
%     S0_estimates = accepted_params(:,2);
%     plot(1:accepted_count, S0_estimates, 'bo-', 'LineWidth', 2); hold on;
%     plot([1, accepted_count], [S0_true, S0_true], 'r-', 'LineWidth', 2);
%     xlabel('Accepted Solution #'); 
%     %ylabel('S(0) estimate');
%     title('Initial Susceptible S(0)', 'Interpreter', 'latex');
% end
% grid on;
% 
% subplot(3,4,10)
% if accepted_count > 0
%     I0_estimates = accepted_params(:,3);
%     plot(1:accepted_count, I0_estimates, 'bo-', 'LineWidth', 2); hold on;
%     plot([1, accepted_count], [I0_true, I0_true], 'r-', 'LineWidth', 2);
%     xlabel('Accepted Solution #'); 
%     %ylabel('I(0) estimate');
%     title('Initial Infected I(0)', 'Interpreter', 'latex');
% end
% grid on;
% 
% subplot(3,4,11)
% if accepted_count > 0
%     R0_estimates = accepted_params(:,4);
%     plot(1:accepted_count, R0_estimates, 'bo-', 'LineWidth', 2); hold on;
%     plot([1, accepted_count], [R0_true, R0_true], 'r-', 'LineWidth', 2);
%     xlabel('Accepted Solution #'); 
%     %ylabel('R(0) estimate');
%     title('Initial Recovered R(0)', 'Interpreter', 'latex');
% end
% grid on;




% Summary and legend box
% nexttile
% axis off;
% % Legend
% text(0.5, 0.7, 'LEGEND', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(0.5, 0.55, '- -  Observation', 'Units', 'normalized', 'FontSize', 12, 'Color', 'magenta', 'HorizontalAlignment', 'center');
% text(0.5, 0.45, '- - True SIR curves', 'Units', 'normalized', 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
% text(0.5, 0.35, '—  Estimated curves', 'Units', 'normalized', 'FontSize', 12, 'Color', [0.5 0.5 0.5], 'HorizontalAlignment', 'center');
% text(0.5, 0.25, 'o  True values', 'Units', 'normalized', 'FontSize', 12, 'Color', 'red', 'HorizontalAlignment', 'center');
% text(0.5, 0.15, 'o  Estimated values', 'Units', 'normalized', 'FontSize', 12, 'Color', 'blue', 'HorizontalAlignment', 'center');
% text(0.5, 0.05, '- -  SSE Threshold', 'Units', 'normalized', 'FontSize', 12, 'Color', 'green', 'HorizontalAlignment', 'center');
% 

fprintf('Plot completed successfully!\n');

% paste this in command window
%print(gcf, 'figure1.png', '-dpng', '-r300');
%print(gcf, 'figure1.pdf', '-dpdf', '-r300');