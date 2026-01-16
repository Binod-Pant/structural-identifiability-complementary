clc; clear; close all;

%% User Settings
max_curves_to_plot = inf; % Set this to desired number, or inf for all curves
% inf for all curves

%% Load saved data
load('fitted_SIR_curves_MORTALITY.mat');
load('fitted_parameters_MORTALITY.mat');
load('Id_data.mat');

%% Define true parameters
beta_true = 0.2;
gamma_true = 0.1;
phi_d_true = 1/10;
N = 1;
I0 = 0.05*N;
I0 = 5/10000*N;  % This becomes 0.0005
R0 = 0;
S0 = N - (I0+R0);

%% Generate true curves
y0_true = [S0; I0; R0; 0];
[t, y_true] = ode45(@(t, y) model_ode(t, y, beta_true, gamma_true, phi_d_true, N), tspan, y0_true);

%% Get data dimensions and apply limit
count_total = size(accepted_params, 1);
curves_available = size(SIR_curves, 3);

% Use minimum of available data
count = min([count_total, curves_available]);

% Apply user limit
if ~isinf(max_curves_to_plot)
    count = min(count, max_curves_to_plot);
end

fprintf('Total solutions available: %d\n', count_total);
fprintf('Total curves available: %d\n', curves_available);
fprintf('Plotting %d solutions\n', count);

%% Subset the data
if count > 0
    accepted_params_subset = accepted_params(1:count, :);
    SIR_curves_subset = SIR_curves(:, :, 1:count);
end

%% Calculate R_naught values
R_naught_true = (beta_true / gamma_true) * (S0 / N);
if count > 0
    R_naught_fitted = squeeze(SIR_curves_subset(1,6,:));
    % Ensure R_naught_fitted is the right size
    if size(R_naught_fitted, 1) ~= count
        R_naught_fitted = R_naught_fitted';
    end
end

%% Create extensive (3,4) plot - Figure 1
figure('Position', [100, 100, 1400, 1000]);
tiledlayout(3,4,"TileSpacing","compact");


% Tile 1: Daily Deaths
nexttile;
if count > 0
    h1 = plot(tspan, squeeze(SIR_curves_subset(:,1,:)/N), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.4); 
    hold on;
end
h2 = plot(t, Id_true/N, 'r--', 'LineWidth', 3);
title('$(a)\ \mathrm{Daily\ Deaths}$','FontSize', 14, 'Interpreter', 'latex');
if count > 0
    legend([h2, h1(1)], {'Observation', 'Fitted'}, 'Location', 'best','FontSize', 11);
else
    legend(h2, {'Observation'}, 'Location', 'best');
end
grid on;

% Tile 2: R_naught
nexttile;
if count > 0
    h1 = plot(1:count, R_naught_fitted, 'o', 'Color', [0.7 0.7 0.7], 'MarkerSize', 6); 
    hold on;
    h2 = plot(1:count, R_naught_true * ones(1, count), 'g-', 'LineWidth', 3);
    legend([h2, h1], {'True', 'Estimated'}, 'Location', 'best','FontSize', 11);
else
    h2 = plot(1, R_naught_true, 'g-', 'LineWidth', 3);
    legend(h2, {'True'}, 'Location', 'best');
end
title('$(b) \frac{\beta}{\widetilde{\gamma}}\,\frac{S(0)}{N(0)}$','FontSize', 14, 'Interpreter', 'latex');
xlabel('Accepted solution number');
grid on;
ylim([R_naught_true*0.999999, max(R_naught_fitted)*1.000001])




% Tile 3: phi_d vs (N - R(0)) with theoretical line
nexttile;
if count > 0
    h1 = scatter(N - accepted_params_subset(:, 5), accepted_params_subset(:, 2), 50, [0.5, 0.5, 0.5], 'filled', 'MarkerFaceAlpha', 0.6);
    hold on;
end
NminusR0_range = linspace(min(N - accepted_params_subset(:,5))*0.5, max(N - accepted_params_subset(:,5))*1.5, 100);
phi_d_theoretical = (phi_d_true * (N - R0)) ./ NminusR0_range;
h3 = plot(NminusR0_range, phi_d_theoretical, 'b-', 'LineWidth', 2);
h2 = plot(N - R0, phi_d_true, 'g*', 'MarkerSize', 12, 'LineWidth', 2);
title('$(c)\ \phi_d$ vs $(N(0) - R(0))$','FontSize', 14, 'Interpreter', 'latex');
xlabel('$N - R(0)$', 'Interpreter', 'latex');
ylabel('$\phi_d$', 'Interpreter', 'latex');
if count > 0
    legend([h2, h1, h3], {'True', 'Estimated', 'Theoretical'}, 'Location', 'best','FontSize', 11);
else
    legend([h2, h3], {'True', 'Theoretical'}, 'Location', 'best','FontSize', 11);
end
grid on;
xlim([0 1])
ylim([0 0.5])

% Skip Tiles 1 and 2 (leave blank)
nexttile; axis off;

% % Tile 3: beta vs 1/S(0)
% nexttile;
% if count > 0
%     h1 = scatter(1./accepted_params_subset(:, 3), accepted_params_subset(:, 1), 50, [0.5, 0.5, 0.5], 'filled', 'MarkerFaceAlpha', 0.6); 
%     hold on;
% end
% h2 = plot(1/S0, beta_true, 'g*', 'MarkerSize', 12, 'LineWidth', 2);
% title('$(c)\ \beta$ vs $1/S(0)$','FontSize', 14, 'Interpreter', 'latex');
% xlabel('$1/S(0)$', 'Interpreter', 'latex');
% ylabel('$\beta$', 'Interpreter', 'latex');
% if count > 0
%     legend([h2, h1], {'True', 'Estimated'}, 'Location', 'best');
% else
%     legend(h2, {'True'}, 'Location', 'best');
% end
% grid on;

% % Tile 4: beta over gamma
% nexttile;
% if count > 0
%     r_naught_fitted_mod = (accepted_params_subset(:,1) ./ gamma_true);
%     h1 = plot(1:count, r_naught_fitted_mod, 'o', 'Color', [0.7 0.7 0.7], 'MarkerSize', 6); 
%     hold on;
%     h2 = plot(1:count, R_naught_true * ones(1, count), 'g-', 'LineWidth', 3);
%     legend([h2, h1], {'True', 'Estimated'}, 'Location', 'best');
% else
%     h2 = plot(1, R_naught_true, 'g-', 'LineWidth', 3);
%     legend(h2, {'True'}, 'Location', 'best');
% end
% title('$(d) \frac{\beta}{\widetilde{\gamma}}$','FontSize', 14, 'Interpreter', 'latex');
% xlabel('Accepted solution number');
% grid on;

% Tile 5: phi_d vs Beta with theoretical line
nexttile;
h1 = scatter(accepted_params_subset(:, 1), accepted_params_subset(:, 2), 50, [0.5, 0.5, 0.5], 'filled', 'MarkerFaceAlpha', 0.6);
hold on;
beta_range = linspace(min(accepted_params_subset(:,1))*0.5, max(accepted_params_subset(:,1))*2, 100);
phi_d_theoretical = (phi_d_true/beta_true) * beta_range;
h3 = plot(beta_range, phi_d_theoretical, 'b-', 'LineWidth', 2);
h2 = plot(beta_true, phi_d_true, 'g*', 'MarkerSize', 12, 'LineWidth', 2);
title('$(d)\ \phi_d$ vs $\beta$','FontSize', 14, 'Interpreter', 'latex');
xlabel('$\beta$', 'Interpreter', 'latex');
ylabel('$\phi_d$', 'Interpreter', 'latex');
legend([h2, h1, h3], {'True', 'Estimated', 'Theoretical'}, 'Location', 'best','FontSize', 11)
grid on;
xlim([0.2 1])
ylim([0 0.5])

% Tile 6: phi_d vs S(0) with theoretical line
nexttile;
if count > 0
    h1 = scatter(accepted_params_subset(:, 3), accepted_params_subset(:, 2), 50, [0.5, 0.5, 0.5], 'filled', 'MarkerFaceAlpha', 0.6);
    hold on;
end
S0_range = linspace(min(accepted_params_subset(:,3))*0.5, max(accepted_params_subset(:,3))*2, 100);
phi_d_theoretical = (phi_d_true * S0) ./ S0_range;
h3 = plot(S0_range, phi_d_theoretical, 'b-', 'LineWidth', 2);
h2 = plot(S0, phi_d_true, 'g*', 'MarkerSize', 12, 'LineWidth', 2);
title('$(e)\ \phi_d$ vs $S(0)$','FontSize', 14, 'Interpreter', 'latex');
xlabel('$S(0)$', 'Interpreter', 'latex');
ylabel('$\phi_d$', 'Interpreter', 'latex');
if count > 0
    legend([h2, h1, h3], {'True', 'Estimated', 'Theoretical'}, 'Location', 'best','FontSize', 11);
else
    legend([h2, h3], {'True', 'Theoretical'}, 'Location', 'best','FontSize', 11);
end
grid on;
xlim([0 1])
ylim([0 0.5])



% Tile 7: phi_d vs I(0) with theoretical line
nexttile;
if count > 0
    h1 = scatter(accepted_params_subset(:, 4), accepted_params_subset(:, 2), 50, [0.5, 0.5, 0.5], 'filled', 'MarkerFaceAlpha', 0.6);
    hold on;
end
I0_range = linspace(min(accepted_params_subset(:,4))*0.5, max(accepted_params_subset(:,4))*2, 100);
phi_d_theoretical = (phi_d_true * I0) ./ I0_range;
h3 = plot(I0_range, phi_d_theoretical, 'b-', 'LineWidth', 2);
h2 = plot(I0, phi_d_true, 'g*', 'MarkerSize', 12, 'LineWidth', 2);
title('$(f)\ \phi_d$ vs $I(0)$','FontSize', 14, 'Interpreter', 'latex');
xlabel('$I(0)$', 'Interpreter', 'latex');
ylabel('$\phi_d$', 'Interpreter', 'latex');
if count > 0
    legend([h2, h1, h3], {'True', 'Estimated', 'Theoretical'}, 'Location', 'best','FontSize', 11);
else
    legend([h2, h3], {'True', 'Theoretical'}, 'Location', 'best','FontSize', 11);
end
grid on;
xlim([0 I0*1.01])
ylim([0 0.5])

% Tile 8: phi_d vs R(0)
nexttile;
if count > 0
    h1 = scatter(accepted_params_subset(:, 5), accepted_params_subset(:, 2), 50, [0.5, 0.5, 0.5], 'filled', 'MarkerFaceAlpha', 0.6);
    hold on;
end
h2 = plot(R0, phi_d_true, 'g*', 'MarkerSize', 12, 'LineWidth', 2);
title('$(g)\ \phi_d$ vs $R(0)$','FontSize', 14, 'Interpreter', 'latex');
xlabel('$R(0)$', 'Interpreter', 'latex');
ylabel('$\phi_d$', 'Interpreter', 'latex');
if count > 0
    legend([h2, h1], {'True', 'Estimated'}, 'Location', 'best','FontSize', 11);
else
    legend(h2, {'True'}, 'Location', 'best');
end
grid on;
xlim([0 1])
ylim([0 0.5])

% Tile 9: Susceptible
nexttile;
if count > 0
    h1 = plot(tspan, squeeze(SIR_curves_subset(:,2,:)/N), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.4); 
    hold on;
end
h2 = plot(t, y_true(:,1)/N, 'g--', 'LineWidth', 3);
title('$(h)\ \mathrm{Susceptible\ (S)}$','FontSize', 14, 'Interpreter', 'latex');
if count > 0
    legend([h2, h1(1)], {'True', 'Estimated'}, 'Location', 'best','FontSize', 11);
else
    legend(h2, {'True'}, 'Location', 'best');
end
grid on;

% Tile 10: Infectious
nexttile;
if count > 0
    h1 = plot(tspan, squeeze(SIR_curves_subset(:,3,:)/N), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.4); 
    hold on;
end
h2 = plot(t, y_true(:,2)/N, 'g--', 'LineWidth', 3);
title('$(i)\ \mathrm{Infectious\ (I)}$','FontSize', 14, 'Interpreter', 'latex');
if count > 0
    legend([h2, h1(1)], {'True', 'Estimated'}, 'Location', 'best','FontSize', 11);
else
    legend(h2, {'True'}, 'Location', 'best');
end
grid on;

% Tile 11: Recovered
nexttile;
if count > 0
    h1 = plot(tspan, squeeze(SIR_curves_subset(:,4,:)/N), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.4); 
    hold on;
end
h2 = plot(t, y_true(:,3)/N, 'g--', 'LineWidth', 3);
title('$(j)\ \mathrm{Recovered\ (R)}$','FontSize', 14, 'Interpreter', 'latex');
if count > 0
    legend([h2, h1(1)], {'True', 'Estimated'}, 'Location', 'best','FontSize', 11);
else
    legend(h2, {'True'}, 'Location', 'best');
end
grid on;

% Tile 12: Cumulative Deaths
nexttile;
if count > 0
    h1 = plot(tspan, squeeze(SIR_curves_subset(:,5,:)/N), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.4); 
    hold on;
end
h2 = plot(t, y_true(:,4)/N, 'g--', 'LineWidth', 3);
title('$(k)\ \mathrm{Cumulative\ Deaths\ (D)}$','FontSize', 14, 'Interpreter', 'latex');
if count > 0
    legend([h2, h1(1)], {'True', 'Estimated'}, 'Location', 'best','FontSize', 11);
else
    legend(h2, {'True'}, 'Location', 'best');
end
grid on;

% %% Parameter distributions plot - Figure 2
% if count > 0
%     figure('Position', [200, 200, 1000, 1000]);
%     n_params = 5;
%     param_names = {'$\beta$', '$\phi_d$', '$S(0)$', '$I(0)$', '$R(0)$'};
%     true_values = [beta_true, phi_d_true, S0, I0, R0];
%     tiledlayout(n_params, n_params, 'TileSpacing', 'compact', 'Padding', 'compact');
% 
%     for i = 1:n_params
%         for j = 1:n_params
%             nexttile;
%             if i > j
%                 scatter(accepted_params_subset(:, j), accepted_params_subset(:, i), '.', 'MarkerEdgeAlpha', 0.5);
%                 hold on;
%                 plot(true_values(j), true_values(i), 'g*', 'MarkerSize', 10, 'LineWidth', 2);
%                 hold off;
%                 xlabel(param_names{j}, 'Interpreter', 'latex');
%                 ylabel(param_names{i}, 'Interpreter', 'latex');
%                 title([param_names{i}, ' vs ', param_names{j}], 'Interpreter', 'latex');
%             else
%                 axis off;
%             end
%         end
%     end
% else
%     fprintf('No parameter data available for distribution plots.\n');
% end

%% Performance and completion info
if count < count_total
    fprintf('Performance: Reduced plotting time by limiting to %d/%d curves\n', count, count_total);
end
fprintf('Plotted %d accepted solutions\n', count);

%% Model ODE Function
function dydt = model_ode(~, y, beta, gamma, phi_d, N)
    S = y(1); I = y(2); R = y(3); D = y(4);
    
    dS = -beta * S * I / N;
    dI = beta * S * I / N - gamma * I;
    dR = (1 - phi_d) * gamma * I;
    dD = phi_d * gamma * I;
    
    dydt = [dS; dI; dR; dD];
end

% Uncomment to save the figure
% print('SIRmortalityv6high.png', '-dpng', '-r300');