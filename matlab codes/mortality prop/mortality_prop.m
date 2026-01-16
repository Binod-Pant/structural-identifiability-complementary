clc; clear; close all;

%% True Parameters
beta_true = 0.2;
gamma_true = 0.1;
phi_d_true = 1/10;  % Mortality fraction
N = 1;  % Small population

%% True Initial Conditions
I0 = 5/10000*N;  % This becomes 0.0005
R0 = 0;
D0 = 0;  % Initial deaths
S0 = N - (I0+R0);
y0_true = [S0; I0; R0; D0];

%% Time Span for Simulation
tspan = 0:1:200;

%% Solve ODE for True Curve with better tolerances
options_ode = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[t, y_true] = ode45(@(t, y) model_ode(t, y, beta_true, gamma_true, phi_d_true, N), tspan, y0_true, options_ode);

% Extract true daily mortality from cumulative deaths
D_cumulative_true = y_true(:, 4);
Id_true = [D_cumulative_true(1); diff(D_cumulative_true)];

%% Save Generated Data
save('Id_data.mat', 'Id_true');

%% Data Fitting and Plotting

% Plot true state variables
figure;
tile1 = tiledlayout(2,3,"TileSpacing","compact");

nexttile;
plot(t, Id_true, 'r', 'LineWidth', 3); hold on;
title('$(a)\ \mathrm{Daily\ Deaths}$','FontSize', 14, 'Interpreter', 'latex');
grid on;

nexttile;
plot(t, y_true(:,1), 'k', 'LineWidth', 3); hold on;
title('$(b)\ \mathrm{Susceptible\ (S)}$','FontSize', 14, 'Interpreter', 'latex');
grid on;

nexttile;
plot(t, y_true(:,2), 'k', 'LineWidth', 3); hold on;
title('$(c)\ \mathrm{Infectious\ (I)}$','FontSize', 14, 'Interpreter', 'latex');
grid on;

nexttile;
plot(t, y_true(:,3), 'k', 'LineWidth', 3); hold on;
title('$(d)\ \mathrm{Recovered\ (R)}$','FontSize', 14, 'Interpreter', 'latex');
grid on;

nexttile;
plot(t, y_true(:,4), 'k', 'LineWidth', 3); hold on;
title('$(e)\ \mathrm{Cumulative\ Deaths\ (D)}$','FontSize', 14, 'Interpreter', 'latex');
grid on;

nexttile;
R_naught_true = (beta_true / gamma_true) * (S0 / N);
plot(t, R_naught_true * ones(size(t)), 'k', 'LineWidth', 3); hold on;
title('$(f)\ R_0$','FontSize', 14, 'Interpreter', 'latex');
grid on;

ylabel(tile1, '$\mathrm{Number\ of\ individuals}$', 'FontSize', 18, 'Interpreter', 'latex');
xlabel(tile1, '$\mathrm{Time\ (in\ days)}$', 'FontSize', 18, 'Interpreter', 'latex');

count = 0;
total_iterations = 10000;  % Increased for N=1
SIR_curves = [];
accepted_params = [];

% Use relative error threshold for N=1
max_Id = max(Id_true);
if max_Id > 0
    relative_threshold = 0.05;  % 5% relative error
    sse_threshold = relative_threshold^2 * sum(Id_true.^2);
else
    sse_threshold = 1e-12;  % Very small absolute threshold if no signal
end

fprintf('Using SSE threshold: %.2e\n', sse_threshold);
fprintf('Max daily deaths: %.6f\n', max_Id);
count = 0;

for i = 1:total_iterations
    count = count + 1
    if mod(i, 100) == 0
        fprintf('Iteration %d/%d, accepted: %d\n', i, total_iterations, count);
    end
    
    % Better initial guesses for N=1 - focus around true values with noise
    beta_guess = rand;  % 0.5x to 1.5x true value
    phi_d_guess = rand;
    I0_guess = I0 * (0.1 + 1.8*rand);  % 0.1x to 1.9x true I0
    R0_guess = rand;  % Small variation around true R0
    
    x0 = [beta_guess, phi_d_guess, I0_guess, R0_guess];
    
    % Tighter bounds appropriate for N=1
    lb = [0, 0.0, 1e-6, 0];      % More restrictive lower bounds
    ub = [1, 1, 0.1, 1];            % More restrictive upper bounds
    
    % Define fitting function with better ODE options
    fit_fun = @(x, t) solve_model(x(1), x(2), gamma_true, N, [N - x(3) - x(4); x(3); x(4); 0], t);
    
    % Perform fitting with better options
    options = optimoptions('lsqcurvefit', 'Display', 'off', ...
        'FunctionTolerance', 1e-12, 'OptimalityTolerance', 1e-12, ...
        'StepTolerance', 1e-12);
    
    try
        x_fit = lsqcurvefit(fit_fun, x0, tspan, Id_true, lb, ub, options);
        
        % Extract fitted parameters
        beta_fit = x_fit(1);
        phi_d_fit = x_fit(2);
        I0_fit = x_fit(3);
        R0_fit = x_fit(4);
        
        % Solve ODE with fitted parameters
        S0_fit = N - I0_fit - R0_fit;
        
        % Check if initial conditions are physically reasonable
        if S0_fit < 0 || I0_fit < 0 || R0_fit < 0 || S0_fit > N
            continue;
        end
        
        y0_fit = [S0_fit; I0_fit; R0_fit; 0];
        [t_fit, y_fit] = ode45(@(t, y) model_ode(t, y, beta_fit, gamma_true, phi_d_fit, N), ...
            tspan, y0_fit, options_ode);
        
        % Compute daily mortality from fitted model
        D_fit = y_fit(:, 4);
        Id_fit = [D_fit(1); diff(D_fit)];
        
        % Compute sum squared error
        sse = sum((Id_fit - Id_true).^2);
        
        % Also check relative error for additional quality control
        if sum(Id_true.^2) > 0
            relative_error = sse / sum(Id_true.^2);
        else
            relative_error = sse;
        end
        
        % Accept if error is below threshold AND relative error is reasonable
        if sse <= sse_threshold && relative_error <= 0.05
            count = count + 1;
            accepted_params = [accepted_params; beta_fit, phi_d_fit, S0_fit, I0_fit, R0_fit];
            
            % Calculate R_naught for this parameter set
            R_naught_fit = (beta_fit / gamma_true) * (S0_fit / N);
            
            % Store curves: (1) fit (2) S (3) I (4) R (5) D (6) R_naught
            temp_data = zeros(length(tspan), 6);
            temp_data(:,1) = Id_fit;
            temp_data(:,2) = y_fit(:,1);
            temp_data(:,3) = y_fit(:,2);
            temp_data(:,4) = y_fit(:,3);
            temp_data(:,5) = y_fit(:,4);
            temp_data(:,6) = R_naught_fit;
            
            SIR_curves(:,:,count) = temp_data;

            % Plot fitted results
            nexttile(1); 
            plot(tspan, Id_fit, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.4);
            
            nexttile(2); 
            plot(tspan, y_fit(:,1), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.4);
            
            nexttile(3); 
            plot(tspan, y_fit(:,2), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.4);
            
            nexttile(4); 
            plot(tspan, y_fit(:,3), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.4);
            
            nexttile(5);
            plot(tspan, y_fit(:,4), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.4);
            
            nexttile(6);
            plot(tspan, R_naught_fit * ones(size(tspan)), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.4);
        end
        
    catch ME
        % Skip this iteration if optimization fails
        continue;
    end
end

% After the loop ends, re-plot the true curves on top
nexttile(1);
plot(t, Id_true, 'r', 'LineWidth', 3);
nexttile(2);
plot(t, y_true(:,1), 'k', 'LineWidth', 3);
nexttile(3);
plot(t, y_true(:,2), 'k', 'LineWidth', 3);
nexttile(4);
plot(t, y_true(:,3), 'k', 'LineWidth', 3);
nexttile(5);
plot(t, y_true(:,4), 'k', 'LineWidth', 3);
nexttile(6);
plot(t, R_naught_true * ones(size(t)), 'k', 'LineWidth', 3);

success_rate = count/total_iterations;
fprintf('Success rate: %.4f (%d accepted out of %d)\n', success_rate, count, total_iterations);

if count > 0
    fprintf('\nParameter Estimation Results:\n');
    fprintf('β: %.4f ± %.4f (true: %.4f)\n', mean(accepted_params(:,1)), std(accepted_params(:,1)), beta_true);
    fprintf('φ_d: %.4f ± %.4f (true: %.4f)\n', mean(accepted_params(:,2)), std(accepted_params(:,2)), phi_d_true);
    fprintf('S(0): %.6f ± %.6f (true: %.6f)\n', mean(accepted_params(:,3)), std(accepted_params(:,3)), S0);
    fprintf('I(0): %.6f ± %.6f (true: %.6f)\n', mean(accepted_params(:,4)), std(accepted_params(:,4)), I0);
    fprintf('R(0): %.6f ± %.6f (true: %.6f)\n', mean(accepted_params(:,5)), std(accepted_params(:,5)), R0);
end

% Save fitted parameters and curves
save('fitted_SIR_curves_MORTALITY.mat', 'SIR_curves', 'tspan');
save('fitted_parameters_MORTALITY.mat', 'accepted_params');

%% Plot Parameter Distributions and Joint Distributions
if count > 0
    figure;
    n_params = 5;
    param_names = {'$\beta$', '$\phi_d$', '$S(0)$', '$I(0)$', '$R(0)$'};
    true_values = [beta_true, phi_d_true, S0, I0, R0];
    tiledlayout(n_params, n_params, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    for i = 1:n_params
        for j = 1:n_params
            nexttile;
            if i > j
                scatter(accepted_params(:, j), accepted_params(:, i), '.', 'MarkerEdgeAlpha', 0.5);
                hold on;
                plot(true_values(j), true_values(i), 'r*', 'MarkerSize', 10, 'LineWidth', 2);
                hold off;
                xlabel(param_names{j}, 'Interpreter', 'latex');
                ylabel(param_names{i}, 'Interpreter', 'latex');
                title([param_names{i}, ' vs ', param_names{j}], 'Interpreter', 'latex');
            else
                axis off;
            end
        end
    end
    
    fig2 = gcf;
    set(fig2, 'PaperUnits', 'inches', 'PaperPosition', [0 0 12 12]);
    set(fig2, 'Position', [100, 100, 600, 600]);
    print(fig2, 'parameter_distributions', '-dpng', '-r300');
else
    fprintf('No accepted parameter sets - cannot create parameter distribution plot\n');
end

%% Model ODE Function
function dydt = model_ode(~, y, beta, gamma, phi_d, N)
    S = y(1); I = y(2); R = y(3); D = y(4);
    
    dS = -beta * S * I / N;
    dI = beta * S * I / N - gamma * I;
    dR = (1 - phi_d) * gamma * I;
    dD = phi_d * gamma * I;
    
    dydt = [dS; dI; dR; dD];
end

%% Function to Solve Model for Fitting
function Id_sim = solve_model(beta, phi_d, gamma, N, y0, tspan)
    options_ode = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    [~, y] = ode45(@(t, y) model_ode(t, y, beta, gamma, phi_d, N), tspan, y0, options_ode);
    D_sim = y(:, 4);
    Id_sim = [D_sim(1); diff(D_sim)];
end