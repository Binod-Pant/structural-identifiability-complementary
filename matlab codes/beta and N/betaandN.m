% SIR Model Parameter Estimation with Accept/Reject
clear; clc; close all;

%% True Parameters
beta_true = 0.3;
gamma_true = 0.1;
N_true = 1;
I0_true = 1/100 * N_true;
R0_true = 1/100 * N_true;
S0_true = N_true - (I0_true + R0_true);

% Define R_naught true
R_naught_true = beta_true * S0_true / N_true * (1/gamma_true);

%% Generate True Data
tspan = 0:1:100;  % Time span
y0_true = [S0_true, I0_true, R0_true];

% Solve ODE system
[t, y_true] = ode45(@(t,y) sir_model(t, y, beta_true, gamma_true, N_true), tspan, y0_true);
S_true = y_true(:,1);
I_true = y_true(:,2);
R_true = y_true(:,3);

% Generate incidence data (beta*S*I/N)
incidence_true = beta_true * S_true .* I_true / N_true;

%% Set SSE Threshold
% Use 1% of the total sum of squares of the true incidence as threshold
total_ss = sum(incidence_true.^2);
sse_threshold = 0.00001 * total_ss;  % Accept if SSE < 0.0011% of total sum of squares
fprintf('SSE threshold: %.6f\n', sse_threshold);

%% Parameter Estimation Setup
max_iterations = 1500;

% Storage for ACCEPTED results only
accepted_params = [];
accepted_S = {};
accepted_I = {};
accepted_R = {};
accepted_incidence = {};
accepted_errors = [];
accepted_R_naught = [];  % Add R_naught storage
accepted_count = 0;

fprintf('Starting parameter estimation...\n');

%% Iterative Parameter Estimation with Accept/Reject
for iter = 1:max_iterations
    fprintf('Iteration %d: ', iter);
    
    % Random initial guess
    beta_guess = rand();
    S0_guess = rand();
    I0_guess = rand();
    R0_guess = rand();
    
    params0 = [beta_guess, S0_guess, I0_guess, R0_guess];
    
    % Optimization
    options = optimoptions('fminunc', 'Display', 'off', 'MaxIterations', 1000);
    
    try
        [params_opt, ~] = fminunc(@(params) objective_function(params, t, incidence_true, gamma_true), params0, options);
        
        % Extract parameters
        beta_est = params_opt(1);
        S0_est = params_opt(2);
        I0_est = params_opt(3);
        R0_est = params_opt(4);
        N_est = S0_est + I0_est + R0_est;
        
        % Check if parameters are positive
        if beta_est > 0 && S0_est > 0 && I0_est > 0 && R0_est > 0
            % Generate fitted curves
            y0_fit = [S0_est, I0_est, R0_est];
            [~, y_fit] = ode45(@(t,y) sir_model(t, y, beta_est, gamma_true, N_est), tspan, y0_fit);
            S_fit = y_fit(:,1);
            I_fit = y_fit(:,2);
            R_fit = y_fit(:,3);
            incidence_fit = beta_est * S_fit .* I_fit / N_est;
            
            % Calculate SSE between true and fitted incidence
            sse = sum((incidence_true - incidence_fit).^2);
            
            % Calculate R_naught estimate
            R_naught_est = beta_est * S0_est / N_est * (1/gamma_true);
            
            % Accept/Reject based on SSE threshold
            %if sse < sse_threshold
            if sse < sse_threshold && R0_est < 1.5
                accepted_count = accepted_count + 1;
                accepted_params(accepted_count, :) = params_opt;
                accepted_S{accepted_count} = S_fit;
                accepted_I{accepted_count} = I_fit;
                accepted_R{accepted_count} = R_fit;
                accepted_incidence{accepted_count} = incidence_fit;
                accepted_errors(accepted_count) = sse;
                accepted_R_naught(accepted_count) = R_naught_est;  % Save R_naught
                
                fprintf('ACCEPTED - SSE = %.6f, beta = %.4f, N = %.4f\n', sse, beta_est, N_est);
            else
                fprintf('REJECTED - SSE = %.6f (threshold = %.6f)\n', sse, sse_threshold);
            end
        else
            fprintf('REJECTED - Negative parameters\n');
        end
        
    catch
        fprintf('REJECTED - Optimization failed\n');
    end
end

fprintf('\nTotal accepted: %d out of %d iterations\n', accepted_count, max_iterations);

%% Plotting
figure('Position', [100, 100, 1200, 1200]);

% First row: Plot all accepted gray curves first, then true curves on top
subplot(2,4,1)
for i = 1:accepted_count
    plot(t, accepted_incidence{i}, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.7); hold on;
end
plot(t, incidence_true, 'r-', 'LineWidth', 2);
xlabel('Time'); ylabel('Incidence');
title('$\beta \times S \times I/N$', 'Interpreter', 'latex');
grid on;

% Third row: Initial conditions S(0), I(0), R(0)
subplot(3,4,9)
if accepted_count > 0
    S0_estimates = accepted_params(:,2);
    plot(1:accepted_count, S0_estimates, 'bo-', 'LineWidth', 2); hold on;
    plot([1, accepted_count], [S0_true, S0_true], 'r-', 'LineWidth', 2);
    xlabel('Accepted Solution #'); ylabel('S(0) estimate');
    title('Initial Susceptible S(0)', 'Interpreter', 'latex');
    legend('Estimated', 'True', 'Location', 'best');
end
grid on;

subplot(3,4,10)
if accepted_count > 0
    I0_estimates = accepted_params(:,3);
    plot(1:accepted_count, I0_estimates, 'bo-', 'LineWidth', 2); hold on;
    plot([1, accepted_count], [I0_true, I0_true], 'r-', 'LineWidth', 2);
    xlabel('Accepted Solution #'); ylabel('I(0) estimate');
    title('Initial Infected I(0)', 'Interpreter', 'latex');
    legend('Estimated', 'True', 'Location', 'best');
end
grid on;

subplot(3,4,11)
if accepted_count > 0
    R0_estimates = accepted_params(:,4);
    plot(1:accepted_count, R0_estimates, 'bo-', 'LineWidth', 2); hold on;
    plot([1, accepted_count], [R0_true, R0_true], 'r-', 'LineWidth', 2);
    xlabel('Accepted Solution #'); ylabel('R(0) estimate');
    title('Initial Recovered R(0)', 'Interpreter', 'latex');
    legend('Estimated', 'True', 'Location', 'best');
end
grid on;

subplot(3,4,12)
% Leave empty or add another plot if needed

subplot(3,4,2)
for i = 1:accepted_count
    plot(t, accepted_S{i}, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.7); hold on;
end
plot(t, S_true, 'k-', 'LineWidth', 2);
xlabel('Time'); ylabel('S(t)');
title('Susceptible (S)', 'Interpreter', 'latex');
grid on;

subplot(3,4,3)
for i = 1:accepted_count
    plot(t, accepted_I{i}, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.7); hold on;
end
plot(t, I_true, '--k', 'LineWidth', 2);
xlabel('Time'); ylabel('I(t)');
title('Infected (I)', 'Interpreter', 'latex');
grid on;

subplot(3,4,4)
for i = 1:accepted_count
    plot(t, accepted_R{i}, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.7); hold on;
end
plot(t, R_true, 'k-', 'LineWidth', 2);
xlabel('Time'); ylabel('R(t)');
title('Recovered (R)', 'Interpreter', 'latex');
grid on;

% Second row: Parameter estimation results (only accepted)
subplot(2,4,5)
if accepted_count > 0
    beta_estimates = accepted_params(:,1);
    plot(1:accepted_count, beta_estimates, 'bo-', 'LineWidth', 2); hold on;
    plot([1, accepted_count], [beta_true, beta_true], 'r-', 'LineWidth', 2);
    xlabel('Accepted Solution #'); ylabel('$\beta$ estimate', 'Interpreter', 'latex');
    title('Beta Estimation (Accepted)');
    legend('Estimated', 'True', 'Location', 'best');
end
grid on;

subplot(3,4,6)
if accepted_count > 0
    N_estimates = sum(accepted_params(:,2:4), 2); % S0 + I0 + R0
    plot(1:accepted_count, N_estimates, 'bo-', 'LineWidth', 2); hold on;
    plot([1, accepted_count], [N_true, N_true], 'r-', 'LineWidth', 2);
    xlabel('Accepted Solution #'); ylabel('N estimate');
    title('Population Size (Accepted)');
    legend('Estimated', 'True', 'Location', 'best');
end
grid on;

subplot(3,4,7)
if accepted_count > 0
    beta_estimates = accepted_params(:,1);
    N_estimates = sum(accepted_params(:,2:4), 2);
    plot(N_estimates, beta_estimates, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'blue'); hold on;
    plot(N_true, beta_true, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'red');
    xlabel('N estimate'); ylabel('$\beta$ estimate', 'Interpreter', 'latex');
    title('$\beta$ vs N (Accepted)', 'Interpreter', 'latex');
    legend('Estimated', 'True', 'Location', 'best');
end
grid on;

subplot(3,4,8)
if accepted_count > 0
    semilogy(1:accepted_count, accepted_errors, 'bo-', 'LineWidth', 2);
    hold on;
    plot([1, accepted_count], [sse_threshold, sse_threshold], 'r--', 'LineWidth', 2);
    xlabel('Accepted Solution #'); ylabel('SSE');
    title('SSE (Accepted Solutions)');
    legend('SSE', 'Threshold', 'Location', 'best');
end
grid on;

sgtitle(sprintf('SIR Model Parameter Estimation (%d/%d Accepted)', accepted_count, max_iterations));

%% Plot R_naught
figure;
histogram(accepted_R_naught, 'FaceColor', 'green'); hold on;
xline(R_naught_true, 'r-', 'LineWidth', 2);
xline(1, 'k--', 'LineWidth', 2);
xlabel('R_0'); ylabel('Frequency');
title('R_0 Distribution');
legend('Estimates', 'True', 'Threshold');

%% Save Results to .mat file
results.true_params = struct('beta', beta_true, 'gamma', gamma_true, 'N', N_true, ...
                            'S0', S0_true, 'I0', I0_true, 'R0', R0_true);
results.time = t;
results.true_curves = struct('S', S_true, 'I', I_true, 'R', R_true, 'incidence', incidence_true);
results.sse_threshold = sse_threshold;
results.accepted_count = accepted_count;
results.total_iterations = max_iterations;
results.accepted_params = accepted_params;
results.accepted_curves = struct('S', {accepted_S}, 'I', {accepted_I}, 'R', {accepted_R}, ...
                                'incidence', {accepted_incidence});
results.accepted_errors = accepted_errors;
results.accepted_R_naught = accepted_R_naught;

% Save with timestamp
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
filename = sprintf('SIR_estimation_results.mat', timestamp);
save(filename, 'results');
fprintf('\nResults saved to: %s\n', filename);

%% Functions
function dydt = sir_model(t, y, beta, gamma, N)
    S = y(1);
    I = y(2);
    R = y(3);
    
    dSdt = -beta * S * I / N;
    dIdt = beta * S * I / N - gamma * I;
    dRdt = gamma * I;
    
    dydt = [dSdt; dIdt; dRdt];
end

function error = objective_function(params, t, incidence_data, gamma)
    beta = params(1);
    S0 = params(2);
    I0 = params(3);
    R0 = params(4);
    N = S0 + I0 + R0;
    
    % Ensure positive parameters
    if beta <= 0 || S0 <= 0 || I0 <= 0 || R0 <= 0 || N <= 0
        error = 1e10;
        return;
    end
    
    try
        % Solve ODE
        y0 = [S0, I0, R0];
        [~, y] = ode45(@(t,y) sir_model(t, y, beta, gamma, N), t, y0);
        
        % Calculate predicted incidence
        S_pred = y(:,1);
        I_pred = y(:,2);
        incidence_pred = beta * S_pred .* I_pred / N;
        
        % Calculate error (sum of squared differences)
        error = sum((incidence_data - incidence_pred).^2);
        
    catch
        error = 1e10;
    end
end