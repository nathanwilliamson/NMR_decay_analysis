function fit = analyze(t, I, model, baseline, number_of_fits, number_of_mc_fits)

% The number of components in the model (excluding the baseline).
number_of_components = numel(model);

% Optimization algorithm settings. Using optimoptions and setting two sets
% of options should make this work for Matlab R2013a and above.
fprintf('Setting optimization options...\n');

options = optimoptions('fmincon');
options.Algorithm = 'sqp';
options.Display = 'off';

options.MaxFunctionEvaluations = 10000; % Modern settings.
options.MaxIterations = 1000;
options.OptimalityTolerance = 1e-8;
options.ConstraintTolerance = 1e-8;

options.MaxFunEvals = 10000; % Legacy settings.
options.MaxIter = 1000;
options.TolFun = 1e-8;
options.TolX = 1e-8;

% Bounds and constraints for parameters.
[lb, ub, Aeq, beq] = bounds_and_constraints(model, baseline);

% Fit model.
ss = inf;
fprintf('Fitting...\n');
for current_fit = 1:number_of_fits
    fprintf('   Fit %d out of %d...\n', current_fit, number_of_fits);
    
    % Generate initial values of parameters.
    param_guess = random_initial_guess(model, lb, ub, b, I);
    
    % Fit model.
    paramhat_ = fmincon(@(param)sumofsquares(b, I, model, param), param_guess, [], [], Aeq, beq, lb, ub, [], options);
    
    % Compute residual sum of squares.
    Imodel_ = signal(b, model, paramhat_);
    ss_ = sum( (I - Imodel_).^2 );
    
    if ss_ < ss
        paramhat = paramhat_;
        ss = ss_;
        Imodel = Imodel_;
    end
end

% Error analysis.
sigma_residual = sqrt(ss / numel(b));
if number_of_mc_fits < 2
    paramhat_MC = nan(size(paramhat));
else
    paramhat_MC = zeros(number_of_mc_fits, numel(paramhat));
    fprintf('Error analysis...\n');
    for current_fit = 1:number_of_mc_fits
        fprintf('   Fit %d out of %d...\n', current_fit, number_of_mc_fits);
        I_MC = I + sigma_residual * randn(size(I));
        paramhat_MC(current_fit, :) = fmincon(@(param)sumofsquares(b, I_MC, model, param), paramhat, [], [], Aeq, beq, lb, ub, [], options);
    end
end

% Structure output.
fprintf('Create output...\n');
fit = struct();
fit.ss = ss;

fit.I0 = paramhat(end);
fit.std_I0 = std(paramhat_MC(:, end), [], 1);

theta = paramhat(end-number_of_components:end-1);
theta_MC = paramhat_MC(:, end-number_of_components:end-1);
std_theta = std(theta_MC, [], 1);

fit.components = cell(number_of_components, 1);
for current_component = 1:number_of_components
    fit.components{current_component} = struct();
end

ind = 1;
for current_component = 1:number_of_components
    
    fit.components{current_component}.model = model{current_component}{1};
    
    switch model{current_component}{1}
        case 'exponential'
            D = paramhat(ind) / tmax;
            D_MC = paramhat_MC(:, ind) / tmax;
            
            fit.components{current_component}.D = D;
            fit.components{current_component}.std_D = std(D_MC, [], 1);
            fit.components{current_component}.D_MC = D_MC;
            
            fit.components{current_component}.meanD = D;
            fit.components{current_component}.std_meanD = std(D_MC, [], 1);
            fit.components{current_component}.meanD_MC = D_MC;
            
            fit.components{current_component}.stdD = 0;
            fit.components{current_component}.std_stdD = 0;
            fit.components{current_component}.stdD_MC = zeros(number_of_mc_fits, 1);
            
            fit.components{current_component}.spreadD = 0;
            fit.components{current_component}.std_spreadD = 0;
            fit.components{current_component}.spreadD_MC = zeros(number_of_mc_fits, 1);
            
            fit.components{current_component}.modeD = D;
            fit.components{current_component}.std_modeD = std(D_MC, [], 1);
            fit.components{current_component}.modeD_MC = D_MC;
            
            fit.components{current_component}.theta = theta(current_component);
            fit.components{current_component}.std_theta = std_theta(current_component);
            fit.components{current_component}.theta_MC = theta_MC(:, current_component);
            
            ind = ind + 1;
        case 'stretchedexponential'
            D = paramhat(ind) / tmax;
            D_MC = paramhat_MC(:, ind) / tmax;
            
            beta = paramhat(ind + 1);
            beta_MC = paramhat_MC(:, ind + 1);
            
            fit.components{current_component}.D = D;
            fit.components{current_component}.std_D = std(D_MC, [], 1);
            fit.components{current_component}.D_MC = D_MC;
            
            fit.components{current_component}.beta = beta;
            fit.components{current_component}.std_beta = std(beta_MC, [], 1);
            fit.components{current_component}.beta_MC = beta_MC;
            
            fit.components{current_component}.meanD = nan;
            fit.components{current_component}.std_meanD = nan;
            fit.components{current_component}.meanD_MC = nan(number_of_mc_fits, 1);
            
            fit.components{current_component}.stdD = nan;
            fit.components{current_component}.std_stdD = nan;
            fit.components{current_component}.stdD_MC = nan(number_of_mc_fits, 1);
            
            fit.components{current_component}.spreadD = nan;
            fit.components{current_component}.std_spreadD = nan;
            fit.components{current_component}.spreadD_MC = nan(number_of_mc_fits, 1);
            
            fit.components{current_component}.modeD = nan;
            fit.components{current_component}.std_modeD = nan;
            fit.components{current_component}.modeD_MC = nan(number_of_mc_fits, 1);
            
            fit.components{current_component}.theta = theta(current_component);
            fit.components{current_component}.std_theta = std_theta(current_component);
            fit.components{current_component}.theta_MC = theta_MC(:, current_component);
            
            ind = ind + 2;
        case 'lognormal'
            mu = paramhat(ind) - log(tmax);
            mu_MC = paramhat_MC(:, ind) - log(tmax);
            
            sigma = paramhat(ind + 1);
            sigma_MC = paramhat_MC(:, ind + 1);
            
            fit.components{current_component}.mu = mu;
            fit.components{current_component}.std_mu = std(mu_MC, [], 1);
            fit.components{current_component}.mu_MC = mu_MC;
            
            fit.components{current_component}.sigma = sigma;
            fit.components{current_component}.std_sigma = std(sigma_MC, [], 1);
            fit.components{current_component}.sigma_MC = sigma_MC;
            
            meanD = exp(mu + sigma^2 / 2);
            meanD_MC = exp(mu_MC + sigma_MC.^2 / 2);
            
            stdD = sqrt(exp(sigma^2) - 1) * exp(mu + sigma^2 / 2);
            stdD_MC = sqrt(exp(sigma_MC.^2) - 1) .* exp(mu_MC + sigma_MC.^2 / 2);
            
            spreadD = sqrt(exp(sigma^2) - 1);
            spreadD_MC = sqrt(exp(sigma_MC.^2) - 1);
            
            modeD = exp(mu - sigma^2);
            modeD_MC = exp(mu_MC - sigma_MC.^2);
            
            fit.components{current_component}.meanD = meanD;
            fit.components{current_component}.std_meanD = std(meanD_MC, [], 1);
            fit.components{current_component}.meanD_MC = meanD_MC;
            
            fit.components{current_component}.stdD = stdD;
            fit.components{current_component}.std_stdD = std(stdD_MC, [], 1);
            fit.components{current_component}.stdD_MC = stdD_MC;
            
            fit.components{current_component}.spreadD = spreadD;
            fit.components{current_component}.std_spreadD = std(spreadD_MC, [], 1);
            fit.components{current_component}.spreadD_MC = spreadD_MC;
            
            fit.components{current_component}.modeD = modeD;
            fit.components{current_component}.std_modeD = std(modeD_MC, [], 1);
            fit.components{current_component}.modeD_MC = modeD_MC;
            
            fit.components{current_component}.theta = theta(current_component);
            fit.components{current_component}.std_theta = std_theta(current_component);
            fit.components{current_component}.theta_MC = theta_MC(:, current_component);
            
            ind = ind + 2;
        case 'gamma'
            alpha = paramhat(ind);
            alpha_MC = paramhat_MC(:, ind);
            
            beta = paramhat(ind + 1) * tmax;
            beta_MC = paramhat_MC(:, ind + 1) * tmax;
            
            fit.components{current_component}.alpha = alpha;
            fit.components{current_component}.std_alpha = std(alpha_MC, [], 1);
            fit.components{current_component}.alpha_MC = alpha_MC;
            
            fit.components{current_component}.beta = beta;
            fit.components{current_component}.std_beta = std(beta_MC, [], 1);
            fit.components{current_component}.beta_MC = beta_MC;
            
            meanD = alpha / beta;
            meanD_MC = alpha_MC ./ beta_MC;
            
            stdD = sqrt(alpha) / beta;
            stdD_MC = sqrt(alpha_MC) ./ beta_MC;
            
            spreadD = 1 / sqrt(alpha);
            spreadD_MC = 1 ./ sqrt(alpha_MC);
            
            modeD = (alpha - 1) / beta;
            modeD_MC = (alpha_MC - 1) ./ beta_MC;
            
            fit.components{current_component}.meanD = meanD;
            fit.components{current_component}.std_meanD = std(meanD_MC, [], 1);
            fit.components{current_component}.meanD_MC = meanD_MC;
            
            fit.components{current_component}.stdD = stdD;
            fit.components{current_component}.std_stdD = std(stdD_MC, [], 1);
            fit.components{current_component}.stdD_MC = stdD_MC;
            
            fit.components{current_component}.spreadD = spreadD;
            fit.components{current_component}.std_spreadD = std(spreadD_MC, [], 1);
            fit.components{current_component}.spreadD_MC = spreadD_MC;
            
            fit.components{current_component}.modeD = modeD;
            fit.components{current_component}.std_modeD = std(modeD_MC, [], 1);
            fit.components{current_component}.modeD_MC = modeD_MC;
            
            fit.components{current_component}.theta = theta(current_component);
            fit.components{current_component}.std_theta = std_theta(current_component);
            fit.components{current_component}.theta_MC = theta_MC(:, current_component);
            
            ind = ind + 2;
    end
end
fit.baseline = paramhat(ind);
fit.std_baseline = std(paramhat_MC(:, ind), [], 1);
fit.baseline_MC = paramhat_MC(:, ind);

fit.Imodel = Imodel;
fit.residuals = I - Imodel;

fit.paramhat_MC = paramhat_MC;

fprintf('Done.\n');

end
