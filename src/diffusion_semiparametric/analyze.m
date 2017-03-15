function fit = analyze(b, I, model, baseline, number_of_fits, number_of_mc_fits)

% The number of components in the model (excluding the baseline).
number_of_components = numel(model);

% Rescaling for numerical stability.
bmax = max(b);
b = b / bmax;

% Optimization algorithm settings. Using optimoptions and setting two sets
% of options should make this work for Matlab R2013a and above.
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
[lb, ub, Aeq, beq] = bounds_and_constraints(bmax, model, baseline);

% Fit model.
ss = inf;

for current_fit = 1:number_of_fits
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
    
    for current_fit = 1:number_of_mc_fits
        disp(current_fit)
        I_MC = I + sigma_residual * randn(size(I));
        paramhat_MC(current_fit, :) = fmincon(@(param)sumofsquares(b, I_MC, model, param), paramhat, [], [], Aeq, beq, lb, ub, [], options);
    end
end

% Structure output.
fit = struct();
fit.ss = ss;

fit.I0 = paramhat(end);
fit.std_I0 = std(paramhat_MC(:, end), [], 1);

theta = paramhat(end-number_of_components:end-1);
theta_MC = paramhat_MC(:, end-number_of_components:end-1);
std_theta = std(theta_MC, [], 1);

fit.components = cell(number_of_components, 1);
for currentComponent = 1:number_of_components
    fit.components{currentComponent} = struct();
end

ind = 1;
for currentComponent = 1:number_of_components
    
    fit.components{currentComponent}.model = model{currentComponent}{1};
    
    switch model{currentComponent}{1}
        case 'exponential'
            D = paramhat(ind) / bmax;
            D_MC = paramhat_MC(:, ind) / bmax;
            
            fit.components{currentComponent}.D = D;
            fit.components{currentComponent}.std_D = std(D_MC, [], 1);
            fit.components{currentComponent}.D_MC = D_MC;
            
            fit.components{currentComponent}.meanD = D;
            fit.components{currentComponent}.std_meanD = std(D_MC, [], 1);
            fit.components{currentComponent}.meanD_MC = D_MC;
            
            fit.components{currentComponent}.stdD = 0;
            fit.components{currentComponent}.std_stdD = 0;
            fit.components{currentComponent}.stdD_MC = zeros(number_of_mc_fits, 1);
            
            fit.components{currentComponent}.spreadD = 0;
            fit.components{currentComponent}.std_spreadD = 0;
            fit.components{currentComponent}.spreadD_MC = zeros(number_of_mc_fits, 1);
            
            fit.components{currentComponent}.modeD = D;
            fit.components{currentComponent}.std_modeD = std(D_MC, [], 1);
            fit.components{currentComponent}.modeD_MC = D_MC;
            
            fit.components{currentComponent}.theta = theta(currentComponent);
            fit.components{currentComponent}.std_theta = std_theta(currentComponent);
            fit.components{currentComponent}.theta_MC = theta_MC(:, currentComponent);
            
            ind = ind + 1;
        case 'stretchedexponential'
            D = paramhat(ind) / bmax;
            D_MC = paramhat_MC(:, ind) / bmax;
            
            beta = paramhat(ind + 1);
            beta_MC = paramhat_MC(:, ind + 1);
            
            fit.components{currentComponent}.D = D;
            fit.components{currentComponent}.std_D = std(D_MC, [], 1);
            fit.components{currentComponent}.D_MC = D_MC;
            
            fit.components{currentComponent}.beta = beta;
            fit.components{currentComponent}.std_beta = std(beta_MC, [], 1);
            fit.components{currentComponent}.beta_MC = beta_MC;
            
            fit.components{currentComponent}.meanD = nan;
            fit.components{currentComponent}.std_meanD = nan;
            fit.components{currentComponent}.meanD_MC = nan(number_of_mc_fits, 1);
            
            fit.components{currentComponent}.stdD = nan;
            fit.components{currentComponent}.std_stdD = nan;
            fit.components{currentComponent}.stdD_MC = nan(number_of_mc_fits, 1);
            
            fit.components{currentComponent}.spreadD = nan;
            fit.components{currentComponent}.std_spreadD = nan;
            fit.components{currentComponent}.spreadD_MC = nan(number_of_mc_fits, 1);
            
            fit.components{currentComponent}.modeD = nan;
            fit.components{currentComponent}.std_modeD = nan;
            fit.components{currentComponent}.modeD_MC = nan(number_of_mc_fits, 1);
            
            fit.components{currentComponent}.theta = theta(currentComponent);
            fit.components{currentComponent}.std_theta = std_theta(currentComponent);
            fit.components{currentComponent}.theta_MC = theta_MC(:, currentComponent);
            
            ind = ind + 2;
        case 'lognormal'
            mu = paramhat(ind) - log(bmax);
            mu_MC = paramhat_MC(:, ind) - log(bmax);
            
            sigma = paramhat(ind + 1);
            sigma_MC = paramhat_MC(:, ind + 1);
            
            fit.components{currentComponent}.mu = mu;
            fit.components{currentComponent}.std_mu = std(mu_MC, [], 1);
            fit.components{currentComponent}.mu_MC = mu_MC;
            
            fit.components{currentComponent}.sigma = sigma;
            fit.components{currentComponent}.std_sigma = std(sigma_MC, [], 1);
            fit.components{currentComponent}.sigma_MC = sigma_MC;
            
            meanD = exp(mu + sigma^2 / 2);
            meanD_MC = exp(mu_MC + sigma_MC.^2 / 2);
            
            stdD = sqrt(exp(sigma^2) - 1) * exp(mu + sigma^2 / 2);
            stdD_MC = sqrt(exp(sigma_MC.^2) - 1) .* exp(mu_MC + sigma_MC.^2 / 2);
            
            spreadD = sqrt(exp(sigma^2) - 1);
            spreadD_MC = sqrt(exp(sigma_MC.^2) - 1);
            
            modeD = exp(mu - sigma^2);
            modeD_MC = exp(mu_MC - sigma_MC.^2);
            
            fit.components{currentComponent}.meanD = meanD;
            fit.components{currentComponent}.std_meanD = std(meanD_MC, [], 1);
            fit.components{currentComponent}.meanD_MC = meanD_MC;
            
            fit.components{currentComponent}.stdD = stdD;
            fit.components{currentComponent}.std_stdD = std(stdD_MC, [], 1);
            fit.components{currentComponent}.stdD_MC = stdD_MC;
            
            fit.components{currentComponent}.spreadD = spreadD;
            fit.components{currentComponent}.std_spreadD = std(spreadD_MC, [], 1);
            fit.components{currentComponent}.spreadD_MC = spreadD_MC;
            
            fit.components{currentComponent}.modeD = modeD;
            fit.components{currentComponent}.std_modeD = std(modeD_MC, [], 1);
            fit.components{currentComponent}.modeD_MC = modeD_MC;
            
            fit.components{currentComponent}.theta = theta(currentComponent);
            fit.components{currentComponent}.std_theta = std_theta(currentComponent);
            fit.components{currentComponent}.theta_MC = theta_MC(:, currentComponent);
            
            ind = ind + 2;
        case 'gamma'
            alpha = paramhat(ind);
            alpha_MC = paramhat_MC(:, ind);
            
            beta = paramhat(ind + 1) * bmax;
            beta_MC = paramhat_MC(:, ind + 1) * bmax;
            
            fit.components{currentComponent}.alpha = alpha;
            fit.components{currentComponent}.std_alpha = std(alpha_MC, [], 1);
            fit.components{currentComponent}.alpha_MC = alpha_MC;
            
            fit.components{currentComponent}.beta = beta;
            fit.components{currentComponent}.std_beta = std(beta_MC, [], 1);
            fit.components{currentComponent}.beta_MC = beta_MC;
            
            meanD = alpha / beta;
            meanD_MC = alpha_MC ./ beta_MC;
            
            stdD = sqrt(alpha) / beta;
            stdD_MC = sqrt(alpha_MC) ./ beta_MC;
            
            spreadD = 1 / sqrt(alpha);
            spreadD_MC = 1 ./ sqrt(alpha_MC);
            
            modeD = (alpha - 1) / beta;
            modeD_MC = (alpha_MC - 1) ./ beta_MC;
            
            fit.components{currentComponent}.meanD = meanD;
            fit.components{currentComponent}.std_meanD = std(meanD_MC, [], 1);
            fit.components{currentComponent}.meanD_MC = meanD_MC;
            
            fit.components{currentComponent}.stdD = stdD;
            fit.components{currentComponent}.std_stdD = std(stdD_MC, [], 1);
            fit.components{currentComponent}.stdD_MC = stdD_MC;
            
            fit.components{currentComponent}.spreadD = spreadD;
            fit.components{currentComponent}.std_spreadD = std(spreadD_MC, [], 1);
            fit.components{currentComponent}.spreadD_MC = spreadD_MC;
            
            fit.components{currentComponent}.modeD = modeD;
            fit.components{currentComponent}.std_modeD = std(modeD_MC, [], 1);
            fit.components{currentComponent}.modeD_MC = modeD_MC;
            
            fit.components{currentComponent}.theta = theta(currentComponent);
            fit.components{currentComponent}.std_theta = std_theta(currentComponent);
            fit.components{currentComponent}.theta_MC = theta_MC(:, currentComponent);
            
            ind = ind + 2;
    end
end
fit.baseline = paramhat(ind);
fit.std_baseline = std(paramhat_MC(:, ind), [], 1);
fit.baseline_MC = paramhat_MC(:, ind);

fit.Imodel = Imodel;
fit.residuals = I - Imodel;

fit.paramhat_MC = paramhat_MC;

end
