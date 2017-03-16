function [muhat, sigmahat, alphahat, betahat, I0hat, Ibhat, ss] = analyze(b, I, number_of_convolutions, baseline, number_of_fits)

% Rescaling for numerical reasons.
bmax = max(b);
b = b / bmax;

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

% Lower and upper bounds for parameters.
lb_m = zeros(1, number_of_convolutions);
ub_m = inf(1, number_of_convolutions);
lb_s = zeros(1, number_of_convolutions);
ub_s = inf(1, number_of_convolutions);

lb_I0 = 0;
ub_I0 = inf;

if baseline
    lb_Ib = 0;
    ub_Ib = inf;
else
    lb_Ib = 0;
    ub_Ib = 0;
end

lb = [lb_m lb_s lb_I0 lb_Ib];
ub = [ub_m ub_s ub_I0 ub_Ib];
            
% Fit model.
ss = inf;

for current_fit = 1:number_of_fits
    disp(current_fit)
    mean0 = rand_mean_D(b);
    std0 = rand_mean_D(b);
    
    mu0 = rand(1, number_of_convolutions);
    mu0 = mu0 / sum(mu0) * mean0;
    
    var0 = rand(1, number_of_convolutions);
    var0 = var0 / sum(var0) * std0^2;
    sigma0 = sqrt(var0);
    
    I00 = 1;
    if baseline
        Ib0 = abs(I(end)) * rand();
    else
        Ib0 = 0;
    end

    param_guess = [mu0 sigma0 I00 Ib0];
    
    paramhat_ = fmincon(@(param)sumofsquares(b, I, param(1:number_of_convolutions), param(number_of_convolutions+1:2*number_of_convolutions), param(2*number_of_convolutions+1), param(2*number_of_convolutions+2)), param_guess, [], [], [], [], lb, ub, [], options);
    
    Imodel_ = signal(b, paramhat_(1:number_of_convolutions), paramhat_(number_of_convolutions+1:2*number_of_convolutions), paramhat_(2*number_of_convolutions+1), paramhat_(2*number_of_convolutions+2));
    ss_ = sum( (I-Imodel_).^2 );
    
    if ss_ < ss
        paramhat = paramhat_;
        ss = ss_;
    end
end

muhat = paramhat(1:number_of_convolutions) / bmax;
sigmahat = paramhat(number_of_convolutions+1:2*number_of_convolutions) / bmax;
I0hat = paramhat(2*number_of_convolutions+1);
Ibhat = paramhat(2*number_of_convolutions+2);

alphahat = muhat.^2 ./ sigmahat.^2;
betahat = muhat ./ sigmahat.^2;

end
