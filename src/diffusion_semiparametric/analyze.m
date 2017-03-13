function fit = analyze(b,I,model,baseline,nFits,nMC)

nComponents                     = numel(model);

% Rescaling for numerical reasons.
kmax                            = max(b);
b                               = b / kmax;

% Optimization algorithm settings.
options                         = optimoptions('fmincon');
options.Algorithm               = 'sqp';
options.Display                 = 'off';
options.MaxFunEvals             = 10000;
options.MaxIter                 = 1000;
options.TolFun                  = 1e-8;
options.TolX                    = 1e-8;

% Lower and upper bounds for parameters.
lb                              = [];
ub                              = [];

lb_theta                        = zeros(1,nComponents);
ub_theta                        = ones(1,nComponents);

for currentComponent = 1:nComponents
    switch model{currentComponent}{1}
        case 'exponential'
            lb_D = 0;
            ub_D = inf;
            
            if numel(model{currentComponent}) > 1
                for i = 2:2:numel(model{currentComponent})-1
                    switch model{currentComponent}{i}
                        case 'D'
                            switch numel(model{currentComponent}{i+1})
                                case 1
                                    lb_D = model{currentComponent}{i+1}(1);
                                    ub_D = model{currentComponent}{i+1}(1);
                                case 2
                                    lb_D = model{currentComponent}{i+1}(1);
                                    ub_D = model{currentComponent}{i+1}(2);
                            end
                        case 'theta'
                            switch numel(model{currentComponent}{i+1})
                                case 1
                                    lb_theta(currentComponent) = model{currentComponent}{i+1}(1);
                                    ub_theta(currentComponent) = model{currentComponent}{i+1}(1);
                                case 2
                                    lb_theta(currentComponent) = model{currentComponent}{i+1}(1);
                                    ub_theta(currentComponent) = model{currentComponent}{i+1}(2);
                            end
                    end
                end
            end
            lb                      = [lb lb_D*kmax];
            ub                      = [ub ub_D*kmax];
        case 'stretchedexponential'
            lb_D = 0;
            ub_D = inf;
            lb_beta = 0;
            ub_beta = 1;
            
            if numel(model{currentComponent}) > 1
                for i = 2:2:numel(model{currentComponent})-1
                    switch model{currentComponent}{i}
                        case 'D'
                            switch numel(model{currentComponent}{i+1})
                                case 1
                                    lb_D = model{currentComponent}{i+1}(1);
                                    ub_D = model{currentComponent}{i+1}(1);
                                case 2
                                    lb_D = model{currentComponent}{i+1}(1);
                                    ub_D = model{currentComponent}{i+1}(2);
                            end
                        case 'beta'
                            switch numel(model{currentComponent}{i+1})
                                case 1
                                    lb_beta = model{currentComponent}{i+1}(1);
                                    ub_beta = model{currentComponent}{i+1}(1);
                                case 2
                                    lb_beta = model{currentComponent}{i+1}(1);
                                    ub_beta = model{currentComponent}{i+1}(2);
                            end
                        case 'theta'
                            switch numel(model{currentComponent}{i+1})
                                case 1
                                    lb_theta(currentComponent) = model{currentComponent}{i+1}(1);
                                    ub_theta(currentComponent) = model{currentComponent}{i+1}(1);
                                case 2
                                    lb_theta(currentComponent) = model{currentComponent}{i+1}(1);
                                    ub_theta(currentComponent) = model{currentComponent}{i+1}(2);
                            end
                    end
                end
            end
            lb                      = [lb lb_D*kmax lb_beta];
            ub                      = [ub ub_D*kmax ub_beta];
        case 'lognormal'
            cv_min                  = 0.01;
            lb_mu = -inf;
            ub_mu = inf;
            lb_sigma = sqrt(log(1+cv_min^2));
            ub_sigma = inf;
            
            if numel(model{currentComponent}) > 1
                for i = 2:2:numel(model{currentComponent})-1
                    switch model{currentComponent}{i}
                        case 'mu'
                            switch numel(model{currentComponent}{i+1})
                                case 1
                                    lb_mu = model{currentComponent}{i+1}(1);
                                    ub_mu = model{currentComponent}{i+1}(1);
                                case 2
                                    lb_mu = model{currentComponent}{i+1}(1);
                                    ub_mu = model{currentComponent}{i+1}(2);
                            end
                        case 'sigma'
                            switch numel(model{currentComponent}{i+1})
                                case 1
                                    lb_sigma = model{currentComponent}{i+1}(1);
                                    ub_sigma = model{currentComponent}{i+1}(1);
                                case 2
                                    lb_sigma = model{currentComponent}{i+1}(1);
                                    ub_sigma = model{currentComponent}{i+1}(2);
                            end
                        case 'theta'
                            switch numel(model{currentComponent}{i+1})
                                case 1
                                    lb_theta(currentComponent) = model{currentComponent}{i+1}(1);
                                    ub_theta(currentComponent) = model{currentComponent}{i+1}(1);
                                case 2
                                    lb_theta(currentComponent) = model{currentComponent}{i+1}(1);
                                    ub_theta(currentComponent) = model{currentComponent}{i+1}(2);
                            end
                    end
                end
            end
            lb                      = [lb lb_mu+log(kmax) lb_sigma];
            ub                      = [ub ub_mu+log(kmax) ub_sigma];
        case 'gamma'
            lb_alpha = 0;
            ub_alpha = inf;
            lb_beta = 0;
            ub_beta = inf;
            
            if numel(model{currentComponent}) > 1
                for i = 2:2:numel(model{currentComponent})-1
                    switch model{currentComponent}{i}
                        case 'alpha'
                            switch numel(model{currentComponent}{i+1})
                                case 1
                                    lb_alpha = model{currentComponent}{i+1}(1);
                                    ub_alpha = model{currentComponent}{i+1}(1);
                                case 2
                                    lb_alpha = model{currentComponent}{i+1}(1);
                                    ub_alpha = model{currentComponent}{i+1}(2);
                            end
                        case 'beta'
                            switch numel(model{currentComponent}{i+1})
                                case 1
                                    lb_beta = model{currentComponent}{i+1}(1);
                                    ub_beta = model{currentComponent}{i+1}(1);
                                case 2
                                    lb_beta = model{currentComponent}{i+1}(1);
                                    ub_beta = model{currentComponent}{i+1}(2);
                            end
                        case 'theta'
                            switch numel(model{currentComponent}{i+1})
                                case 1
                                    lb_theta(currentComponent) = model{currentComponent}{i+1}(1);
                                    ub_theta(currentComponent) = model{currentComponent}{i+1}(1);
                                case 2
                                    lb_theta(currentComponent) = model{currentComponent}{i+1}(1);
                                    ub_theta(currentComponent) = model{currentComponent}{i+1}(2);
                            end
                    end
                end
            end
            lb                      = [lb lb_alpha lb_beta/kmax];
            ub                      = [ub ub_alpha ub_beta/kmax];
    end
end

if baseline
    lb_b = -inf;
    ub_b = inf;
else
    lb_b = 0;
    ub_b = 0;
end
lb                              = [lb lb_b];
ub                              = [ub ub_b];
            
lb                              = [lb lb_theta]; 
ub                              = [ub ub_theta]; 

lb_I0                           = 0;
ub_I0                           = inf;
lb                              = [lb lb_I0];
ub                              = [ub ub_I0];

Aeq                             = zeros(size(lb));
Aeq(end-nComponents:end-1)      = 1;
beq                             = 1;

% Fit model.
ss                              = inf;

nComponents                     = numel(model);

for i = 1:nFits
    param0                          = [];
    
    % Generate initial values of parameters.
    ind = 1;
    for currentComponent = 1:nComponents
        switch model{currentComponent}{1}
            case 'exponential'
                Drnd                    = randmeanD(b,I);
                Drnd                    = max(lb(ind),Drnd);
                Drnd                    = min(ub(ind),Drnd);
                param0(ind)             = Drnd;
                ind                     = ind + 1;
            case 'stretchedexponential'
                Drnd                    = randmeanD(b,I);
                Drnd                    = max(lb(ind),Drnd);
                Drnd                    = min(ub(ind),Drnd);
                param0(ind)             = Drnd;
                ind                     = ind + 1;
                
                betarnd                 = rand();
                betarnd                 = max(lb(ind),betarnd);
                betarnd                 = min(ub(ind),betarnd);
                param0(ind)             = betarnd;
                ind                     = ind + 1;
            case 'lognormal'
                Mrnd                    = randmeanD(b,I);
                CVrnd                   = 0.05 + 1 * rand();
                
                murnd                   = log(Mrnd) - 1/2*log(1+CVrnd^2);
                murnd                   = max(lb(ind),murnd);
                murnd                   = min(ub(ind),murnd);
                param0(ind)             = murnd;
                ind                     = ind + 1;                
                
                sigmarnd                = sqrt(log(1+CVrnd^2));
                sigmarnd                = max(lb(ind),sigmarnd);
                sigmarnd                = min(ub(ind),sigmarnd);
                param0(ind)             = sigmarnd;
                ind                     = ind + 1;   
            case 'gamma'
                Mrnd                    = randmeanD(b,I);
                CVrnd                   = 0.05 + 0.75 * rand();
                
                alpharnd                = 1/CVrnd^2;
                alpharnd                = max(lb(ind),alpharnd);
                alpharnd                = min(ub(ind),alpharnd);
                param0(ind)             = alpharnd;
                ind                     = ind + 1;                
                
                betarnd                 = 1/(Mrnd*CVrnd^2);
                betarnd                 = max(lb(ind),betarnd);
                betarnd                 = min(ub(ind),betarnd);
                param0(ind)             = betarnd;
                ind                     = ind + 1;  
        end
    end
    brnd                            = abs(I(end))*randn();
    brnd                            = max(lb(ind),brnd);
    brnd                            = min(ub(ind),brnd);
    param0(ind)                     = brnd;
    
    theta0                          = rand(1,nComponents);
    theta0                          = theta0/sum(theta0); % This can be out of bounds, but will be corrected by the algorithm.
    param0                          = [param0 theta0];
    param0                          = [param0 max(I)*rand()];
%     return

    % Fit model.
    paramhat_                       = fmincon(@(param)sumofsquares(b,I,model,param),param0,[],[],Aeq,beq,lb,ub,[],options);
    
    % Compute residual sum of squares.
    Imodel_                         = signal(b,model,paramhat_);
    ss_                             = sum( (I-Imodel_).^2 );
    
    if ss_ < ss
        paramhat                    = paramhat_;
        ss                          = ss_;
        Imodel                      = Imodel_;
    end
end

% Error analysis.
sigma_residual = sqrt(ss/numel(b));
if nMC < 2
    paramhat_MC = nan(size(paramhat));
else
    paramhat_MC = zeros(nMC,numel(paramhat));
    
    for i = 1:nMC
        disp(i)
        I_MC                        = I + sigma_residual*randn(size(I));
        paramhat_MC(i,:)            = fmincon(@(param)sumofsquares(b,I_MC,model,param),paramhat,[],[],Aeq,beq,lb,ub,[],options);
    end
end

% Structure output.
fit                         = struct();
fit.ss                      = ss;

fit.I0                      = paramhat(end);
fit.std_I0                  = std(paramhat_MC(:,end),[],1);

theta                       = paramhat(end-nComponents:end-1);
theta_MC                    = paramhat_MC(:,end-nComponents:end-1);
std_theta                   = std(theta_MC,[],1);

fit.components              = cell(nComponents,1);
for currentComponent = 1:nComponents
    fit.components{currentComponent} = struct();
end

ind                         = 1;
for currentComponent = 1:nComponents
    
    fit.components{currentComponent}.model       = model{currentComponent}{1};
    
    switch model{currentComponent}{1}
        case 'exponential'
            D                                               = paramhat(ind) / kmax;
            D_MC                                            = paramhat_MC(:,ind) / kmax;
            
            fit.components{currentComponent}.D              = D;
            fit.components{currentComponent}.std_D          = std(D_MC,[],1);
            fit.components{currentComponent}.D_MC           = D_MC;
            
            fit.components{currentComponent}.meanD          = D;
            fit.components{currentComponent}.std_meanD      = std(D_MC,[],1);
            fit.components{currentComponent}.meanD_MC       = D_MC;
            
            fit.components{currentComponent}.stdD           = 0;
            fit.components{currentComponent}.std_stdD       = 0;
            fit.components{currentComponent}.stdD_MC        = zeros(nMC,1);
            
            fit.components{currentComponent}.spreadD        = 0;
            fit.components{currentComponent}.std_spreadD    = 0;
            fit.components{currentComponent}.spreadD_MC     = zeros(nMC,1);
            
            fit.components{currentComponent}.modeD          = D;
            fit.components{currentComponent}.std_modeD      = std(D_MC,[],1);
            fit.components{currentComponent}.modeD_MC       = D_MC;
            
            fit.components{currentComponent}.theta          = theta(currentComponent);
            fit.components{currentComponent}.std_theta      = std_theta(currentComponent);
            fit.components{currentComponent}.theta_MC       = theta_MC(:,currentComponent);  
            
            ind                                             = ind + 1;
        case 'stretchedexponential'
            D                                               = paramhat(ind) / kmax;
            D_MC                                            = paramhat_MC(:,ind) / kmax;
            
            beta                                            = paramhat(ind+1);
            beta_MC                                         = paramhat_MC(:,ind+1);
                        
            fit.components{currentComponent}.D              = D;
            fit.components{currentComponent}.std_D          = std(D_MC,[],1);
            fit.components{currentComponent}.D_MC           = D_MC;
            
            fit.components{currentComponent}.beta           = beta;
            fit.components{currentComponent}.std_beta       = std(beta_MC,[],1);
            fit.components{currentComponent}.beta_MC        = beta_MC;
            
            fit.components{currentComponent}.meanD          = nan;
            fit.components{currentComponent}.std_meanD      = nan;
            fit.components{currentComponent}.meanD_MC       = nan(nMC,1);
            
            fit.components{currentComponent}.stdD           = nan;
            fit.components{currentComponent}.std_stdD       = nan;
            fit.components{currentComponent}.stdD_MC        = nan(nMC,1);
            
            fit.components{currentComponent}.spreadD        = nan;
            fit.components{currentComponent}.std_spreadD    = nan;
            fit.components{currentComponent}.spreadD_MC     = nan(nMC,1);
            
            fit.components{currentComponent}.modeD          = nan;
            fit.components{currentComponent}.std_modeD      = nan;
            fit.components{currentComponent}.modeD_MC       = nan(nMC,1);
            
            fit.components{currentComponent}.theta          = theta(currentComponent);
            fit.components{currentComponent}.std_theta      = std_theta(currentComponent);
            fit.components{currentComponent}.theta_MC       = theta_MC(:,currentComponent); 
            
            ind                                             = ind + 2;
        case 'lognormal'
            mu                                              = paramhat(ind) - log(kmax);
            mu_MC                                           = paramhat_MC(:,ind)- log(kmax);
            
            sigma                                           = paramhat(ind+1);
            sigma_MC                                        = paramhat_MC(:,ind+1);

            fit.components{currentComponent}.mu             = mu;
            fit.components{currentComponent}.std_mu         = std(mu_MC,[],1);
            fit.components{currentComponent}.mu_MC          = mu_MC;
            
            fit.components{currentComponent}.sigma          = sigma;
            fit.components{currentComponent}.std_sigma      = std(sigma_MC,[],1);
            fit.components{currentComponent}.sigma_MC       = sigma_MC;
            
            meanD                                           = exp(mu+sigma^2/2);
            meanD_MC                                        = exp(mu_MC+sigma_MC.^2/2);
            
            stdD                                            = sqrt(exp(sigma^2)-1)*exp(mu+sigma^2/2); 
            stdD_MC                                         = sqrt(exp(sigma_MC.^2)-1).*exp(mu_MC+sigma_MC.^2/2);
            
            spreadD                                         = sqrt(exp(sigma^2)-1);
            spreadD_MC                                      = sqrt(exp(sigma_MC.^2)-1);
            
            modeD                                           = exp(mu-sigma^2);
            modeD_MC                                        = exp(mu_MC-sigma_MC.^2);
            
            fit.components{currentComponent}.meanD          = meanD; 
            fit.components{currentComponent}.std_meanD      = std(meanD_MC,[],1);
            fit.components{currentComponent}.meanD_MC       = meanD_MC;
            
            fit.components{currentComponent}.stdD           = stdD;
            fit.components{currentComponent}.std_stdD       = std(stdD_MC,[],1);
            fit.components{currentComponent}.stdD_MC        = stdD_MC;
            
            fit.components{currentComponent}.spreadD        = spreadD;
            fit.components{currentComponent}.std_spreadD    = std(spreadD_MC,[],1);
            fit.components{currentComponent}.spreadD_MC     = spreadD_MC;
            
            fit.components{currentComponent}.modeD          = modeD;
            fit.components{currentComponent}.std_modeD      = std(modeD_MC,[],1);
            fit.components{currentComponent}.modeD_MC       = modeD_MC;
            
            fit.components{currentComponent}.theta          = theta(currentComponent);
            fit.components{currentComponent}.std_theta      = std_theta(currentComponent);
            fit.components{currentComponent}.theta_MC       = theta_MC(:,currentComponent);
            
            ind                                             = ind + 2;
        case 'gamma'
            alpha                                           = paramhat(ind);
            alpha_MC                                        = paramhat_MC(:,ind);
            
            beta                                            = paramhat(ind+1) * kmax;
            beta_MC                                         = paramhat_MC(:,ind+1) * kmax;
                        
            fit.components{currentComponent}.alpha          = alpha;
            fit.components{currentComponent}.std_alpha      = std(alpha_MC,[],1);
            fit.components{currentComponent}.alpha_MC       = alpha_MC;
            
            fit.components{currentComponent}.beta           = beta;
            fit.components{currentComponent}.std_beta       = std(beta_MC,[],1);
            fit.components{currentComponent}.beta_MC        = beta_MC;
            
            meanD                                           = alpha/beta;
            meanD_MC                                        = alpha_MC./beta_MC;
            
            stdD                                            = sqrt(alpha)/beta;
            stdD_MC                                         = sqrt(alpha_MC)./beta_MC;
            
            spreadD                                         = 1/sqrt(alpha);
            spreadD_MC                                      = 1./sqrt(alpha_MC);

            modeD                                           = (alpha-1)/beta;
            modeD_MC                                        = (alpha_MC-1)./beta_MC;
            
            fit.components{currentComponent}.meanD          = meanD; 
            fit.components{currentComponent}.std_meanD      = std(meanD_MC,[],1);
            fit.components{currentComponent}.meanD_MC       = meanD_MC;
            
            fit.components{currentComponent}.stdD           = stdD;
            fit.components{currentComponent}.std_stdD       = std(stdD_MC,[],1);
            fit.components{currentComponent}.stdD_MC        = stdD_MC;
            
            fit.components{currentComponent}.spreadD        = spreadD;
            fit.components{currentComponent}.std_spreadD    = std(spreadD_MC,[],1);
            fit.components{currentComponent}.spreadD_MC     = spreadD_MC;
            
            fit.components{currentComponent}.modeD          = modeD;
            fit.components{currentComponent}.std_modeD      = std(modeD_MC,[],1);
            fit.components{currentComponent}.modeD_MC       = modeD_MC;
            
            fit.components{currentComponent}.theta          = theta(currentComponent);
            fit.components{currentComponent}.std_theta      = std_theta(currentComponent);
            fit.components{currentComponent}.theta_MC       = theta_MC(:,currentComponent);
            
            ind                                             = ind + 2;
    end
end
fit.baseline                                                = paramhat(ind);
fit.std_baseline                                            = std(paramhat_MC(:,ind),[],1);
fit.baseline_MC                                             = paramhat_MC(:,ind);

fit.Imodel                                                  = Imodel;
fit.residuals                                               = I - Imodel;

fit.paramhat_MC                                             = paramhat_MC;

end
