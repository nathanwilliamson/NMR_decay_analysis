function fit = analyze(t,I,type,baseline,nFits,nMC)

TolFun                          = 1e-8;
TolX                            = 1e-8;

% Fit model.
ss                              = inf;

nComponents                     = numel(type);

for i = 1:nFits
    param0                          = [];
    % Generate initial values of parameters.
    for currentComponent = 1:nComponents
        switch type{currentComponent}
            case 'exponential'
                taurnd                  = randmeantau(t);
                param0                  = [param0 taurnd];
            case 'stretchedexponential'
                taurnd                  = randmeantau(t);
                betarnd                 = rand();
                param0                  = [param0 taurnd betarnd];
            case 'lognormal'
                Mrnd                    = randmeantau(t);
                CVrnd                   = 0.05 + 1 * rand();
                murnd                   = log(Mrnd) - 1/2*log(1+CVrnd^2);
                sigmarnd                = sqrt(log(1+CVrnd^2));
                param0                  = [param0 murnd sigmarnd];    
            case 'inversegamma'
                Mrnd                    = randmeantau(t);
                CVrnd                   = 0.05 + 0.75 * rand();
                alpharnd                = 2 + 1/CVrnd^2;
                betarnd                 = Mrnd*(alpharnd-1);
                param0                  = [param0 alpharnd betarnd];
        end
    end
    if baseline
        param0                      = [param0 abs(I(end))*randn()];
    end
    theta0                          = rand(1,nComponents);
    theta0                          = theta0/sum(theta0);
    param0                          = [param0 theta0];
    param0                          = [param0 max(I)*rand()];

    % Fit model.
    paramhat_                       = estimate(t,I,type,baseline,param0,TolFun,TolX);
    
    % Compute residual sum of squares.
    Imodel_                         = signal(t,type,baseline,paramhat_);
    ss_                             = sum( (I-Imodel_).^2 );
    
    if ss_ < ss
        paramhat                    = paramhat_;
        ss                          = ss_;
        Imodel                      = Imodel_;
    end
end

% Error analysis.
if nMC < 2
    paramhat_MC = nan(size(paramhat));
else
    paramhat_MC = zeros(nMC,numel(paramhat));
    
    for i = 1:nMC
        disp(i)
        sigma_residual              = sqrt(ss/numel(t));
        I_MC                        = I + sigma_residual*randn(size(I));
        paramhat_MC(i,:)            = estimate(t,I_MC,type,baseline,paramhat,TolFun,TolX);
    end
end

% Structure output.
fit                         = struct();
fit.ss                      = ss;

fit.I0                      = paramhat(end);
fit.std_I0                  = std(paramhat_MC(:,end),[],1);

weights                     = paramhat(end-nComponents:end-1);
weights_MC                  = paramhat_MC(:,end-nComponents:end-1);
std_weights                 = std(weights_MC,[],1);

fit.components              = cell(nComponents,1);
for currentComponent = 1:nComponents
    fit.components{currentComponent} = struct();
end

ind                         = 1;
for currentComponent = 1:nComponents
    
    fit.components{currentComponent}.type       = type{currentComponent};
    
    switch type{currentComponent}
        case 'exponential'
            tau                                             = paramhat(ind);
            tau_MC                                          = paramhat_MC(:,ind);
            
            fit.components{currentComponent}.tau            = tau;
            fit.components{currentComponent}.std_tau        = std(tau_MC,[],1);
            fit.components{currentComponent}.tau_MC         = tau_MC;
            
            fit.components{currentComponent}.meantau        = tau;
            fit.components{currentComponent}.std_meantau    = std(tau_MC,[],1);
            fit.components{currentComponent}.meantau_MC     = tau_MC;
            
            fit.components{currentComponent}.stdtau         = 0;
            fit.components{currentComponent}.std_stdtau     = 0;
            fit.components{currentComponent}.stdtau_MC      = zeros(nMC,1);
            
            fit.components{currentComponent}.spreadtau      = 0;
            fit.components{currentComponent}.std_spreadtau  = 0;
            fit.components{currentComponent}.spreadtau_MC   = zeros(nMC,1);
            
            fit.components{currentComponent}.modetau        = tau;
            fit.components{currentComponent}.std_modetau    = std(tau_MC,[],1);
            fit.components{currentComponent}.modetau_MC     = tau_MC;
            
            fit.components{currentComponent}.weight         = weights(currentComponent);
            fit.components{currentComponent}.std_weight     = std_weights(currentComponent);
            fit.components{currentComponent}.weight_MC      = weights_MC(:,currentComponent);  
            
            ind                                             = ind + 1;
        case 'stretchedexponential'
            tau                                             = paramhat(ind);
            tau_MC                                          = paramhat_MC(:,ind);
            
            beta                                            = paramhat(ind+1);
            beta_MC                                         = paramhat_MC(:,ind+1);
                        
            fit.components{currentComponent}.tau            = tau;
            fit.components{currentComponent}.std_tau        = std(tau_MC,[],1);
            fit.components{currentComponent}.tau_MC         = tau_MC;
            
            fit.components{currentComponent}.beta           = beta;
            fit.components{currentComponent}.std_beta       = std(beta_MC,[],1);
            fit.components{currentComponent}.beta_MC        = beta_MC;
            
            fit.components{currentComponent}.meantau        = nan;
            fit.components{currentComponent}.std_meantau    = nan;
            fit.components{currentComponent}.meantau_MC     = nan(nMC,1);
            
            fit.components{currentComponent}.stdtau         = nan;
            fit.components{currentComponent}.std_stdtau     = nan;
            fit.components{currentComponent}.stdtau_MC      = nan(nMC,1);
            
            fit.components{currentComponent}.spreadtau      = nan;
            fit.components{currentComponent}.std_spreadtau  = nan;
            fit.components{currentComponent}.spreadtau_MC   = nan(nMC,1);
            
            fit.components{currentComponent}.modetau        = nan;
            fit.components{currentComponent}.std_modetau    = nan;
            fit.components{currentComponent}.modetau_MC     = nan(nMC,1);
            
            fit.components{currentComponent}.weight         = weights(currentComponent);
            fit.components{currentComponent}.std_weight     = std_weights(currentComponent);
            fit.components{currentComponent}.weight_MC      = weights_MC(:,currentComponent); 
            
            ind                                             = ind + 2;
        case 'lognormal'
            mu                                              = paramhat(ind);
            mu_MC                                           = paramhat_MC(:,ind);
            
            sigma                                           = paramhat(ind+1);
            sigma_MC                                        = paramhat_MC(:,ind+1);

            fit.components{currentComponent}.mu             = mu;
            fit.components{currentComponent}.std_mu         = std(mu_MC,[],1);
            fit.components{currentComponent}.mu_MC          = mu_MC;
            
            fit.components{currentComponent}.sigma          = sigma;
            fit.components{currentComponent}.std_sigma      = std(sigma_MC,[],1);
            fit.components{currentComponent}.sigma_MC       = sigma_MC;
            
            meantau                                         = exp(mu+sigma^2/2);
            meantau_MC                                      = exp(mu_MC+sigma_MC.^2/2);
            
            stdtau                                          = sqrt(exp(sigma^2)-1)*exp(mu+sigma^2/2); 
            stdtau_MC                                       = sqrt(exp(sigma_MC.^2)-1).*exp(mu_MC+sigma_MC.^2/2);
            
            spreadtau                                       = sqrt(exp(sigma^2)-1);
            spreadtau_MC                                    = sqrt(exp(sigma_MC.^2)-1);
            
            modetau                                         = exp(mu-sigma^2);
            modetau_MC                                      = exp(mu_MC-sigma_MC.^2);
            
            fit.components{currentComponent}.meantau        = meantau; 
            fit.components{currentComponent}.std_meantau    = std(meantau_MC,[],1);
            fit.components{currentComponent}.meantau_MC     = meantau_MC;
            
            fit.components{currentComponent}.stdtau         = stdtau;
            fit.components{currentComponent}.std_stdtau     = std(stdtau_MC,[],1);
            fit.components{currentComponent}.stdtau_MC      = stdtau_MC;
            
            fit.components{currentComponent}.spreadtau      = spreadtau;
            fit.components{currentComponent}.std_spreadtau  = std(spreadtau_MC,[],1);
            fit.components{currentComponent}.spreadtau_MC   = spreadtau_MC;
            
            fit.components{currentComponent}.modetau        = modetau;
            fit.components{currentComponent}.std_modetau    = std(modetau_MC,[],1);
            fit.components{currentComponent}.modetau_MC     = modetau_MC;
            
            fit.components{currentComponent}.weight         = weights(currentComponent);
            fit.components{currentComponent}.std_weight     = std_weights(currentComponent);
            fit.components{currentComponent}.weight_MC      = weights_MC(:,currentComponent);
            
            ind                                             = ind + 2;
        case 'inversegamma'
            alpha                                           = paramhat(ind);
            alpha_MC                                        = paramhat_MC(:,ind);
            
            beta                                            = paramhat(ind+1);
            beta_MC                                         = paramhat_MC(:,ind+1);
            
            fit.components{currentComponent}.alpha          = alpha;
            fit.components{currentComponent}.std_alpha      = std(alpha_MC,[],1);
            fit.components{currentComponent}.alpha_MC       = alpha_MC;
            
            fit.components{currentComponent}.beta           = beta;
            fit.components{currentComponent}.std_beta       = std(beta_MC,[],1);
            fit.components{currentComponent}.beta_MC        = beta_MC;
            
            meantau                                         = beta./(alpha-1);
            meantau_MC                                      = beta_MC./(alpha_MC-1);
            
            stdtau                                          = beta./(alpha-1)./sqrt(alpha-2);
            stdtau_MC                                       = beta_MC./(alpha_MC-1)./sqrt(alpha_MC-2);
            
            spreadtau                                       = 1./sqrt(alpha-2);
            spreadtau_MC                                    = 1./sqrt(alpha_MC-2);

            modetau                                         = beta./(alpha+1);
            modetau_MC                                      = beta_MC./(alpha_MC+1);
            
            fit.components{currentComponent}.meantau        = meantau; 
            fit.components{currentComponent}.std_meantau    = std(meantau_MC,[],1);
            fit.components{currentComponent}.meantau_MC     = meantau_MC;
            
            fit.components{currentComponent}.stdtau         = stdtau;
            fit.components{currentComponent}.std_stdtau     = std(stdtau_MC,[],1);
            fit.components{currentComponent}.stdtau_MC      = stdtau_MC;
            
            fit.components{currentComponent}.spreadtau      = spreadtau;
            fit.components{currentComponent}.std_spreadtau  = std(spreadtau_MC,[],1);
            fit.components{currentComponent}.spreadtau_MC   = spreadtau_MC;
            
            fit.components{currentComponent}.modetau        = modetau;
            fit.components{currentComponent}.std_modetau    = std(modetau_MC,[],1);
            fit.components{currentComponent}.modetau_MC     = modetau_MC;
            
            fit.components{currentComponent}.weight         = weights(currentComponent);
            fit.components{currentComponent}.std_weight     = std_weights(currentComponent);
            fit.components{currentComponent}.weight_MC      = weights_MC(:,currentComponent);
            
            ind                                             = ind + 2;
    end
end
if baseline
    fit.baseline                                            = paramhat(ind);
    fit.std_baseline                                        = std(paramhat_MC(:,ind),[],1);
    fit.baseline_MC                                         = paramhat_MC(:,ind);
end 

fit.model                                                   = Imodel;
fit.residuals                                               = I - Imodel;

fit.paramhat_MC                                             = paramhat_MC;

end
