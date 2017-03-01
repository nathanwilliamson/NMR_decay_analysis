function [muhat,sigmahat,I0hat,Ibhat,ss] = analyze_lognormal(b,I,baseline,nFits)

% Rescaling for numerical reasons.
bmax = max(b);
b = b / bmax;

% Optimization algorithm settings.
options = optimset( 'Display','off',...
                    'Algorithm','sqp',...
                    'GradObj','off',...
                    'DerivativeCheck','off',...
                    'MaxFunEvals',100000,...
                    'MaxIter',1000,...
                    'TolFun',1e-8,...
                    'TolX',1e-8);

% Lower and upper bounds for parameters.
lb_mu = -inf;
ub_mu = inf;
lb_sigma = zeros;
ub_sigma = inf;

lb_I0 = 0;
ub_I0 = inf;

if baseline
    lb_Ib = 0;
    ub_Ib = inf;
else
    lb_Ib = 0;
    ub_Ib = 0;
end

lb = [lb_mu lb_sigma lb_I0 lb_Ib];
ub = [ub_mu ub_sigma ub_I0 ub_Ib];
            
% Fit model.
ss = inf;

for i = 1:nFits
    mean0 = randmeanD(b,I);
    spread0 = 0.05 + rand();
    
    mu0 = log(mean0) - 1/2*log(1+spread0^2);
    sigma0 = sqrt(log(1+spread0^2));
    
    I00 = 1;
    if baseline
        Ib0 = abs(I(end))*rand();
    else
        Ib0 = 0;
    end

    param0 = [mu0 sigma0 I00 Ib0];
    
    paramhat_ = fmincon(@(param)sumofsquares_lognormal(b,I,param(1),param(2),param(3),param(4)),param0,[],[],[],[],lb,ub,[],options);
    
    Imodel_ = signal_lognormal(b,paramhat_(1),paramhat_(2),paramhat_(3),paramhat_(4));
    ss_  = sum( (I-Imodel_).^2 );
    
    if ss_ < ss
        paramhat = paramhat_;
        ss = ss_;
    end
end

muhat = paramhat(1) - log(bmax);
sigmahat = paramhat(2);
I0hat = paramhat(3);
Ibhat = paramhat(4);

end
