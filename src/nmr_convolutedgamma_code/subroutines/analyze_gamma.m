function [muhat,sigmahat,alphahat,betahat,I0hat,Ibhat,ss] = analyze_gamma(b,I,baseline,nFits)

% Rescaling for numerical reasons.
bmax = max(b);
b = b / bmax;

% Optimization algorithm settings.
options = optimset( 'Display','off',...
                    'Algorithm','sqp',...
                    'GradObj','on',...
                    'DerivativeCheck','off',...
                    'MaxFunEvals',100000,...
                    'MaxIter',1000,...
                    'TolFun',1e-8,...
                    'TolX',1e-8);

% Lower and upper bounds for parameters.
lb_m = 0;
ub_m = inf;
lb_s = 0;
ub_s = inf;

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

for i = 1:nFits
    mu0 = randmeanD(b,I);
    sigma0 = randmeanD(b,I);
        
    I00 = 1;
    if baseline
        Ib0 = abs(I(end))*rand();
    else
        Ib0 = 0;
    end

    param0 = [mu0 sigma0 I00 Ib0];
    
    [paramhat_,ss_] = fmincon(@(param)sumofsquares_gamma(b,I,param(1),param(2),param(3),param(4)),param0,[],[],[],[],lb,ub,[],options);
    
    if ss_ < ss
        paramhat = paramhat_;
        ss = ss_;
    end
end

muhat = paramhat(1) / bmax;
sigmahat = paramhat(2) / bmax;
I0hat = paramhat(3);
Ibhat = paramhat(4);

alphahat = muhat.^2./sigmahat.^2;
betahat = muhat./sigmahat.^2;

end
