function [muhat,sigmahat,alphahat,betahat,I0hat,Ibhat,ss] = analyze_convolutedgamma(b,I,nConv,baseline,nFits)

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
lb_m = zeros(1,nConv);
ub_m = inf(1,nConv);
lb_s = zeros(1,nConv);
ub_s = inf(1,nConv);

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
    mean0 = randmeanD(b,I);
    std0 = randmeanD(b,I);
    
    mu0 = rand(1,nConv);
    mu0 = mu0 / sum(mu0) * mean0;
    
    var0 = rand(1,nConv);
    var0 = var0 / sum(var0) * std0^2;
    sigma0 = sqrt(var0);
    
    I00 = 1;
    if baseline
        Ib0 = abs(I(end))*rand();
    else
        Ib0 = 0;
    end

    param0 = [mu0 sigma0 I00 Ib0];
    
    paramhat_ = fmincon(@(param)sumofsquares_convolutedgamma(b,I,param(1:nConv),param(nConv+1:2*nConv),param(2*nConv+1),param(2*nConv+2)),param0,[],[],[],[],lb,ub,[],options);
    
    Imodel_ = signal_convolutedgamma(b,paramhat_(1:nConv),paramhat_(nConv+1:2*nConv),paramhat_(2*nConv+1),paramhat_(2*nConv+2));
    ss_ = sum( (I-Imodel_).^2 );
    
    if ss_ < ss
        paramhat = paramhat_;
        ss = ss_;
    end
end

muhat = paramhat(1:nConv) / bmax;
sigmahat = paramhat(nConv+1:2*nConv) / bmax;
I0hat = paramhat(2*nConv+1);
Ibhat = paramhat(2*nConv+2);

alphahat = muhat.^2./sigmahat.^2;
betahat = muhat./sigmahat.^2;

end
