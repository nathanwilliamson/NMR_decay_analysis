function paramhat = estimate(t,I,type,baseline,param0,TolFun,TolX)

options = optimset( 'Display','off',...
                    'Algorithm','sqp',...
                    'GradObj','off',...
                    'DerivativeCheck','off',...
                    'MaxFunEvals',10000,...
                    'MaxIter',1000,...
                    'TolFun',TolFun,...
                    'TolX',TolX);

nComponents                     = numel(type);

lb                              = [];
ub                              = [];

for currentComponent = 1:nComponents
    switch type{currentComponent}
        case 'exponential'
            lb                      = [lb 0];
            %lb                      = [lb 1];
            ub                      = [ub inf];
        case 'stretchedexponential'
            lb                      = [lb 0 0];
            ub                      = [ub inf 1];
        case 'lognormal'
            cv_min                  = 0.01;
            lb                      = [lb -inf sqrt(log(1+cv_min^2))];
            ub                      = [ub inf inf];
        case 'inversegamma'
            lb                      = [lb 0 0];
            ub                      = [ub inf inf];
    end
end

if baseline
    lb                      = [lb -inf];
    ub                      = [ub inf];
end 
            
% theta.
lb                              = [lb zeros(1,nComponents)]; 
ub                              = [ub ones(1,nComponents)]; 
%lb                              = [lb 0.8176 0.1881];  %for Song T1, fixing the weights
%ub                              = [ub 0.8178 0.1883]; %for Song T1, fixing the weights
%lb                              = [lb 0.70 0.289];  %for Song T1, fixing the weights
%ub                              = [ub 0.71 0.291]; %for Song T1, fixing the weights
% I0.
lb                              = [lb 0];
ub                              = [ub inf];

% theta constraints.
Aeq                             = zeros(size(lb));
Aeq(end-nComponents:end-1)      = 1;
beq                             = 1;
    
% Estimate.
paramhat                        = fmincon(@(param)sumofsquares(t,I,type,baseline,param),param0,[],[],Aeq,beq,lb,ub,[],options);

end
