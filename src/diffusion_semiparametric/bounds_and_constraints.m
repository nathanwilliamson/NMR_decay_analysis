function [lb, ub, Aeq, beq] = bounds_and_constraints(bmax, model, baseline)

lb = [];
ub = [];

lb_theta = zeros(1, number_of_components);
ub_theta = ones(1, number_of_components);

for currentComponent = 1:number_of_components
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
            lb                      = [lb lb_D*bmax];
            ub                      = [ub ub_D*bmax];
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
            lb                      = [lb lb_D*bmax lb_beta];
            ub                      = [ub ub_D*bmax ub_beta];
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
            lb                      = [lb lb_mu+log(bmax) lb_sigma];
            ub                      = [ub ub_mu+log(bmax) ub_sigma];
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
            lb                      = [lb lb_alpha lb_beta/bmax];
            ub                      = [ub ub_alpha ub_beta/bmax];
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
Aeq(end-number_of_components:end-1)      = 1;
beq                             = 1;