function [] = print_results(fit)

disp('==============================================================================================');
disp('Results')
disp('==============================================================================================');
disp('')
disp(['Fit (residual sum of squares):   ' num2str(fit.ss)])
disp('')
disp(['I0:                              ' num2str(fit.I0)])
disp('')
disp(['baseline:                        ' num2str(fit.baseline)])
disp('')
for currentComponent = 1:numel(fit.components)
    disp('----------------------------------------------------------------------------------------------');
    disp(['Component ' num2str(currentComponent) ': ' fit.components{currentComponent}.model])
    disp('')
    
    switch fit.components{currentComponent}.model
        case 'exponential'
            disp('Decay model:          I(b) = exp(-b*D)')
            disp('')
            disp('Distribution model:   f(D) = delta(D-D0)')
            disp('----------------------------------------------------------------------------------------------');
            
            value = fit.components{currentComponent}.D;
            stddev = fit.components{currentComponent}.std_D;
            ci = [ fit.components{currentComponent}.D - 1.960*fit.components{currentComponent}.std_D , ...
                fit.components{currentComponent}.D + 1.960*fit.components{currentComponent}.std_D];
            disp('D')
            disp(['      - Value:                   ' num2str(value)])
            disp(['      - Std:                     ' num2str(stddev)])
            disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            
            value = fit.components{currentComponent}.theta;
            stddev = fit.components{currentComponent}.std_theta;
            ci = [ fit.components{currentComponent}.theta - 1.960*fit.components{currentComponent}.std_theta , ...
                fit.components{currentComponent}.theta + 1.960*fit.components{currentComponent}.std_theta];
            disp('theta')
            disp(['      - Value:                   ' num2str(value)])
            disp(['      - Std:                     ' num2str(stddev)])
            disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
        case 'stretchedexponential'
            disp('Decay model:          I(b) = exp(-(b*D)^beta)')
            disp('')
            disp('Distribution model:   f(D) unknown')
            disp('----------------------------------------------------------------------------------------------');
            
            value = fit.components{currentComponent}.D;
            stddev = fit.components{currentComponent}.std_D;
            ci = [ fit.components{currentComponent}.D - 1.960*fit.components{currentComponent}.std_D , ...
                fit.components{currentComponent}.D + 1.960*fit.components{currentComponent}.std_D];
            disp('D')
            disp(['      - Value:                   ' num2str(value)])
            disp(['      - Std:                     ' num2str(stddev)])
            disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            
            value = fit.components{currentComponent}.beta;
            stddev = fit.components{currentComponent}.std_beta;
            ci = [ fit.components{currentComponent}.beta - 1.960*fit.components{currentComponent}.std_beta , ...
                fit.components{currentComponent}.beta + 1.960*fit.components{currentComponent}.std_beta];
            disp('beta')
            disp(['      - Value:                   ' num2str(value)])
            disp(['      - Std:                     ' num2str(stddev)])
            disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            
            value = fit.components{currentComponent}.theta;
            stddev = fit.components{currentComponent}.std_theta;
            ci = [ fit.components{currentComponent}.theta - 1.960*fit.components{currentComponent}.std_theta , ...
                fit.components{currentComponent}.theta + 1.960*fit.components{currentComponent}.std_theta];
            disp('theta')
            disp(['      - Value:                   ' num2str(value)])
            disp(['      - Std:                     ' num2str(stddev)])
            disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
        case 'lognormal'
            disp('Decay model:          I(b) = integral(f(D),0,inf) (numerical approx.)')
            disp('')
            disp('Distribution model:   f(D) = 1/(D*sigma*sqrt(2*pi))*exp(-(log(D)-mu)^2/(2*sigma^2))')
            disp('----------------------------------------------------------------------------------------------');
            
            value = fit.components{currentComponent}.mu;
            stddev = fit.components{currentComponent}.std_mu;
            ci = [ fit.components{currentComponent}.mu - 1.960*fit.components{currentComponent}.std_mu , ...
                fit.components{currentComponent}.mu + 1.960*fit.components{currentComponent}.std_mu];
            disp('mu')
            disp(['      - Value:                   ' num2str(value)])
            disp(['      - Std:                     ' num2str(stddev)])
            disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            
            value = fit.components{currentComponent}.sigma;
            stddev = fit.components{currentComponent}.std_sigma;
            ci = [ fit.components{currentComponent}.sigma - 1.960*fit.components{currentComponent}.std_sigma , ...
                fit.components{currentComponent}.sigma + 1.960*fit.components{currentComponent}.std_sigma];
            disp('sigma')
            disp(['      - Value:                   ' num2str(value)])
            disp(['      - Std:                     ' num2str(stddev)])
            disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            
            value = fit.components{currentComponent}.meanD;
            stddev = fit.components{currentComponent}.std_meanD;
            ci = [ fit.components{currentComponent}.meanD - 1.960*fit.components{currentComponent}.std_meanD , ...
                fit.components{currentComponent}.meanD + 1.960*fit.components{currentComponent}.std_meanD];
            disp('meanD')
            disp(['      - Value:                   ' num2str(value)])
            disp(['      - Std:                     ' num2str(stddev)])
            disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            
            value = fit.components{currentComponent}.stdD;
            stddev = fit.components{currentComponent}.std_stdD;
            ci = [ fit.components{currentComponent}.stdD - 1.960*fit.components{currentComponent}.std_stdD , ...
                fit.components{currentComponent}.stdD + 1.960*fit.components{currentComponent}.std_stdD];
            disp('stdD')
            disp(['      - Value:                   ' num2str(value)])
            disp(['      - Std:                     ' num2str(stddev)])
            disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            
            value = fit.components{currentComponent}.spreadD;
            stddev = fit.components{currentComponent}.std_spreadD;
            ci = [ fit.components{currentComponent}.spreadD - 1.960*fit.components{currentComponent}.std_spreadD , ...
                fit.components{currentComponent}.spreadD + 1.960*fit.components{currentComponent}.std_spreadD];
            disp('spreadD (stdD/meanD, coefficient of variation)')
            disp(['      - Value:                   ' num2str(value)])
            disp(['      - Std:                     ' num2str(stddev)])
            disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            
            value = fit.components{currentComponent}.modeD;
            stddev = fit.components{currentComponent}.std_modeD;
            ci = [ fit.components{currentComponent}.modeD - 1.960*fit.components{currentComponent}.std_modeD , ...
                fit.components{currentComponent}.modeD + 1.960*fit.components{currentComponent}.std_modeD];
            disp('modeD (peak of distribution)')
            disp(['      - Value:                   ' num2str(value)])
            disp(['      - Std:                     ' num2str(stddev)])
            disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            
            value = fit.components{currentComponent}.theta;
            stddev = fit.components{currentComponent}.std_theta;
            ci = [ fit.components{currentComponent}.theta - 1.960*fit.components{currentComponent}.std_theta , ...
                fit.components{currentComponent}.theta + 1.960*fit.components{currentComponent}.std_theta];
            disp('theta')
            disp(['      - Value:                   ' num2str(value)])
            disp(['      - Std:                     ' num2str(stddev)])
            disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
        case 'gamma'
            disp('Decay model:          I(b) = (beta/(beta+b))^alpha')
            disp('')
            disp('Distribution model:   f(D) = beta^alpha/gamma(alpha)*D^(alpha-1)*exp(-beta*D)')
            disp('----------------------------------------------------------------------------------------------');
            
            value = fit.components{currentComponent}.alpha;
            stddev = fit.components{currentComponent}.std_alpha;
            ci = [ fit.components{currentComponent}.alpha - 1.960*fit.components{currentComponent}.std_alpha , ...
                fit.components{currentComponent}.alpha + 1.960*fit.components{currentComponent}.std_alpha];
            disp('alpha')
            disp(['      - Value:                   ' num2str(value)])
            disp(['      - Std:                     ' num2str(stddev)])
            disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            
            value = fit.components{currentComponent}.beta;
            stddev = fit.components{currentComponent}.std_beta;
            ci = [ fit.components{currentComponent}.beta - 1.960*fit.components{currentComponent}.std_beta , ...
                fit.components{currentComponent}.beta + 1.960*fit.components{currentComponent}.std_beta];
            disp('beta')
            disp(['      - Value:                   ' num2str(value)])
            disp(['      - Std:                     ' num2str(stddev)])
            disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            
            value = fit.components{currentComponent}.meanD;
            stddev = fit.components{currentComponent}.std_meanD;
            ci = [ fit.components{currentComponent}.meanD - 1.960*fit.components{currentComponent}.std_meanD , ...
                fit.components{currentComponent}.meanD + 1.960*fit.components{currentComponent}.std_meanD];
            disp('meanD')
            disp(['      - Value:                   ' num2str(value)])
            disp(['      - Std:                     ' num2str(stddev)])
            disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            
            value = fit.components{currentComponent}.stdD;
            stddev = fit.components{currentComponent}.std_stdD;
            ci = [ fit.components{currentComponent}.stdD - 1.960*fit.components{currentComponent}.std_stdD , ...
                fit.components{currentComponent}.stdD + 1.960*fit.components{currentComponent}.std_stdD];
            disp('stdD')
            disp(['      - Value:                   ' num2str(value)])
            disp(['      - Std:                     ' num2str(stddev)])
            disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            
            value = fit.components{currentComponent}.spreadD;
            stddev = fit.components{currentComponent}.std_spreadD;
            ci = [ fit.components{currentComponent}.spreadD - 1.960*fit.components{currentComponent}.std_spreadD , ...
                fit.components{currentComponent}.spreadD + 1.960*fit.components{currentComponent}.std_spreadD];
            disp('spreadD (stdD/meanD, coefficient of variation)')
            disp(['      - Value:                   ' num2str(value)])
            disp(['      - Std:                     ' num2str(stddev)])
            disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            
            value = fit.components{currentComponent}.modeD;
            stddev = fit.components{currentComponent}.std_modeD;
            ci = [ fit.components{currentComponent}.modeD - 1.960*fit.components{currentComponent}.std_modeD , ...
                fit.components{currentComponent}.modeD + 1.960*fit.components{currentComponent}.std_modeD];
            disp('modeD (peak of distribution)')
            disp(['      - Value:                   ' num2str(value)])
            disp(['      - Std:                     ' num2str(stddev)])
            disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            
            value = fit.components{currentComponent}.theta;
            stddev = fit.components{currentComponent}.std_theta;
            ci = [ fit.components{currentComponent}.theta - 1.960*fit.components{currentComponent}.std_theta , ...
                fit.components{currentComponent}.theta + 1.960*fit.components{currentComponent}.std_theta];
            disp('theta')
            disp(['      - Value:                   ' num2str(value)])
            disp(['      - Std:                     ' num2str(stddev)])
            disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
    end
end
end
