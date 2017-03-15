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
            disp('Decay model:          I(t) = exp(-t/T)')
            disp('')
            disp('Distribution model:   f(T) = delta(T-T0)')
            disp('----------------------------------------------------------------------------------------------');
            
            value = fit.components{currentComponent}.T;
            stddev = fit.components{currentComponent}.std_T;
            ci = [ fit.components{currentComponent}.T - 1.960*fit.components{currentComponent}.std_T , ...
                fit.components{currentComponent}.T + 1.960*fit.components{currentComponent}.std_T];
            disp('T')
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
            disp('Decay model:          I(t) = exp(-(t/T)^beta)')
            disp('')
            disp('Distribution model:   f(T) unknown')
            disp('----------------------------------------------------------------------------------------------');
            
            value = fit.components{currentComponent}.T;
            stddev = fit.components{currentComponent}.std_T;
            ci = [ fit.components{currentComponent}.T - 1.960*fit.components{currentComponent}.std_T , ...
                fit.components{currentComponent}.T + 1.960*fit.components{currentComponent}.std_T];
            disp('T')
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
            disp('Decay model:          I(t) = integral(f(T),0,inf) (numerical approx.)')
            disp('')
            disp('Distribution model:   f(T) = 1/(T*sigma*sqrt(2*pi))*exp(-(log(T)-mu)^2/(2*sigma^2))')
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
            
            value = fit.components{currentComponent}.meanT;
            stddev = fit.components{currentComponent}.std_meanT;
            ci = [ fit.components{currentComponent}.meanT - 1.960*fit.components{currentComponent}.std_meanT , ...
                fit.components{currentComponent}.meanT + 1.960*fit.components{currentComponent}.std_meanT];
            disp('meanT')
            disp(['      - Value:                   ' num2str(value)])
            disp(['      - Std:                     ' num2str(stddev)])
            disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            
            value = fit.components{currentComponent}.stdT;
            stddev = fit.components{currentComponent}.std_stdT;
            ci = [ fit.components{currentComponent}.stdT - 1.960*fit.components{currentComponent}.std_stdT , ...
                fit.components{currentComponent}.stdT + 1.960*fit.components{currentComponent}.std_stdT];
            disp('stdT')
            disp(['      - Value:                   ' num2str(value)])
            disp(['      - Std:                     ' num2str(stddev)])
            disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            
            value = fit.components{currentComponent}.spreadT;
            stddev = fit.components{currentComponent}.std_spreadT;
            ci = [ fit.components{currentComponent}.spreadT - 1.960*fit.components{currentComponent}.std_spreadT , ...
                fit.components{currentComponent}.spreadT + 1.960*fit.components{currentComponent}.std_spreadT];
            disp('spreadT (stdT/meanT, coefficient of variation)')
            disp(['      - Value:                   ' num2str(value)])
            disp(['      - Std:                     ' num2str(stddev)])
            disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            
            value = fit.components{currentComponent}.modeT;
            stddev = fit.components{currentComponent}.std_modeT;
            ci = [ fit.components{currentComponent}.modeT - 1.960*fit.components{currentComponent}.std_modeT , ...
                fit.components{currentComponent}.modeT + 1.960*fit.components{currentComponent}.std_modeT];
            disp('modeT (peak of distribution)')
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
        case 'inversegamma'
            disp('Decay model:          I(t) = (beta/(beta+t))^alpha')
            disp('')
            disp('Distribution model:   f(T) = beta^alpha/gamma(alpha)*T^(-alpha-1)*exp(-beta/T)')
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
            
            value = fit.components{currentComponent}.meanT;
            stddev = fit.components{currentComponent}.std_meanT;
            ci = [ fit.components{currentComponent}.meanT - 1.960*fit.components{currentComponent}.std_meanT , ...
                fit.components{currentComponent}.meanT + 1.960*fit.components{currentComponent}.std_meanT];
            disp('meanT')
            disp(['      - Value:                   ' num2str(value)])
            disp(['      - Std:                     ' num2str(stddev)])
            disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            
            value = fit.components{currentComponent}.stdT;
            stddev = fit.components{currentComponent}.std_stdT;
            ci = [ fit.components{currentComponent}.stdT - 1.960*fit.components{currentComponent}.std_stdT , ...
                fit.components{currentComponent}.stdT + 1.960*fit.components{currentComponent}.std_stdT];
            disp('stdT')
            disp(['      - Value:                   ' num2str(value)])
            disp(['      - Std:                     ' num2str(stddev)])
            disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            
            value = fit.components{currentComponent}.spreadT;
            stddev = fit.components{currentComponent}.std_spreadT;
            ci = [ fit.components{currentComponent}.spreadT - 1.960*fit.components{currentComponent}.std_spreadT , ...
                fit.components{currentComponent}.spreadT + 1.960*fit.components{currentComponent}.std_spreadT];
            disp('spreadT (stdT/meanT, coefficient of variation)')
            disp(['      - Value:                   ' num2str(value)])
            disp(['      - Std:                     ' num2str(stddev)])
            disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            
            value = fit.components{currentComponent}.modeT;
            stddev = fit.components{currentComponent}.std_modeT;
            ci = [ fit.components{currentComponent}.modeT - 1.960*fit.components{currentComponent}.std_modeT , ...
                fit.components{currentComponent}.modeT + 1.960*fit.components{currentComponent}.std_modeT];
            disp('modeT (peak of distribution)')
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
