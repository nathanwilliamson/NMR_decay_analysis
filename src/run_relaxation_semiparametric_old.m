%parametric model fitting routine for relaxation data, v 1.0
%written by Magnus Roding
% clear
clc
close all hidden

addpath('subroutines');

%% Initiate global random number stream.
seed    = round(sum(1e6*clock()));
s       = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(s);

%% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO BE EDITED BY USER.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1st option for data preparation
% filename_dx                = 't1t2.dx';
% 
% experiment                = readT1data(filename_dx);
%  
% t=experiment.t{1}
% I=experiment.I{1}
 
%% 2nd option for data preparation

I=load('I.txt');
t=load('t1.txt'); 

I=I/I(1)

% I = - I;
% I = I + max(abs(I));

%%
% Number of fittings (with random initializations).
nFits                       = 10;

% Number of Monte Carlo repetitions (0 = no error analysis).
nMC                         = 10;

% Type of model (combine exponential, stretched exponential, lognormal, and
% gamma freely). Specify more than one model for joint analysis and 
% comparison.
%type = {{'inversegamma'}};
type   = {{'inversegamma'},{'lognormal'},{'exponential','exponential'}};
%type   = {{'lognormal','lognormal'},{'inversegamma'},{'inversegamma','inversegamma'},{'lognormal'}};
%type = {{'exponential','exponential'}};
%type = {{'exponential'}};
% type = {{'stretchedexponential'}};

% Baseline or no baseline included in model?
baseline                    = false;

% Plot histograms of all the Monte Carlo values of parameters (to manually
% check for Gaussianity).
plotMChist                  = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF PARAMETERS TO BE EDITED BY USER.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Estimate.
nModels                     = numel(type);

fit                        = cell(nModels,1);
for currentModel = 1:nModels
    fit_                        = analyze(t,I,type{currentModel},baseline,nFits,nMC);
    fit{currentModel}           = fit_;
end

%% Plot fits and residuals.
close all

COLORS                      = get(groot,'DefaultAxesColorOrder');

fig                         = figure();
fig.Units                   = 'centimeters';
fig.PaperUnits              = 'centimeters';
fig.Position                = [0 0 12 9];
fig.PaperPosition           = [0 0 12 9];

ax1                         = subplot('Position',[0.15 0.45 0.8 0.50]);
ax1.FontSize                = 12;
ax1.Box                     = 'on';
ax1.YLabel.String           = 'Intensity';
% ax1.YScale                  = 'log';
ax1.XTickLabel              = {};
ax1.YLabel.Units            = 'centimeters';
ax1.YLim                    = [0 1.25*max(I)];
ax1.XLim                    = [0 max(t)];

ax1.YScale = 'log';

hl                          = line(t,I);
hl.Marker                   = 'o';
hl.Color                    = [0 0 0];
hl.LineStyle                = 'none';

legendStr                   = cell(nModels,1);
legendEntries               = [];
for currentModel = 1:nModels
    hl                          = line(t,fit{currentModel}.model);
    hl.Color                    = COLORS(currentModel,:);
    
    if numel(type{currentModel}) == 1
        legendStr{currentModel}     = type{currentModel}{1};
    else
        legendStr{currentModel}     = [];
        for i = 1:numel(type{currentModel})
            legendStr{currentModel} = [legendStr{currentModel} type{currentModel}{i} '+'];
        end
        legendStr{currentModel} = legendStr{currentModel}(1:end-1);
    end
    legendEntries(currentModel) = hl;
end
legend(legendEntries,legendStr)

ax2                         = subplot('Position',[0.15 0.20 0.8 0.15]);
ax2.FontSize                = 12;
ax2.Box                     = 'on';
ax2.XLabel.String           = 't (s)';
ax2.YLabel.String           = 'Residual';
ax2.YLabel.Units            = 'centimeters';
ax2.YLabel.Position(1)      = ax1.YLabel.Position(1);

hl                          = line(t,zeros(size(t)));
hl.Color                    = [0 0 0];

for currentModel = 1:nModels
    hl                          = line(t,fit{currentModel}.residuals);
    hl.Color                    = COLORS(currentModel,:);
    hl.Marker                   = 'o';
    hl.LineStyle                = 'none';
end

%% Plot distributions of relaxation times.
COLORS                      = get(groot,'DefaultAxesColorOrder');

fig                         = figure();
fig.Units                   = 'centimeters';
fig.PaperUnits              = 'centimeters';
fig.Position                = [0 0 12 9];
fig.PaperPosition           = [0 0 12 9];

ax1                         = subplot('Position',[0.15 0.20 0.8 0.75]);
ax1.FontSize                = 12;
ax1.Box                     = 'on';
ax1.XLabel.String           = 'T_1/T_2 (s)';
ax1.YLabel.String           = 'Probability';
% ax1.XScale                  = 'log';
ax1.XLim                    = [0 max(t)];

ax1.XScale = 'log';
tautau                      = linspace(0,max(t),10000);
% y                           = zeros(size(tautau));

legendStr                   = cell(nModels,1);
legendEntries               = [];
for currentModel = 1:nModels
    y = zeros(size(tautau));
    for currentComponent = 1:numel(type{currentModel})
        switch type{currentModel}{currentComponent}
            case 'exponential'
                hl = line(ones(1,2)*fit{currentModel}.components{currentComponent}.tau,[0 fit{currentModel}.components{currentComponent}.weight]);
                hl.Color                    = COLORS(currentModel,:);
            case 'stretchedexponential'
                hl = line(ones(1,2)*fit{currentModel}.components{currentComponent}.tau,[0 fit{currentModel}.components{currentComponent}.weight]);
                hl.Color                    = COLORS(currentModel,:);
            case 'lognormal'
                y = y + fit{currentModel}.components{currentComponent}.weight * lognpdf(tautau,fit{currentModel}.components{currentComponent}.mu,fit{currentModel}.components{currentComponent}.sigma);
            case 'inversegamma'
                p = fit{currentModel}.components{currentComponent}.beta^fit{currentModel}.components{currentComponent}.alpha./gamma(fit{currentModel}.components{currentComponent}.alpha).*tautau.^(-fit{currentModel}.components{currentComponent}.alpha-1).*exp(-fit{currentModel}.components{currentComponent}.beta./tautau);
                p(isnan(p)) = 0;
                y = y + fit{currentModel}.components{currentComponent}.weight * p;
        end
    end
    hl                          = line(tautau,y/max(y));
    hl.Color                    = COLORS(currentModel,:);
    
    if numel(type{currentModel}) == 1
        legendStr{currentModel}     = type{currentModel}{1};
    else
        legendStr{currentModel}     = [];
        for i = 1:numel(type{currentModel})
            legendStr{currentModel} = [legendStr{currentModel} type{currentModel}{i} '+'];
        end
        legendStr{currentModel} = legendStr{currentModel}(1:end-1);
    end
    legendEntries(currentModel) = hl;
end
legend(legendEntries,legendStr)

%% Plot histograms of Monte Carlo values for all parameters to check manually for Gaussianity.

if plotMChist
    for currentModel = 1:nModels
        for i = 1:size(fit{currentModel}.paramhat_MC,2)
            figure, hist(fit{currentModel}.paramhat_MC(:,i))
        end
    end
end

%% Present report of result in text form.

for currentModel = 1:nModels
    disp('==============================================================');
    if numel(type{currentModel}) == 1
        str     = type{currentModel}{1};
    else
        str     = [];
        for i = 1:numel(type{currentModel})
            str = [str type{currentModel}{i} '+'];
        end
        str = str(1:end-1);
    end
    disp(['MODEL ' num2str(currentModel) ': ' str])
    disp('==============================================================');
    disp('')
    disp(['Fit (residual sum of squares):   ' num2str(fit{currentModel}.ss)])
    disp('')
    disp(['I0:                              ' num2str(fit{currentModel}.I0)])
    disp('')
    if baseline 
        disp(['baseline:                        ' num2str(fit{currentModel}.baseline)])
        disp('')
    end
    for currentComponent = 1:numel(fit{currentModel}.components)
        disp('--------------------------------------------------------------');
        disp(['Component ' num2str(currentComponent) ': ' fit{currentModel}.components{currentComponent}.type])
        disp('')

        switch fit{currentModel}.components{currentComponent}.type
            case 'exponential'
                disp('Decay model:          I(t) = exp(-t/tau)')
                disp('')
                disp('Distribution model:   f(tau) = delta(tau-tau0)')
                disp('--------------------------------------------------------------')
                
                value               = fit{currentModel}.components{currentComponent}.tau;
                stddev              = fit{currentModel}.components{currentComponent}.std_tau;
                ci                  = [ fit{currentModel}.components{currentComponent}.tau - 1.960*fit{currentModel}.components{currentComponent}.std_tau , ...
                                        fit{currentModel}.components{currentComponent}.tau + 1.960*fit{currentModel}.components{currentComponent}.std_tau];
                disp('tau')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fit{currentModel}.components{currentComponent}.weight;
                stddev              = fit{currentModel}.components{currentComponent}.std_weight;
                ci                  = [ fit{currentModel}.components{currentComponent}.weight - 1.960*fit{currentModel}.components{currentComponent}.std_weight , ...
                                        fit{currentModel}.components{currentComponent}.weight + 1.960*fit{currentModel}.components{currentComponent}.std_weight];
                disp('weight')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            case 'stretchedexponential'
                disp('Decay model:          I(t) = exp(-(t/tau)^beta)')
                disp('')
                disp('Distribution model:   f(tau) unknown')
                disp('--------------------------------------------------------------')
                
                value               = fit{currentModel}.components{currentComponent}.tau;
                stddev              = fit{currentModel}.components{currentComponent}.std_tau;
                ci                  = [ fit{currentModel}.components{currentComponent}.tau - 1.960*fit{currentModel}.components{currentComponent}.std_tau , ...
                                        fit{currentModel}.components{currentComponent}.tau + 1.960*fit{currentModel}.components{currentComponent}.std_tau];
                disp('tau')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fit{currentModel}.components{currentComponent}.beta;
                stddev              = fit{currentModel}.components{currentComponent}.std_beta;
                ci                  = [ fit{currentModel}.components{currentComponent}.beta - 1.960*fit{currentModel}.components{currentComponent}.std_beta , ...
                                        fit{currentModel}.components{currentComponent}.beta + 1.960*fit{currentModel}.components{currentComponent}.std_beta];
                disp('beta')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fit{currentModel}.components{currentComponent}.weight;
                stddev              = fit{currentModel}.components{currentComponent}.std_weight;
                ci                  = [ fit{currentModel}.components{currentComponent}.weight - 1.960*fit{currentModel}.components{currentComponent}.std_weight , ...
                                        fit{currentModel}.components{currentComponent}.weight + 1.960*fit{currentModel}.components{currentComponent}.std_weight];
                disp('weight')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            case 'lognormal'
                disp('Decay model:          I(t) = integral(f(tau),0,inf) (numerical approx.)')
                disp('')
                disp('Distribution model:   f(tau) = 1/(tau*sigma*sqrt(2*pi))*exp(-(log(tau)-mu)^2/(2*sigma^2))')
                disp('--------------------------------------------------------------')
                
                value               = fit{currentModel}.components{currentComponent}.mu;
                stddev              = fit{currentModel}.components{currentComponent}.std_mu;
                ci                  = [ fit{currentModel}.components{currentComponent}.mu - 1.960*fit{currentModel}.components{currentComponent}.std_mu , ...
                                        fit{currentModel}.components{currentComponent}.mu + 1.960*fit{currentModel}.components{currentComponent}.std_mu];
                disp('mu')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fit{currentModel}.components{currentComponent}.sigma;
                stddev              = fit{currentModel}.components{currentComponent}.std_sigma;
                ci                  = [ fit{currentModel}.components{currentComponent}.sigma - 1.960*fit{currentModel}.components{currentComponent}.std_sigma , ...
                                        fit{currentModel}.components{currentComponent}.sigma + 1.960*fit{currentModel}.components{currentComponent}.std_sigma];
                disp('sigma')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fit{currentModel}.components{currentComponent}.meantau;
                stddev              = fit{currentModel}.components{currentComponent}.std_meantau;
                ci                  = [ fit{currentModel}.components{currentComponent}.meantau - 1.960*fit{currentModel}.components{currentComponent}.std_meantau , ...
                                        fit{currentModel}.components{currentComponent}.meantau + 1.960*fit{currentModel}.components{currentComponent}.std_meantau];
                disp('meantau')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fit{currentModel}.components{currentComponent}.stdtau;
                stddev              = fit{currentModel}.components{currentComponent}.std_stdtau;
                ci                  = [ fit{currentModel}.components{currentComponent}.stdtau - 1.960*fit{currentModel}.components{currentComponent}.std_stdtau , ...
                                        fit{currentModel}.components{currentComponent}.stdtau + 1.960*fit{currentModel}.components{currentComponent}.std_stdtau];
                disp('stdtau')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fit{currentModel}.components{currentComponent}.spreadtau;
                stddev              = fit{currentModel}.components{currentComponent}.std_spreadtau;
                ci                  = [ fit{currentModel}.components{currentComponent}.spreadtau - 1.960*fit{currentModel}.components{currentComponent}.std_spreadtau , ...
                                        fit{currentModel}.components{currentComponent}.spreadtau + 1.960*fit{currentModel}.components{currentComponent}.std_spreadtau];
                disp('spreadtau (meantau/stdtau, coefficient of variation)')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fit{currentModel}.components{currentComponent}.modetau;
                stddev              = fit{currentModel}.components{currentComponent}.std_modetau;
                ci                  = [ fit{currentModel}.components{currentComponent}.modetau - 1.960*fit{currentModel}.components{currentComponent}.std_modetau , ...
                                        fit{currentModel}.components{currentComponent}.modetau + 1.960*fit{currentModel}.components{currentComponent}.std_modetau];
                disp('modetau (peak of distribution)')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fit{currentModel}.components{currentComponent}.weight;
                stddev              = fit{currentModel}.components{currentComponent}.std_weight;
                ci                  = [ fit{currentModel}.components{currentComponent}.weight - 1.960*fit{currentModel}.components{currentComponent}.std_weight , ...
                                        fit{currentModel}.components{currentComponent}.weight + 1.960*fit{currentModel}.components{currentComponent}.std_weight];
                disp('weight')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            case 'inversegamma'
                disp('Decay model:          I(t) = (beta/(beta+t))^alpha')
                disp('')
                disp('Distribution model:   f(tau) = beta^alpha/gamma(alpha)*tau^(-alpha-1)*exp(-beta/tau)')
                disp('--------------------------------------------------------------')
                
                value               = fit{currentModel}.components{currentComponent}.alpha;
                stddev              = fit{currentModel}.components{currentComponent}.std_alpha;
                ci                  = [ fit{currentModel}.components{currentComponent}.alpha - 1.960*fit{currentModel}.components{currentComponent}.std_alpha , ...
                                        fit{currentModel}.components{currentComponent}.alpha + 1.960*fit{currentModel}.components{currentComponent}.std_alpha];
                disp('alpha')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fit{currentModel}.components{currentComponent}.beta;
                stddev              = fit{currentModel}.components{currentComponent}.std_beta;
                ci                  = [ fit{currentModel}.components{currentComponent}.beta - 1.960*fit{currentModel}.components{currentComponent}.std_beta , ...
                                        fit{currentModel}.components{currentComponent}.beta + 1.960*fit{currentModel}.components{currentComponent}.std_beta];
                disp('beta')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                if fit{currentModel}.components{currentComponent}.alpha > 1
                    value               = fit{currentModel}.components{currentComponent}.meantau;
                    stddev              = fit{currentModel}.components{currentComponent}.std_meantau;
                    ci                  = [ fit{currentModel}.components{currentComponent}.meantau - 1.960*fit{currentModel}.components{currentComponent}.std_meantau , ...
                                            fit{currentModel}.components{currentComponent}.meantau + 1.960*fit{currentModel}.components{currentComponent}.std_meantau];
                    disp('meantau')
                    disp(['      - Value:                   ' num2str(value)])
                    disp(['      - Std:                     ' num2str(stddev)])
                    disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                else
                    disp('meantau')
                    disp(['      - Value:                   ' '(alpha <= 1 -> Undefined)'])
                    disp(['      - Std:                     ' '(alpha <= 1 -> Undefined)'])
                    disp(['      - 95 % CI:                 ' '(alpha <= 1 -> Undefined)'])
                end
                
                if fit{currentModel}.components{currentComponent}.alpha > 2
                    value               = fit{currentModel}.components{currentComponent}.stdtau;
                    stddev              = fit{currentModel}.components{currentComponent}.std_stdtau;
                    ci                  = [ fit{currentModel}.components{currentComponent}.stdtau - 1.960*fit{currentModel}.components{currentComponent}.std_stdtau , ...
                        fit{currentModel}.components{currentComponent}.stdtau + 1.960*fit{currentModel}.components{currentComponent}.std_stdtau];
                    disp('stdtau')
                    disp(['      - Value:                   ' num2str(value)])
                    disp(['      - Std:                     ' num2str(stddev)])
                    disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                else
                    disp('stdtau')
                    disp(['      - Value:                   ' '(alpha <= 2 -> Undefined)'])
                    disp(['      - Std:                     ' '(alpha <= 2 -> Undefined)'])
                    disp(['      - 95 % CI:                 ' '(alpha <= 2 -> Undefined)'])
                end
                
                if fit{currentModel}.components{currentComponent}.alpha > 2
                    value               = fit{currentModel}.components{currentComponent}.spreadtau;
                    stddev              = fit{currentModel}.components{currentComponent}.std_spreadtau;
                    ci                  = [ fit{currentModel}.components{currentComponent}.spreadtau - 1.960*fit{currentModel}.components{currentComponent}.std_spreadtau , ...
                        fit{currentModel}.components{currentComponent}.spreadtau + 1.960*fit{currentModel}.components{currentComponent}.std_spreadtau];
                    disp('spreadtau (meantau/stdtau, coefficient of variation)')
                    disp(['      - Value:                   ' num2str(value)])
                    disp(['      - Std:                     ' num2str(stddev)])
                    disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                else
                    disp('stdtau')
                    disp(['      - Value:                   ' '(alpha <= 2 -> Undefined)'])
                    disp(['      - Std:                     ' '(alpha <= 2 -> Undefined)'])
                    disp(['      - 95 % CI:                 ' '(alpha <= 2 -> Undefined)'])
                end
                
                value               = fit{currentModel}.components{currentComponent}.modetau;
                stddev              = fit{currentModel}.components{currentComponent}.std_modetau;
                ci                  = [ fit{currentModel}.components{currentComponent}.modetau - 1.960*fit{currentModel}.components{currentComponent}.std_modetau , ...
                                        fit{currentModel}.components{currentComponent}.modetau + 1.960*fit{currentModel}.components{currentComponent}.std_modetau];
                disp('modetau (peak of distribution)')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fit{currentModel}.components{currentComponent}.weight;
                stddev              = fit{currentModel}.components{currentComponent}.std_weight;
                ci                  = [ fit{currentModel}.components{currentComponent}.weight - 1.960*fit{currentModel}.components{currentComponent}.std_weight , ...
                                        fit{currentModel}.components{currentComponent}.weight + 1.960*fit{currentModel}.components{currentComponent}.std_weight];
                disp('weight')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
        end
    end
end

