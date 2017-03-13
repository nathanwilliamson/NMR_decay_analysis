%% Initialization.
clear
clc
close all hidden

%% Path to subroutines.
addpath('diffusion_semiparametric');

%% Random stream.
random_seed = round( sum( 1e6 * clock() ) );
s = RandStream('mt19937ar','Seed',random_seed);
RandStream.setGlobalStream(s);

%% Simulate data set.



%% end 2nd option
I=I/I(1);
% Select which points in the data set to include in fitting.
pointsInd                   = 1:numel(b);%%[1:2,4:numel(b)];% 
b                           = b(pointsInd);
I                           = I(pointsInd);

% return

% Number of fittings (with random initializations).
nFits                       = 10;

% Number of Monte Carlo repetitions (0 = no error analysis).
nMC                         = 10;

% Type of model (combine exponential, stretched exponential, lognormal, and
% gamma freely). Specify more than one model for joint analysis and 
% comparison.
%models = {    {{'stretchedexponential'}} };
%models = {    {{'exponential'},{'exponential'}}     };
models = {    {{'gamma'}},{{'lognormal'}},{{'exponential'},{'exponential'}}     };
%models = {    {{'lognormal'},{'lognormal'}}     };
%models = {   {{'lognormal'}} };
%models = {   {{'lognormal'},{'exponential','D',[1.7e-9 6e-9]}} };
%models = {   {{'gamma'},{'exponential','D',[1.7e-9 4e-9]}} };
%models = {    {{'lognormal','mu',[-23.17 -23.14],'sigma',[0.31 0.32]},{'exponential'}}     };
%models = {    {{'gamma','alpha',9.217,'beta',2.802E10}}     }; 
% Baseline or no baseline included in model(s)?
baseline                    = false;

% Plot histograms of all the Monte Carlo values of parameters (to manually
% check for Gaussianity).
plotMChist                  = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF PARAMETERS TO BE EDITED BY USER.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Estimate.
nModels                     = numel(models);

fit                        = cell(nModels,1);
for currentModel = 1:nModels
    fit_                        = analyze(b,I,models{currentModel},baseline,nFits,nMC);
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
ax1.YScale                  = 'log';
ax1.XTickLabel              = {};
ax1.YLabel.Units            = 'centimeters';
ax1.YLim                    = [max(1e-4,0.5*min(I)) 1.25*max(I)];

hl                          = line(b,I);
hl.Marker                   = 'o';
hl.Color                    = [0 0 0];
hl.LineStyle                = 'none';

legendStr                   = cell(nModels,1);
legendEntries               = [];
for currentModel = 1:nModels
    hl                          = line(b,fit{currentModel}.Imodel);
    hl.Color                    = COLORS(currentModel,:);
    
    if numel(models{currentModel}) == 1
        legendStr{currentModel}     = models{currentModel}{1}{1};
    else
        legendStr{currentModel}     = [];
        for i = 1:numel(models{currentModel})
            legendStr{currentModel} = [legendStr{currentModel} models{currentModel}{i}{1} '+'];
        end
        legendStr{currentModel} = legendStr{currentModel}(1:end-1);
    end
    legendEntries(currentModel) = hl;
end
legend(legendEntries,legendStr)

ax2                         = subplot('Position',[0.15 0.20 0.8 0.15]);
ax2.FontSize                = 12;
ax2.Box                     = 'on';
ax2.XLabel.String           = 'b (sm^{-2})';
ax2.YLabel.String           = 'Residual';
ax2.YLabel.Units            = 'centimeters';
ax2.YLabel.Position(1)      = ax1.YLabel.Position(1);

hl                          = line(b,zeros(size(b)));
hl.Color                    = [0 0 0];

for currentModel = 1:nModels
    hl                          = line(b,fit{currentModel}.residuals);
    hl.Color                    = COLORS(currentModel,:);
    hl.Marker                   = 'o';
    hl.LineStyle                = 'none';
end

%% Plot distributions of diffusion coefficients.
COLORS                      = get(groot,'DefaultAxesColorOrder');

fig                         = figure();
fig.Units                   = 'centimeters';
fig.PaperUnits              = 'centimeters';
fig.Position                = [0 0 12 9];
fig.PaperPosition           = [0 0 12 9];

ax1                         = subplot('Position',[0.15 0.20 0.8 0.75]);
ax1.FontSize                = 12;
ax1.Box                     = 'on';
ax1.XLabel.String           = 'D (m^2s^{-1})';
ax1.YLabel.String           = 'Probability';
ax1.XScale                  = 'log';

minD                        = inf;
maxD                        = 0;

for currentModel = 1:nModels
    for currentComponent = 1:numel(models{currentModel})
        switch models{currentModel}{currentComponent}{1}
            case 'exponential'
                minD = min(minD,fit{currentModel}.components{currentComponent}.D);
                maxD = max(maxD,fit{currentModel}.components{currentComponent}.D);
            case 'stretchedexponential'
                minD = min(minD,fit{currentModel}.components{currentComponent}.D);
                maxD = max(maxD,fit{currentModel}.components{currentComponent}.D);
            case 'lognormal'
                minD = min(minD, logninv(1e-6,fit{currentModel}.components{currentComponent}.mu,fit{currentModel}.components{currentComponent}.sigma));
                maxD = max(maxD, logninv(1-1e-6,fit{currentModel}.components{currentComponent}.mu,fit{currentModel}.components{currentComponent}.sigma));
            case 'gamma'
                minD = min(minD, gaminv(1e-6,fit{currentModel}.components{currentComponent}.alpha,1/fit{currentModel}.components{currentComponent}.beta));
                maxD = max(maxD, gaminv(1-1e-6,fit{currentModel}.components{currentComponent}.alpha,1/fit{currentModel}.components{currentComponent}.beta));
        end
    end
end

DD                          = logspace(log10(minD),log10(maxD),10000);
y                           = zeros(size(DD));

legendStr                   = cell(nModels,1);
legendEntries               = [];
for currentModel = 1:nModels
    y = zeros(size(DD));
    for currentComponent = 1:numel(models{currentModel})
        switch models{currentModel}{currentComponent}{1}
            case 'exponential'
                hl = line(ones(1,2)*fit{currentModel}.components{currentComponent}.D,[0 fit{currentModel}.components{currentComponent}.theta]);
                hl.Color                    = COLORS(currentModel,:);
            case 'stretchedexponential'
                hl = line(ones(1,2)*fit{currentModel}.components{currentComponent}.D,[0 fit{currentModel}.components{currentComponent}.theta]);
                hl.Color                    = COLORS(currentModel,:);
            case 'lognormal'
                y = y + fit{currentModel}.components{currentComponent}.theta * lognpdf(DD,fit{currentModel}.components{currentComponent}.mu,fit{currentModel}.components{currentComponent}.sigma);
            case 'gamma'
                y = y + fit{currentModel}.components{currentComponent}.theta * gampdf(DD,fit{currentModel}.components{currentComponent}.alpha,1/fit{currentModel}.components{currentComponent}.beta);
        end
    end
    hl                          = line(DD,y/max(y));
    hl.Color                    = COLORS(currentModel,:);
    
    if numel(models{currentModel}) == 1
        legendStr{currentModel}     = models{currentModel}{1}{1};
    else
        legendStr{currentModel}     = [];
        for i = 1:numel(models{currentModel})
            legendStr{currentModel} = [legendStr{currentModel} models{currentModel}{i}{1} '+'];
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
    if numel(models{currentModel}) == 1
        str     = models{currentModel}{1}{1};
    else
        str     = [];
        for i = 1:numel(models{currentModel})
            str = [str models{currentModel}{i}{1} '+'];
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
    disp(['baseline:                        ' num2str(fit{currentModel}.baseline)])
    disp('')
    for currentComponent = 1:numel(fit{currentModel}.components)
        disp('--------------------------------------------------------------');
        disp(['Component ' num2str(currentComponent) ': ' fit{currentModel}.components{currentComponent}.model])
        disp('')

        switch fit{currentModel}.components{currentComponent}.model
            case 'exponential'
                disp('Decay model:          I(b) = exp(-b*D)')
                disp('')
                disp('Distribution model:   f(D) = delta(D-D0)')
                disp('--------------------------------------------------------------')
                
                value               = fit{currentModel}.components{currentComponent}.D;
                stddev              = fit{currentModel}.components{currentComponent}.std_D;
                ci                  = [ fit{currentModel}.components{currentComponent}.D - 1.960*fit{currentModel}.components{currentComponent}.std_D , ...
                                        fit{currentModel}.components{currentComponent}.D + 1.960*fit{currentModel}.components{currentComponent}.std_D];
                disp('D')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fit{currentModel}.components{currentComponent}.theta;
                stddev              = fit{currentModel}.components{currentComponent}.std_theta;
                ci                  = [ fit{currentModel}.components{currentComponent}.theta - 1.960*fit{currentModel}.components{currentComponent}.std_theta , ...
                                        fit{currentModel}.components{currentComponent}.theta + 1.960*fit{currentModel}.components{currentComponent}.std_theta];
                disp('theta')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            case 'stretchedexponential'
                disp('Decay model:          I(b) = exp(-(b*D)^beta)')
                disp('')
                disp('Distribution model:   f(D) unknown')
                disp('--------------------------------------------------------------')
                
                value               = fit{currentModel}.components{currentComponent}.D;
                stddev              = fit{currentModel}.components{currentComponent}.std_D;
                ci                  = [ fit{currentModel}.components{currentComponent}.D - 1.960*fit{currentModel}.components{currentComponent}.std_D , ...
                                        fit{currentModel}.components{currentComponent}.D + 1.960*fit{currentModel}.components{currentComponent}.std_D];
                disp('D')
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
                
                value               = fit{currentModel}.components{currentComponent}.theta;
                stddev              = fit{currentModel}.components{currentComponent}.std_theta;
                ci                  = [ fit{currentModel}.components{currentComponent}.theta - 1.960*fit{currentModel}.components{currentComponent}.std_theta , ...
                                        fit{currentModel}.components{currentComponent}.theta + 1.960*fit{currentModel}.components{currentComponent}.std_theta];
                disp('theta')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            case 'lognormal'
                disp('Decay model:          I(b) = integral(f(D),0,inf) (numerical approx.)')
                disp('')
                disp('Distribution model:   f(D) = 1/(D*sigma*sqrt(2*pi))*exp(-(log(D)-mu)^2/(2*sigma^2))')
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
                
                value               = fit{currentModel}.components{currentComponent}.meanD;
                stddev              = fit{currentModel}.components{currentComponent}.std_meanD;
                ci                  = [ fit{currentModel}.components{currentComponent}.meanD - 1.960*fit{currentModel}.components{currentComponent}.std_meanD , ...
                                        fit{currentModel}.components{currentComponent}.meanD + 1.960*fit{currentModel}.components{currentComponent}.std_meanD];
                disp('meanD')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fit{currentModel}.components{currentComponent}.stdD;
                stddev              = fit{currentModel}.components{currentComponent}.std_stdD;
                ci                  = [ fit{currentModel}.components{currentComponent}.stdD - 1.960*fit{currentModel}.components{currentComponent}.std_stdD , ...
                                        fit{currentModel}.components{currentComponent}.stdD + 1.960*fit{currentModel}.components{currentComponent}.std_stdD];
                disp('stdD')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fit{currentModel}.components{currentComponent}.spreadD;
                stddev              = fit{currentModel}.components{currentComponent}.std_spreadD;
                ci                  = [ fit{currentModel}.components{currentComponent}.spreadD - 1.960*fit{currentModel}.components{currentComponent}.std_spreadD , ...
                                        fit{currentModel}.components{currentComponent}.spreadD + 1.960*fit{currentModel}.components{currentComponent}.std_spreadD];
                disp('spreadD (stdD/meanD, coefficient of variation)')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fit{currentModel}.components{currentComponent}.modeD;
                stddev              = fit{currentModel}.components{currentComponent}.std_modeD;
                ci                  = [ fit{currentModel}.components{currentComponent}.modeD - 1.960*fit{currentModel}.components{currentComponent}.std_modeD , ...
                                        fit{currentModel}.components{currentComponent}.modeD + 1.960*fit{currentModel}.components{currentComponent}.std_modeD];
                disp('modeD (peak of distribution)')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fit{currentModel}.components{currentComponent}.theta;
                stddev              = fit{currentModel}.components{currentComponent}.std_theta;
                ci                  = [ fit{currentModel}.components{currentComponent}.theta - 1.960*fit{currentModel}.components{currentComponent}.std_theta , ...
                                        fit{currentModel}.components{currentComponent}.theta + 1.960*fit{currentModel}.components{currentComponent}.std_theta];
                disp('theta')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            case 'gamma'
                disp('Decay model:          I(b) = (beta/(beta+b))^alpha')
                disp('')
                disp('Distribution model:   f(D) = beta^alpha/gamma(alpha)*D^(alpha-1)*exp(-beta*D)')
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
                
                value               = fit{currentModel}.components{currentComponent}.meanD;
                stddev              = fit{currentModel}.components{currentComponent}.std_meanD;
                ci                  = [ fit{currentModel}.components{currentComponent}.meanD - 1.960*fit{currentModel}.components{currentComponent}.std_meanD , ...
                                        fit{currentModel}.components{currentComponent}.meanD + 1.960*fit{currentModel}.components{currentComponent}.std_meanD];
                disp('meanD')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fit{currentModel}.components{currentComponent}.stdD;
                stddev              = fit{currentModel}.components{currentComponent}.std_stdD;
                ci                  = [ fit{currentModel}.components{currentComponent}.stdD - 1.960*fit{currentModel}.components{currentComponent}.std_stdD , ...
                                        fit{currentModel}.components{currentComponent}.stdD + 1.960*fit{currentModel}.components{currentComponent}.std_stdD];
                disp('stdD')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fit{currentModel}.components{currentComponent}.spreadD;
                stddev              = fit{currentModel}.components{currentComponent}.std_spreadD;
                ci                  = [ fit{currentModel}.components{currentComponent}.spreadD - 1.960*fit{currentModel}.components{currentComponent}.std_spreadD , ...
                                        fit{currentModel}.components{currentComponent}.spreadD + 1.960*fit{currentModel}.components{currentComponent}.std_spreadD];
                disp('spreadD (stdD/meanD, coefficient of variation)')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fit{currentModel}.components{currentComponent}.modeD;
                stddev              = fit{currentModel}.components{currentComponent}.std_modeD;
                ci                  = [ fit{currentModel}.components{currentComponent}.modeD - 1.960*fit{currentModel}.components{currentComponent}.std_modeD , ...
                                        fit{currentModel}.components{currentComponent}.modeD + 1.960*fit{currentModel}.components{currentComponent}.std_modeD];
                disp('modeD (peak of distribution)')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fit{currentModel}.components{currentComponent}.theta;
                stddev              = fit{currentModel}.components{currentComponent}.std_theta;
                ci                  = [ fit{currentModel}.components{currentComponent}.theta - 1.960*fit{currentModel}.components{currentComponent}.std_theta , ...
                                        fit{currentModel}.components{currentComponent}.theta + 1.960*fit{currentModel}.components{currentComponent}.std_theta];
                disp('theta')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
        end
    end
end

