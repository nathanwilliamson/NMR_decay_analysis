%parametric model fitting routine for diffusion data, V 1.0
%written by Magnus Roding
%% Initiate.
clear
clc
close all hidden

addpath('subroutines');

seed    = round(sum(1e6*clock()));
s       = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(s);

%% Data, settings, parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO BE EDITED BY USER.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exp=[5,6,7,8];
lexp=length(exp);
for expi=1:lexp
% % File names.
% % filename_dx                 = 't1t2_sim.dx';
% % filename_xml                = 'diff_sim.xml';
filename_dx                = sprintf('Sample2_checkDeg/%d/pdata/1/t1t2.dx',exp(expi));
filename_xml                = sprintf('%d/diff.xml',exp(expi));
 experiment                  = readdiffdata(filename_dx,filename_xml);
experimentcell{expi,1}=experiment;
 % 
% % Select which peak to analyze.
 peakInd                     = 1;
 k                           = experiment.k{peakInd};
 I                           = experiment.I{peakInd};
% 
% % Select which points in the data set to include in fitting.
 pointsInd                   = 1:numel(k); % 
 k                           = k(pointsInd);
 I                           = I(pointsInd);
kcell{expi,1}=k;
Icell{expi,1}=I;
end
%load('data_ps_alt.mat')
%k =k';
%I=I';
% return

% Number of fittings (with random initializations).
nFits                       = 10;

% Number of Monte Carlo repetitions (0 = no error analysis).
nMC                         = 0;

% Type of model (combine exponential, stretched exponential, lognormal, and
% gamma freely). Specify more than one model for joint analysis and 
% comparison.
% models = {    {{'lognormal','mu',-24},{'exponential'}} };
models = {    {{'exponential'}}      };
%models = { {{'exponential'}} };
%models = {    {{'lognormal'},{'exponential','D',[1.6e-9 1.8e-9]}} };
%models = {    {{'exponential'},{'exponential','D',[1e-9 2e-9]}} };
%models = {    {{'gamma'}}};
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
fitcell=cell(lexp,1);
for expi=1:lexp
fit                        = cell(nModels,1);
for currentModel = 1:nModels
    fit_                        = analyze(kcell{expi,1},Icell{expi,1},models{currentModel},baseline,nFits,nMC);
    fit{currentModel}           = fit_;
end
fitcell{expi,1}=fit;
end
%% Plot fits and residuals.
close all
for expi=1:lexp
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
ax1.YLim                    = [max(1e-4,0.5*min(Icell{expi,1})) 1.25*max(Icell{expi,1})];

hl                          = line(kcell{expi,1},Icell{expi,1});
hl.Marker                   = 'o';
hl.Color                    = [0 0 0];
hl.LineStyle                = 'none';

legendStr                   = cell(nModels,1);
legendEntries               = [];
for currentModel = 1:nModels
    hl                          = line(kcell{expi,1},fitcell{expi,1}{currentModel}.Imodel);
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
title(sprintf('DELTA = %d ms' ,experimentcell{expi}.DELTA * 1000))
legend(legendEntries,legendStr)

ax2                         = subplot('Position',[0.15 0.20 0.8 0.15]);
ax2.FontSize                = 12;
ax2.Box                     = 'on';
ax2.XLabel.String           = 'k (sm^{-2})';
ax2.YLabel.String           = 'Residual';
ax2.YLabel.Units            = 'centimeters';
ax2.YLabel.Position(1)      = ax1.YLabel.Position(1);

hl                          = line(kcell{expi,1},zeros(size(kcell{expi,1})));
hl.Color                    = [0 0 0];

for currentModel = 1:nModels
    hl                          = line(kcell{expi,1},fitcell{expi,1}{currentModel}.residuals);
    hl.Color                    = COLORS(currentModel,:);
    hl.Marker                   = 'o';
    hl.LineStyle                = 'none';
end
end
%% Plot distributions of diffusion coefficients.
for expi=1:lexp
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
                minD = min(minD,fitcell{expi,1}{currentModel}.components{currentComponent}.D);
                maxD = max(maxD,fitcell{expi,1}{currentModel}.components{currentComponent}.D);
            case 'stretchedexponential'
                minD = min(minD,fitcell{expi,1}{currentModel}.components{currentComponent}.D);
                maxD = max(maxD,fitcell{expi,1}{currentModel}.components{currentComponent}.D);
            case 'lognormal'
                minD = min(minD, logninv(1e-6,fitcell{expi,1}{currentModel}.components{currentComponent}.mu,fitcell{expi,1}{currentModel}.components{currentComponent}.sigma));
                maxD = max(maxD, logninv(1-1e-6,fitcell{expi,1}{currentModel}.components{currentComponent}.mu,fitcell{expi,1}{currentModel}.components{currentComponent}.sigma));
            case 'gamma'
                minD = min(minD, gaminv(1e-6,fitcell{expi,1}{currentModel}.components{currentComponent}.alpha,1/fitcell{expi,1}{currentModel}.components{currentComponent}.beta));
                maxD = max(maxD, gaminv(1-1e-6,fitcell{expi,1}{currentModel}.components{currentComponent}.alpha,1/fitcell{expi,1}{currentModel}.components{currentComponent}.beta));
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
                hl = line(ones(1,2)*fitcell{expi,1}{currentModel}.components{currentComponent}.D,[0 fitcell{expi,1}{currentModel}.components{currentComponent}.theta]);
                hl.Color                    = COLORS(currentModel,:);
            case 'stretchedexponential'
                hl = line(ones(1,2)*fitcell{expi,1}{currentModel}.components{currentComponent}.D,[0 fitcell{expi,1}{currentModel}.components{currentComponent}.theta]);
                hl.Color                    = COLORS(currentModel,:);
            case 'lognormal'
                y = y + fitcell{expi,1}{currentModel}.components{currentComponent}.theta * lognpdf(DD,fitcell{expi,1}{currentModel}.components{currentComponent}.mu,fitcell{expi,1}{currentModel}.components{currentComponent}.sigma);
            case 'gamma'
                y = y + fitcell{expi,1}{currentModel}.components{currentComponent}.theta * gampdf(DD,fitcell{expi,1}{currentModel}.components{currentComponent}.alpha,1/fitcell{expi,1}{currentModel}.components{currentComponent}.beta);
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
title(sprintf('DELTA = %d ms' ,experimentcell{expi}.DELTA * 1000))
legend(legendEntries,legendStr)
end
%% Plot histograms of Monte Carlo values for all parameters to check manually for Gaussianity.

if plotMChist
    for currentModel = 1:nModels
        for i = 1:size(fitcell{expi,1}{currentModel}.paramhat_MC,2)
            figure, hist(fitcell{expi,1}{currentModel}.paramhat_MC(:,i))
        end
    end
end
%% Plot attenuation data overlay
minY=min(Icell{1} / fitcell{1}{1}.I0);
maxY=max(Icell{1} / fitcell{1}{1}.I0);

currentModel=1;

for expi=1:lexp
    if min(Icell{expi} / fitcell{expi}{currentModel}.I0) < minY
        minY = min(Icell{expi} / fitcell{expi}{currentModel}.I0);
    end
    if max(Icell{expi} / fitcell{expi}{currentModel}.I0)>maxY
        maxY = max(Icell{expi} / fitcell{expi}{currentModel}.I0);
    end
end

COLORS                      = get(groot,'DefaultAxesColorOrder');

fig                         = figure();
fig.Units                   = 'centimeters';
fig.PaperUnits              = 'centimeters';
fig.Position                = [0 0 12 9];
fig.PaperPosition           = [0 0 12 9];


ax1                         = subplot('Position',[0.15 0.15 0.8 0.75]);
ax1.FontSize                = 12;
ax1.Box                     = 'on';
ax1.YLabel.String           = 'Intensity';
ax1.YScale                  = 'log';
%ax1.XTickLabel              = {};
ax1.YLabel.Units            = 'centimeters';
ax1.YLim                    = [max(1e-4,0.5*minY) 1.25*maxY];
ax1.XLabel.String           = 'b (sm^{-2})';

legendStr                   = cell(lexp,1);
legendEntries               = [];

for expi=1:lexp
hl                          = line(kcell{expi},Icell{expi}/fitcell{expi}{currentModel}.I0);
hl.Marker                   = 'o';
hl.Color                    = COLORS(expi,:);

hl.LineStyle                = 'none';




%h2                          = line(kcell{expi},fitcell{expi}{currentModel}.model/fitcell{expi}{currentModel}.I0);
%h2.Color                    = COLORS(expi,:);
    
legendStr{expi}             = sprintf('%d ms' ,experimentcell{expi}.DELTA * 1000);
legendEntries(expi) = hl;

end
legend(legendEntries,legendStr)
%% Present report of result in text form.
for expi = 1:lexp
    disp('===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ==');
    disp((sprintf('DELTA = %d ms' ,experimentcell{expi}.DELTA * 1000)))
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
    disp(['Fit (residual sum of squares):   ' num2str(fitcell{expi,1}{currentModel}.ss)])
    disp('')
    disp(['I0:                              ' num2str(fitcell{expi,1}{currentModel}.I0)])
    disp('')
    disp(['baseline:                        ' num2str(fitcell{expi,1}{currentModel}.baseline)])
    disp('')
    for currentComponent = 1:numel(fitcell{expi,1}{currentModel}.components)
        disp('--------------------------------------------------------------');
        disp(['Component ' num2str(currentComponent) ': ' fitcell{expi,1}{currentModel}.components{currentComponent}.model])
        disp('')

        switch fitcell{expi,1}{currentModel}.components{currentComponent}.model
            case 'exponential'
                disp('Decay model:          I(k) = exp(-k*D)')
                disp('')
                disp('Distribution model:   f(D) = delta(D-D0)')
                disp('--------------------------------------------------------------')
                
                value               = fitcell{expi,1}{currentModel}.components{currentComponent}.D;
                stddev              = fitcell{expi,1}{currentModel}.components{currentComponent}.std_D;
                ci                  = [ fitcell{expi,1}{currentModel}.components{currentComponent}.D - 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_D , ...
                                        fitcell{expi,1}{currentModel}.components{currentComponent}.D + 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_D];
                disp('D')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fitcell{expi,1}{currentModel}.components{currentComponent}.theta;
                stddev              = fitcell{expi,1}{currentModel}.components{currentComponent}.std_theta;
                ci                  = [ fitcell{expi,1}{currentModel}.components{currentComponent}.theta - 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_theta , ...
                                        fitcell{expi,1}{currentModel}.components{currentComponent}.theta + 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_theta];
                disp('theta')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            case 'stretchedexponential'
                disp('Decay model:          I(k) = exp(-(k*D)^beta)')
                disp('')
                disp('Distribution model:   f(D) unknown')
                disp('--------------------------------------------------------------')
                
                value               = fitcell{expi,1}{currentModel}.components{currentComponent}.D;
                stddev              = fitcell{expi,1}{currentModel}.components{currentComponent}.std_D;
                ci                  = [ fitcell{expi,1}{currentModel}.components{currentComponent}.D - 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_D , ...
                                        fitcell{expi,1}{currentModel}.components{currentComponent}.D + 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_D];
                disp('D')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fitcell{expi,1}{currentModel}.components{currentComponent}.beta;
                stddev              = fitcell{expi,1}{currentModel}.components{currentComponent}.std_beta;
                ci                  = [ fitcell{expi,1}{currentModel}.components{currentComponent}.beta - 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_beta , ...
                                        fitcell{expi,1}{currentModel}.components{currentComponent}.beta + 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_beta];
                disp('beta')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fitcell{expi,1}{currentModel}.components{currentComponent}.theta;
                stddev              = fitcell{expi,1}{currentModel}.components{currentComponent}.std_theta;
                ci                  = [ fitcell{expi,1}{currentModel}.components{currentComponent}.theta - 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_theta , ...
                                        fitcell{expi,1}{currentModel}.components{currentComponent}.theta + 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_theta];
                disp('theta')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            case 'lognormal'
                disp('Decay model:          I(k) = integral(f(D),0,inf) (numerical approx.)')
                disp('')
                disp('Distribution model:   f(D) = 1/(D*sigma*sqrt(2*pi))*exp(-(log(D)-mu)^2/(2*sigma^2))')
                disp('--------------------------------------------------------------')
                
                value               = fitcell{expi,1}{currentModel}.components{currentComponent}.mu;
                stddev              = fitcell{expi,1}{currentModel}.components{currentComponent}.std_mu;
                ci                  = [ fitcell{expi,1}{currentModel}.components{currentComponent}.mu - 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_mu , ...
                                        fitcell{expi,1}{currentModel}.components{currentComponent}.mu + 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_mu];
                disp('mu')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fitcell{expi,1}{currentModel}.components{currentComponent}.sigma;
                stddev              = fitcell{expi,1}{currentModel}.components{currentComponent}.std_sigma;
                ci                  = [ fitcell{expi,1}{currentModel}.components{currentComponent}.sigma - 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_sigma , ...
                                        fitcell{expi,1}{currentModel}.components{currentComponent}.sigma + 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_sigma];
                disp('sigma')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fitcell{expi,1}{currentModel}.components{currentComponent}.meanD;
                stddev              = fitcell{expi,1}{currentModel}.components{currentComponent}.std_meanD;
                ci                  = [ fitcell{expi,1}{currentModel}.components{currentComponent}.meanD - 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_meanD , ...
                                        fitcell{expi,1}{currentModel}.components{currentComponent}.meanD + 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_meanD];
                disp('meanD')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fitcell{expi,1}{currentModel}.components{currentComponent}.stdD;
                stddev              = fitcell{expi,1}{currentModel}.components{currentComponent}.std_stdD;
                ci                  = [ fitcell{expi,1}{currentModel}.components{currentComponent}.stdD - 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_stdD , ...
                                        fitcell{expi,1}{currentModel}.components{currentComponent}.stdD + 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_stdD];
                disp('stdD')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fitcell{expi,1}{currentModel}.components{currentComponent}.spreadD;
                stddev              = fitcell{expi,1}{currentModel}.components{currentComponent}.std_spreadD;
                ci                  = [ fitcell{expi,1}{currentModel}.components{currentComponent}.spreadD - 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_spreadD , ...
                                        fitcell{expi,1}{currentModel}.components{currentComponent}.spreadD + 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_spreadD];
                disp('spreadD (stdD/meanD, coefficient of variation)')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fitcell{expi,1}{currentModel}.components{currentComponent}.modeD;
                stddev              = fitcell{expi,1}{currentModel}.components{currentComponent}.std_modeD;
                ci                  = [ fitcell{expi,1}{currentModel}.components{currentComponent}.modeD - 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_modeD , ...
                                        fitcell{expi,1}{currentModel}.components{currentComponent}.modeD + 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_modeD];
                disp('modeD (peak of distribution)')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fitcell{expi,1}{currentModel}.components{currentComponent}.theta;
                stddev              = fitcell{expi,1}{currentModel}.components{currentComponent}.std_theta;
                ci                  = [ fitcell{expi,1}{currentModel}.components{currentComponent}.theta - 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_theta , ...
                                        fitcell{expi,1}{currentModel}.components{currentComponent}.theta + 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_theta];
                disp('theta')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
            case 'gamma'
                disp('Decay model:          I(k) = (beta/(beta+k))^alpha')
                disp('')
                disp('Distribution model:   f(D) = beta^alpha/gamma(alpha)*D^(alpha-1)*exp(-beta*D)')
                disp('--------------------------------------------------------------')
                
                value               = fitcell{expi,1}{currentModel}.components{currentComponent}.alpha;
                stddev              = fitcell{expi,1}{currentModel}.components{currentComponent}.std_alpha;
                ci                  = [ fitcell{expi,1}{currentModel}.components{currentComponent}.alpha - 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_alpha , ...
                                        fitcell{expi,1}{currentModel}.components{currentComponent}.alpha + 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_alpha];
                disp('alpha')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fitcell{expi,1}{currentModel}.components{currentComponent}.beta;
                stddev              = fitcell{expi,1}{currentModel}.components{currentComponent}.std_beta;
                ci                  = [ fitcell{expi,1}{currentModel}.components{currentComponent}.beta - 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_beta , ...
                                        fitcell{expi,1}{currentModel}.components{currentComponent}.beta + 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_beta];
                disp('beta')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fitcell{expi,1}{currentModel}.components{currentComponent}.meanD;
                stddev              = fitcell{expi,1}{currentModel}.components{currentComponent}.std_meanD;
                ci                  = [ fitcell{expi,1}{currentModel}.components{currentComponent}.meanD - 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_meanD , ...
                                        fitcell{expi,1}{currentModel}.components{currentComponent}.meanD + 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_meanD];
                disp('meanD')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fitcell{expi,1}{currentModel}.components{currentComponent}.stdD;
                stddev              = fitcell{expi,1}{currentModel}.components{currentComponent}.std_stdD;
                ci                  = [ fitcell{expi,1}{currentModel}.components{currentComponent}.stdD - 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_stdD , ...
                                        fitcell{expi,1}{currentModel}.components{currentComponent}.stdD + 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_stdD];
                disp('stdD')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fitcell{expi,1}{currentModel}.components{currentComponent}.spreadD;
                stddev              = fitcell{expi,1}{currentModel}.components{currentComponent}.std_spreadD;
                ci                  = [ fitcell{expi,1}{currentModel}.components{currentComponent}.spreadD - 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_spreadD , ...
                                        fitcell{expi,1}{currentModel}.components{currentComponent}.spreadD + 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_spreadD];
                disp('spreadD (stdD/meanD, coefficient of variation)')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fitcell{expi,1}{currentModel}.components{currentComponent}.modeD;
                stddev              = fitcell{expi,1}{currentModel}.components{currentComponent}.std_modeD;
                ci                  = [ fitcell{expi,1}{currentModel}.components{currentComponent}.modeD - 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_modeD , ...
                                        fitcell{expi,1}{currentModel}.components{currentComponent}.modeD + 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_modeD];
                disp('modeD (peak of distribution)')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
                
                value               = fitcell{expi,1}{currentModel}.components{currentComponent}.theta;
                stddev              = fitcell{expi,1}{currentModel}.components{currentComponent}.std_theta;
                ci                  = [ fitcell{expi,1}{currentModel}.components{currentComponent}.theta - 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_theta , ...
                                        fitcell{expi,1}{currentModel}.components{currentComponent}.theta + 1.960*fitcell{expi,1}{currentModel}.components{currentComponent}.std_theta];
                disp('theta')
                disp(['      - Value:                   ' num2str(value)])
                disp(['      - Std:                     ' num2str(stddev)])
                disp(['      - 95 % CI:                 ' '[' num2str(ci(1)) ', ' num2str(ci(2)) ']'])
        end
    end
end
end
