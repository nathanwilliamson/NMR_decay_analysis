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
%% 1st option for data preparation
% b=load('SampleData/cartilage/b_cartilage.txt');
% I=load('SampleData/cartilage/I_cartilage.txt');

%% 2nd option for data preparation
% File names.

filename_dx                = 'SampleData/PVA/t1t2.dx';
filename_xml                = 'SampleData/PVA/diff.xml';
experiment                  = readdiffdata(filename_dx,filename_xml);

% Select which peak to analyze.

peakInd                     = 1;
 
b                           = experiment.b{peakInd};
I                           = experiment.I{peakInd};

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