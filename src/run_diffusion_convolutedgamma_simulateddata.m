%% Initialization.

clear
clc
close all hidden

%% Path to subroutines.

addpath('diffusion_convolutedgamma');

%% Random stream.

random_seed = round( sum( 1e6 * clock() ) );
s = RandStream('mt19937ar', 'Seed', random_seed);
RandStream.setGlobalStream(s);

%% Simulate data set.

b = linspace(0, 1e2, 64);
b = b(:);

I = signal(b, {{'lognormal'}}, [-0.75, 0.5, 1e-2, 1, 1]);

sigma_error = 0.01;
I = I + sigma_error * randn(size(I));

%% Model and fit parameters.

% Number of fits (with random initializations).
number_of_fits = 3;

% Number of Monte Carlo repetitions (0 = no error analysis).
number_of_mc_fits = 0;

% Type of model (combine exponential, stretched exponential, lognormal, and
% gamma freely)
% model = {{'gamma'}, {'exponential'}};
% model = {{'lognormal'}};
model = {{'gamma'}};

% Baseline toggle.
baseline = true;

%% Fit model and estimate parameters.

fit = analyze(b, I, model, baseline, number_of_fits, number_of_mc_fits);

%% Print results.

print_results(fit);

%% Plot fit and residuals.

plot_fit_and_residuals(b, I, fit);

%% Plot distribution of D.

plot_distribution(fit);


%Gamma convolution model for diffusion data
%written by Magnus Roding
%published in Roding, Williamson, Nyden, JMR 261 (2015)
%% Initiate.

clear
clc
close all hidden

disp('Hello!')

seed = round(sum(1e6*clock()));
s = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(s);

%% 1st option for data preparation

% b=load('SampleData/cartilage/b_cartilage.txt');
% I=load('SampleData/cartilage/I_cartilage.txt');

%% second option for data preparation
% File names.

filename_dx                 = 'data_pva/t1t2.dx';
filename_xml                = 'data_pva/diff.xml';
experiment                  = readdiffdatabruker(filename_dx,filename_xml);
peakInd                     = 3;
b                           = experiment.b{peakInd};
I                           = experiment.I{peakInd};
b                           = b(:)';
I                           = I(:)';
I=I/I(1)

clear filename* experiment peakInd

 
% Just add this to get a positive baseline value, otherwise the logscale 
% plot looks terrible. Does of course not affect any other parameter 
% estimates than that of the baseline itself.
I = I + 0.003;

%% Parameters.

baseline = true; % true = include baseline in models.
nFits = 10; % Number of random initial guesses.

%% Fit gamma model (equivalent to the convoluted gamma model for n = 1).
% This simpler code for only the gamma model is included for clarity.

[muhat_gamma,sigmahat_gamma,alphahat_gamma,betahat_gamma,I0hat_gamma,Ibhat_gamma,ss_gamma] = analyze_gamma(b,I,baseline,nFits);

%% Fit convoluted gamma model.

nConv = 2; %number of convolutions 
[muhat_gamma2,sigmahat_gamma2,alphahat_gamma2,betahat_gamma2,I0hat_gamma2,Ibhat_gamma2,ss_gamma2] = analyze_convolutedgamma(b,I,nConv,baseline,nFits);

%% Fit lognormal model.

[muhat_lognormal,sigmahat_lognormal,I0hat_lognormal,Ibhat_lognormal,ss_lognormal] = analyze_lognormal(b,I,baseline,nFits);

%% Analyze and present results.

mean_gamma2 = sum(muhat_gamma2);
var_gamma2 = sum(sigmahat_gamma2.^2);
std_gamma2 = sqrt(var_gamma2);
spread_gamma2 = std_gamma2/mean_gamma2;

mean_lognormal = exp(muhat_lognormal+1/2*sigmahat_lognormal^2);
var_lognormal = (exp(sigmahat_lognormal^2)-1)*exp(2*muhat_lognormal+sigmahat_lognormal^2);
std_lognormal = sqrt(var_lognormal);
spread_lognormal = std_lognormal/mean_lognormal;

disp('Results:')
disp('---------------------------')
disp('Gamma model')
disp(['   Sum of squares: ' num2str(ss_gamma)]);
disp(['   Mean of D:      ' num2str(muhat_gamma)]);
disp(['   Std of D:       ' num2str(sigmahat_gamma)]);
disp(['   Spread of D:    ' num2str(sigmahat_gamma/muhat_gamma)]);
disp('Convoluted gamma model with n = 2')
disp(['   Sum of squares: ' num2str(ss_gamma2)]);
disp(['   Mean of D:      ' num2str(mean_gamma2)]);
disp(['   Std of D:       ' num2str(std_gamma2)]);
disp(['   Spread of D:    ' num2str(spread_gamma2)]);
disp('Lognormal model')
disp(['   Sum of squares: ' num2str(ss_lognormal)]);
disp(['   Mean of D:      ' num2str(mean_lognormal)]);
disp(['   Std of D:       ' num2str(std_lognormal)]);
disp(['   Spread of D:    ' num2str(spread_lognormal)]);

% Plot figure.

COLORS = get(groot,'DefaultAxesColorOrder');

fig = figure();
fig.Units = 'centimeters';
fig.PaperUnits = 'centimeters';
fig.Position = [0 0 8 10.5];
fig.PaperPosition = fig.Position;
FontName = 'helvetica';
FontSize = 7;
FontWeight = 'normal'; 

h = axes();
h.Units = 'centimeters';
h.FontName = FontName;
h.FontSize = FontSize;
h.FontWeight = FontWeight;
h.Position = [1.25 7.25 6.5 3];
h.Box = 'on';
h.XTickLabel = {};
h.YLabel.String = 'Intensity';
h.YScale = 'log';
h.YLabel.Units = 'centimeters';
h.YLabel.Position(1) = -0.65;
hold on


I_gamma = signal_gamma(b,muhat_gamma,sigmahat_gamma,I0hat_gamma,Ibhat_gamma);
hl = line(b/1e11,I_gamma);
hl.Color = COLORS(1,:);

I_gamma2 = signal_convolutedgamma(b,muhat_gamma2,sigmahat_gamma2,I0hat_gamma2,Ibhat_gamma2);
hl = line(b/1e11,I_gamma2);
hl.Color = COLORS(2,:);

I_lognormal = signal_lognormal(b,muhat_lognormal,sigmahat_lognormal,I0hat_lognormal,Ibhat_lognormal);
hl = line(b/1e11,I_lognormal);
hl.Color = COLORS(3,:);

legend('Gamma','Gamma convolution', 'Lognormal');

hl = line(b/1e11,I);
hl.LineStyle = 'none';
hl.Marker = 'o';
hl.MarkerSize = 3;
hl.Color = [0 0 0];

I_gamma = signal_gamma(b,muhat_gamma,sigmahat_gamma,I0hat_gamma,Ibhat_gamma);
hl = line(b/1e11,I_gamma);
hl.Color = COLORS(1,:);

I_gamma2 = signal_convolutedgamma(b,muhat_gamma2,sigmahat_gamma2,I0hat_gamma2,Ibhat_gamma2);
hl = line(b/1e11,I_gamma2);
hl.Color = COLORS(2,:);

I_lognormal = signal_lognormal(b,muhat_lognormal,sigmahat_lognormal,I0hat_lognormal,Ibhat_lognormal);
hl = line(b/1e11,I_lognormal);
hl.Color = COLORS(3,:);

h = axes();
h.Units = 'centimeters';
h.FontName = FontName;
h.FontSize = FontSize;
h.FontWeight = FontWeight;
h.Position = [1.25 5 6.5 2];
h.Box = 'on';
h.XLabel.String = 'b (\times 10^{11} s/m^2)';
h.YLabel.String = 'Residual (\times 10^{-3})';
h.YLim = [-2.5 2.5];
h.YLabel.Units = 'centimeters';
h.YLabel.Position(1) = -0.65;
hold on

hl = line(b/1e11,zeros(size(b)));
hl.Color = [0 0 0];

hl = line(b/1e11,(I_gamma-I)/1e-3);
hl.LineStyle = 'none';
hl.Marker = 'o';
hl.MarkerSize = 3;
hl.Color = COLORS(1,:);

hl = line(b/1e11,(I_gamma2-I)/1e-3);
hl.LineStyle = 'none';
hl.Marker = 'o';
hl.MarkerSize = 3;
hl.Color = COLORS(2,:);

hl = line(b/1e11,(I_lognormal-I)/1e-3);
hl.LineStyle = 'none';
hl.Marker           = 'o';
hl.MarkerSize       = 3;
hl.Color = COLORS(3,:);

h = axes();
h.Units = 'centimeters';
h.FontName = FontName;
h.FontSize = FontSize;
h.FontWeight = FontWeight;
h.Position = [1.25 1 6.5 3];
h.Box = 'on';
h.XLabel.String = 'D (m^2/s)';
h.YLabel.String = 'Probability';
h.XScale = 'log';
h.YTickLabel = {};
h.YLabel.Units = 'centimeters';
h.YLabel.Position(1) = -0.65;
hold on

D = logspace(-12,-10,10000);
dD = diff(logspace(-12,-10,10001));

bmax = max(b);

f_gamma = gampdf(D*bmax,alphahat_gamma,1/(betahat_gamma/bmax));
hl = line(D,f_gamma);
hl.Color = COLORS(1,:);

nTerms = 1e4;
f_gamma2 = pdf_convolutedgamma(D*bmax,alphahat_gamma2,betahat_gamma2/bmax,nTerms);
hl = line(D,f_gamma2);
hl.Color = COLORS(2,:);

f_lognormal = lognpdf(D*bmax,muhat_lognormal+log(bmax),sigmahat_lognormal);
hl = line(D,f_lognormal);
hl.Color = COLORS(3,:);