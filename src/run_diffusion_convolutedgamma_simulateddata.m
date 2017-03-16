%% Initialization.

clear
clc
close all hidden

%% Path to subroutines.

addpath('diffusion_convolutedgamma');

%% Random stream.

random_seed = round( sum( 1e6 * clock() ) );
random_stream = RandStream('mt19937ar', 'Seed', random_seed);
RandStream.setGlobalStream(random_stream);

%% Simulate data set.

b = linspace(0, 1e2, 64);

mu = [0.35, 0.15];
sigma = [0.2, 0.1];
I0 = 1;
Ib = 0.001;

I = signal(b, mu, sigma, I0, Ib);

sigma_error = 0.00;
I = I + sigma_error * randn(size(I));

% semilogy(b, I);
% return

%% Model and fit parameters.

% Number of gamma convolutions.
number_of_convolutions = 2;

% Number of fits (with random initializations).
number_of_fits = 100;

% Baseline toggle.
baseline = true;

%% Fit model and estimate parameters.

[muhat, sigmahat, alphahat, betahat, I0hat, Ibhat, ss] = analyze(b, I, number_of_convolutions, baseline, number_of_fits);

%% Print results.

disp('Results of fit')
disp(['   Sum of squares: ' num2str(s)]);
disp(['   Mean of D:      ' num2str(mean_gamma2)]);
disp(['   Std of D:       ' num2str(std_gamma2)]);
disp(['   Spread of D:    ' num2str(spread_gamma2)]);

%% Plot fit and residuals.

plot_fit_and_residuals(b, I, muhat, sigmahat, alphahat, betahat, I0hat, Ibhat, ss);

%% Plot distribution of D.

plot_distribution(alphahat, betahat);
% 
% disp('Results:')
% disp('---------------------------')
% disp('Gamma model')
% disp(['   Sum of squares: ' num2str(ss_gamma)]);
% disp(['   Mean of D:      ' num2str(muhat_gamma)]);
% disp(['   Std of D:       ' num2str(sigmahat_gamma)]);
% disp(['   Spread of D:    ' num2str(sigmahat_gamma/muhat_gamma)]);

% disp('Lognormal model')
% disp(['   Sum of squares: ' num2str(ss_lognormal)]);
% disp(['   Mean of D:      ' num2str(mean_lognormal)]);
% disp(['   Std of D:       ' num2str(std_lognormal)]);
% disp(['   Spread of D:    ' num2str(spread_lognormal)]);
% 
% % Plot figure.
% 
% COLORS = get(groot,'DefaultAxesColorOrder');
% 
% fig = figure();
% fig.Units = 'centimeters';
% fig.PaperUnits = 'centimeters';
% fig.Position = [0 0 8 10.5];
% fig.PaperPosition = fig.Position;
% FontName = 'helvetica';
% FontSize = 7;
% FontWeight = 'normal'; 
% 
% h = axes();
% h.Units = 'centimeters';
% h.FontName = FontName;
% h.FontSize = FontSize;
% h.FontWeight = FontWeight;
% h.Position = [1.25 7.25 6.5 3];
% h.Box = 'on';
% h.XTickLabel = {};
% h.YLabel.String = 'Intensity';
% h.YScale = 'log';
% h.YLabel.Units = 'centimeters';
% h.YLabel.Position(1) = -0.65;
% hold on
% 
% 
% I_gamma = signal_gamma(b,muhat_gamma,sigmahat_gamma,I0hat_gamma,Ibhat_gamma);
% hl = line(b/1e11,I_gamma);
% hl.Color = COLORS(1,:);
% 
% I_gamma2 = signal_convolutedgamma(b,muhat_gamma2,sigmahat_gamma2,I0hat_gamma2,Ibhat_gamma2);
% hl = line(b/1e11,I_gamma2);
% hl.Color = COLORS(2,:);
% 
% I_lognormal = signal_lognormal(b,muhat_lognormal,sigmahat_lognormal,I0hat_lognormal,Ibhat_lognormal);
% hl = line(b/1e11,I_lognormal);
% hl.Color = COLORS(3,:);
% 
% legend('Gamma','Gamma convolution', 'Lognormal');
% 
% hl = line(b/1e11,I);
% hl.LineStyle = 'none';
% hl.Marker = 'o';
% hl.MarkerSize = 3;
% hl.Color = [0 0 0];
% 
% I_gamma = signal_gamma(b,muhat_gamma,sigmahat_gamma,I0hat_gamma,Ibhat_gamma);
% hl = line(b/1e11,I_gamma);
% hl.Color = COLORS(1,:);
% 
% I_gamma2 = signal_convolutedgamma(b,muhat_gamma2,sigmahat_gamma2,I0hat_gamma2,Ibhat_gamma2);
% hl = line(b/1e11,I_gamma2);
% hl.Color = COLORS(2,:);
% 
% I_lognormal = signal_lognormal(b,muhat_lognormal,sigmahat_lognormal,I0hat_lognormal,Ibhat_lognormal);
% hl = line(b/1e11,I_lognormal);
% hl.Color = COLORS(3,:);
% 
% h = axes();
% h.Units = 'centimeters';
% h.FontName = FontName;
% h.FontSize = FontSize;
% h.FontWeight = FontWeight;
% h.Position = [1.25 5 6.5 2];
% h.Box = 'on';
% h.XLabel.String = 'b (\times 10^{11} s/m^2)';
% h.YLabel.String = 'Residual (\times 10^{-3})';
% h.YLim = [-2.5 2.5];
% h.YLabel.Units = 'centimeters';
% h.YLabel.Position(1) = -0.65;
% hold on
% 
% hl = line(b/1e11,zeros(size(b)));
% hl.Color = [0 0 0];
% 
% hl = line(b/1e11,(I_gamma-I)/1e-3);
% hl.LineStyle = 'none';
% hl.Marker = 'o';
% hl.MarkerSize = 3;
% hl.Color = COLORS(1,:);
% 
% hl = line(b/1e11,(I_gamma2-I)/1e-3);
% hl.LineStyle = 'none';
% hl.Marker = 'o';
% hl.MarkerSize = 3;
% hl.Color = COLORS(2,:);
% 
% hl = line(b/1e11,(I_lognormal-I)/1e-3);
% hl.LineStyle = 'none';
% hl.Marker           = 'o';
% hl.MarkerSize       = 3;
% hl.Color = COLORS(3,:);
% 
% h = axes();
% h.Units = 'centimeters';
% h.FontName = FontName;
% h.FontSize = FontSize;
% h.FontWeight = FontWeight;
% h.Position = [1.25 1 6.5 3];
% h.Box = 'on';
% h.XLabel.String = 'D (m^2/s)';
% h.YLabel.String = 'Probability';
% h.XScale = 'log';
% h.YTickLabel = {};
% h.YLabel.Units = 'centimeters';
% h.YLabel.Position(1) = -0.65;
% hold on
% 
% D = logspace(-12,-10,10000);
% dD = diff(logspace(-12,-10,10001));
% 
% bmax = max(b);
% 
% f_gamma = gampdf(D*bmax,alphahat_gamma,1/(betahat_gamma/bmax));
% hl = line(D,f_gamma);
% hl.Color = COLORS(1,:);
% 
% nTerms = 1e4;
% f_gamma2 = pdf_convolutedgamma(D*bmax,alphahat_gamma2,betahat_gamma2/bmax,nTerms);
% hl = line(D,f_gamma2);
% hl.Color = COLORS(2,:);
% 
% f_lognormal = lognpdf(D*bmax,muhat_lognormal+log(bmax),sigmahat_lognormal);
% hl = line(D,f_lognormal);
% hl.Color = COLORS(3,:);