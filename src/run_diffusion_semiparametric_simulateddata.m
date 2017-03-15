%% Initialization.

clear
clc
close all hidden

%% Path to subroutines.

addpath('diffusion_semiparametric');

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
number_of_fits = 100;

% Number of Monte Carlo repetitions (0 = no error analysis).
number_of_mc_fits = 0;

% Type of model (combine exponential, stretched exponential, lognormal, and
% gamma freely)
model = {{'gamma'}, {'exponential'}};
% model = {{'lognormal'}};

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