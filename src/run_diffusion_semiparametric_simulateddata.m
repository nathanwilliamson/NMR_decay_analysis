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

b = linspace(0, 1e2, 64);
I = signal(b, {{'gamma'}}, [2, 2, 0, 1, 1]);

b = b(:);
I = I(:);

%% Model and fit parameters.

% Number of fits (with random initializations).
number_of_fits = 10;

% Number of Monte Carlo repetitions (0 = no error analysis).
number_of_mc_fits = 10;

% Type of model (combine exponential, stretched exponential, lognormal, and
% gamma freely)
model = {{'lognormal'}};

% Baseline toggle.
baseline = false;

%% Fit model and estimate parameters.

fit = analyze(b, I, model, baseline, number_of_fits, number_of_mc_fits);