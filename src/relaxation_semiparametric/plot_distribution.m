function [] = plot_distribution(fit)

number_of_components = numel(fit.components);

COLORS = get(groot,'DefaultAxesColorOrder');

fig = figure();
fig.Units = 'centimeters';
fig.PaperUnits = 'centimeters';
fig.Position = [0 0 12 9];
fig.PaperPosition = [0 0 12 9];

ax1 = subplot('Position',[0.15 0.20 0.8 0.75]);
ax1.FontSize = 12;
ax1.Box = 'on';
ax1.XLabel.String = 'T (s)';
ax1.YLabel.String = 'Probability';
ax1.XScale = 'log';

minT = inf;
maxT = 0;

for current_component = 1:number_of_components
    switch fit.components{current_component}.model
        case 'exponential'
            minT = min(minT, fit.components{current_component}.T);
            maxT = max(maxT, fit.components{current_component}.T);
        case 'stretchedexponential'
            minT = min(minT, fit.components{current_component}.T);
            maxT = max(maxT, fit.components{current_component}.T);
        case 'lognormal'
            minT = min(minT, logninv(1e-6, fit.components{current_component}.mu, fit.components{current_component}.sigma));
            maxT = max(maxT, logninv(1-1e-6, fit.components{current_component}.mu, fit.components{current_component}.sigma));
        case 'inversegamma' % NB the 'P' values for the inversion are 'inverted', P -> 1 - P, evaluated for the gamma, not the inversegamma
            minT = min(minT, gaminv(1-1e-6, fit.components{current_component}.alpha, 1 / fit.components{current_component}.beta));
            maxT = max(maxT, gaminv(1e-6, fit.components{current_component}.alpha, 1 / fit.components{current_component}.beta));
    end
end

minT = 0.5 * minT;
maxT = 2.0 * maxT;

TT = logspace(log10(minT),log10(maxT),10000);
y = zeros(size(TT));

for current_component = 1:number_of_components
    switch fit.components{current_component}.model
        case 'exponential'
            hl = line(ones(1,2) * fit.components{current_component}.T,[0 fit.components{current_component}.theta]);
            hl.Color = COLORS(1, :);
        case 'stretchedexponential'
            hl = line(ones(1,2) * fit.components{current_component}.T,[0 fit.components{current_component}.theta]);
            hl.Color = COLORS(1, :);
        case 'lognormal'
            y = y + fit.components{current_component}.theta * lognpdf(TT,fit.components{current_component}.mu,fit.components{current_component}.sigma);
        case 'inversegamma'
            y = y + fit.components{current_component}.theta * inversegamma_pdf(TT, fit.components{current_component}.alpha, fit.components{current_component}.beta);
    end
end
hl = line(TT, y / max(y));
hl.Color = COLORS(1, :);

end