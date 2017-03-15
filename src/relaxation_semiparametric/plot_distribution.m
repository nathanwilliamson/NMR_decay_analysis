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
ax1.XLabel.String = 'D (m^2s^{-1})';
ax1.YLabel.String = 'Probability';
ax1.XScale = 'log';

minD = inf;
maxD = 0;

for current_component = 1:number_of_components
    switch fit.components{current_component}.model
        case 'exponential'
            minD = min(minD, fit.components{current_component}.D);
            maxD = max(maxD, fit.components{current_component}.D);
        case 'stretchedexponential'
            minD = min(minD, fit.components{current_component}.D);
            maxD = max(maxD, fit.components{current_component}.D);
        case 'lognormal'
            minD = min(minD, logninv(1e-6, fit.components{current_component}.mu, fit.components{current_component}.sigma));
            maxD = max(maxD, logninv(1-1e-6, fit.components{current_component}.mu, fit.components{current_component}.sigma));
        case 'gamma'
            minD = min(minD, gaminv(1e-6, fit.components{current_component}.alpha, 1 / fit.components{current_component}.beta));
            maxD = max(maxD, gaminv(1-1e-6, fit.components{current_component}.alpha, 1 / fit.components{current_component}.beta));
    end
end

minD = 0.5 * minD;
maxD = 2.0 * maxD;

DD = logspace(log10(minD),log10(maxD),10000);
y = zeros(size(DD));

for current_component = 1:number_of_components
    switch fit.components{current_component}.model
        case 'exponential'
            hl = line(ones(1,2) * fit.components{current_component}.D,[0 fit.components{current_component}.theta]);
            hl.Color = COLORS(1, :);
        case 'stretchedexponential'
            hl = line(ones(1,2) * fit.components{current_component}.D,[0 fit.components{current_component}.theta]);
            hl.Color = COLORS(1, :);
        case 'lognormal'
            y = y + fit.components{current_component}.theta * lognpdf(DD,fit.components{current_component}.mu,fit.components{current_component}.sigma);
        case 'gamma'
            y = y + fit.components{current_component}.theta * gampdf(DD,fit.components{current_component}.alpha,1/fit.components{current_component}.beta);
    end
end
hl = line(DD, y / max(y));
hl.Color = COLORS(1, :);

end