function [] = plot_fits_and_residuals(fit)

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
    hl                          = line(b,fit.Imodel);
    hl.Color                    = COLORS(currentModel,:);
    
    if numel(models) == 1
        legendStr     = models{1}{1};
    else
        legendStr     = [];
        for i = 1:numel(models)
            legendStr = [legendStr models{i}{1} '+'];
        end
        legendStr = legendStr(1:end-1);
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
    hl                          = line(b,fit.residuals);
    hl.Color                    = COLORS(currentModel,:);
    hl.Marker                   = 'o';
    hl.LineStyle                = 'none';
end



end

