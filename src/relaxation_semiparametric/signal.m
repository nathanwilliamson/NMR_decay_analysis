function I = signal(t, model, param)

I = zeros(size(t));

number_of_components = numel(model);

theta = param(end-number_of_components:end-1);
I0 = param(end);

ind = 1;
for current_component = 1:number_of_components
    switch model{current_component}{1}
        case 'exponential'
            T = param(ind);
            ind = ind + 1;
            
            I = I + theta(current_component) * exp(-t./T);
        case 'stretchedexponential'
            T = param(ind);
            beta = param(ind+1);
            ind = ind + 2;
                       
            I = I + theta(current_component) * exp(-(t./T).^beta);
        case 'lognormal'
            mu = param(ind);
            sigma = param(ind+1);
            ind = ind + 2;
            
            number_of_points = 1000;
            
            Tmax = exp(mu - sigma^2);
            Tspace = logspace(log10(Tmax) - 3 * sigma, log10(Tmax) + 3 * sigma, number_of_points);
            [karray, Darray] = ndgrid(t, Tspace);
            integrandarray = 1 ./ (Darray .* sigma * sqrt(2 * pi)) .* exp(- 1 / 2 * ((log(Darray) - mu) ./ sigma).^2) .* exp(-karray .* Darray);
            meanintegrandarray = (integrandarray(:, 1:(number_of_points - 1)) + integrandarray(:, 2:number_of_points)) / 2;
            dT = Tspace(2:number_of_points) - Tspace(1:(number_of_points - 1));
            integral = meanintegrandarray * dT';
            
            I = I + theta(current_component) * integral;
        case 'inversegamma'
            alpha = param(ind);
            beta = param(ind+1);
            ind = ind + 2;
            
            I = I + theta(current_component) * (beta ./ (beta + t)).^alpha;
    end
end

I = I * I0;

baseline = param(ind);
I = I + baseline;

end