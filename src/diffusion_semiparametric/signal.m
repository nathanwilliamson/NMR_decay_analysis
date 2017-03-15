function I = signal(b, model, param)

I = zeros(size(b));

number_of_components = numel(model);

theta = param(end-number_of_components:end-1);
I0 = param(end);

ind = 1;
for current_component = 1:number_of_components
    switch model{current_component}{1}
        case 'exponential'
            D = param(ind);
            ind = ind + 1;
            
            I = I + theta(current_component) * exp(-b.*D);
        case 'stretchedexponential'
            D = param(ind);
            beta = param(ind+1);
            ind = ind + 2;
                       
            I = I + theta(current_component) * exp(-(b.*D).^beta);
        case 'lognormal'
            mu = param(ind);
            sigma = param(ind+1);
            ind = ind + 2;
            
            number_of_points = 1000;
            
            Dmax = exp(mu - sigma^2);
            Dspace = logspace(log10(Dmax) - 3 * sigma, log10(Dmax) + 3 * sigma, number_of_points);
            [karray, Darray] = ndgrid(b, Dspace);
            integrandarray = 1 ./ (Darray .* sigma * sqrt(2 * pi)) .* exp(- 1 / 2 * ((log(Darray) - mu) ./ sigma).^2) .* exp(-karray .* Darray);
            meanintegrandarray = (integrandarray(:, 1:(number_of_points - 1)) + integrandarray(:, 2:number_of_points)) / 2;
            dD = Dspace(2:number_of_points) - Dspace(1:(number_of_points - 1));
            integral = meanintegrandarray * dD';
            
            I = I + theta(current_component) * integral;
        case 'gamma'
            alpha = param(ind);
            beta = param(ind+1);
            ind = ind + 2;
            
            I = I + theta(current_component) * (beta ./ (beta + b)).^alpha;
    end
end

I = I * I0;

baseline = param(ind);
I = I + baseline;

end