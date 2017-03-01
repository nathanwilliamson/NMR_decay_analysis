function I = signal(b,model,param)

I                               = zeros(size(b));

nComponents                     = numel(model);

theta                           = param(end-nComponents:end-1);
I0                              = param(end);

ind                             = 1;
for currentComponent = 1:nComponents
    switch model{currentComponent}{1}
        case 'exponential'
            D                       = param(ind);
            ind                     = ind + 1;
            
            I                       = I + theta(currentComponent) * exp(-b.*D);
        case 'stretchedexponential'
            D                       = param(ind);
            beta                    = param(ind+1);
            ind                     = ind + 2;
                       
            I                       = I + theta(currentComponent) * exp(-(b.*D).^beta);
        case 'lognormal'
            mu                      = param(ind);
            sigma                   = param(ind+1);
            ind                     = ind + 2;
            
            nPoints                 = 1000;
            
            Dmax                    = exp(mu-sigma^2);
            Dspace                  = logspace(log10(Dmax)-3*sigma,log10(Dmax)+3*sigma,nPoints);
            [karray, Darray]        = ndgrid(b,Dspace);
            integrandarray          = 1./(Darray.*sigma*sqrt(2*pi)).*exp(-1/2*((log(Darray)-mu)./sigma).^2).*exp(-karray.*Darray);
            meanintegrandarray      = (integrandarray(:,1:(nPoints-1)) + integrandarray(:,2:nPoints))/2;
            dD                      = Dspace(2:nPoints) - Dspace(1:(nPoints-1));
            integral                = meanintegrandarray*dD';
            
            I                       = I + theta(currentComponent) * integral;
        case 'gamma'
            alpha                   = param(ind);
            beta                    = param(ind+1);
            ind                     = ind + 2;
            
            I                       = I + theta(currentComponent) * (beta./(beta+b)).^alpha;
    end
end

I                       = I * I0;

b                       = param(ind);
I                       = I + b;

end