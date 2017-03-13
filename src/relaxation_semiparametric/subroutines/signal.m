function I = signal(t,type,baseline,param)

I                               = zeros(size(t));

nComponents                     = numel(type);

theta                           = param(end-nComponents:end-1);
I0                              = param(end);

ind                             = 1;
for currentComponent = 1:nComponents
    switch type{currentComponent}
        case 'exponential'
            tau                     = param(ind);
            ind                     = ind + 1;
            
            I                       = I + theta(currentComponent) * exp(-t/tau);
        case 'stretchedexponential'
            tau                     = param(ind);
            beta                    = param(ind+1);
            ind                     = ind + 2;
                       
            I                       = I + theta(currentComponent) * exp(-(t/tau).^beta);
        case 'lognormal'
            mu                      = param(ind);
            sigma                   = param(ind+1);
            ind                     = ind + 2;
            
            nPoints                 = 1000;
            
            taumax                  = exp(mu-sigma^2);
            tauspace                = logspace(log10(taumax)-3*sigma,log10(taumax)+3*sigma,nPoints);
            [tarray, tauarray]      = ndgrid(t,tauspace);
            integrandarray          = 1./(tauarray.*sigma*sqrt(2*pi)).*exp(-1/2*((log(tauarray)-mu)./sigma).^2).*exp(-tarray./tauarray);
            meanintegrandarray      = (integrandarray(:,1:(nPoints-1)) + integrandarray(:,2:nPoints))/2;
            dtau                    = tauspace(2:nPoints) - tauspace(1:(nPoints-1));
            integral                = meanintegrandarray*dtau';
            
            I                       = I + theta(currentComponent) * integral;
        case 'inversegamma'
            alpha                   = param(ind);
            beta                    = param(ind+1);
            ind                     = ind + 2;
            
            I                       = I + theta(currentComponent) * (beta./(beta+t)).^alpha;
    end
end

I                           = I * I0;

if baseline
    b                       = param(ind);
    I                       = I + b;
end            

end