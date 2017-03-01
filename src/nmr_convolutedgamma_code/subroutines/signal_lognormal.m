function val = signal_lognormal(b,mu,sigma,I0,Ib)

nPoints = 1000;

Dmax = exp(mu-sigma^2);
Dspace = logspace(log10(Dmax)-3*sigma,log10(Dmax)+3*sigma,nPoints);
[karray, Darray] = ndgrid(b,Dspace);
integrandarray = 1./(Darray.*sigma*sqrt(2*pi)).*exp(-1/2*((log(Darray)-mu)./sigma).^2).*exp(-karray.*Darray);
meanintegrandarray = (integrandarray(:,1:(nPoints-1)) + integrandarray(:,2:nPoints))/2;
dD = Dspace(2:nPoints) - Dspace(1:(nPoints-1));
integral = meanintegrandarray*dD';

val = I0 * integral + Ib;

val = reshape(val,size(b));

end

