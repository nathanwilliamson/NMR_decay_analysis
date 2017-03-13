function val = sumofsquares_lognormal(b,I,mu,sigma,I0,Ib)

val_signal = signal_lognormal(b,mu,sigma,I0,Ib);

val = sum( (I - val_signal).^2 );

end