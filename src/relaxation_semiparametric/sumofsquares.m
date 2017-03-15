function val = sumofsquares(t, I, type, param)

Imodel = signal(t, type, param);

val = sum( (I - Imodel).^2 );

end

