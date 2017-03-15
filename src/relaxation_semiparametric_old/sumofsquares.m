function val = sumofsquares(t,I,type,baseline,param)

Imodel          = signal(t,type,baseline,param);

val             = sum( (I - Imodel).^2 );

end

