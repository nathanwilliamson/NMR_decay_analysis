function val = sumofsquares(b,I,type,param)

Imodel          = signal(b, type, param);
b
I
Imodel


val             = sum( (I - Imodel).^2 );

end

