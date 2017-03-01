function val = randmeanD(b,I)

ind = round(numel(b)/4);
Dslope = (log(I(ind))-log(I(1)))/(b(ind)-b(1));
val = 3 * abs(Dslope) * rand();