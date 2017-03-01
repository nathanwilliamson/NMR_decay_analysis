function val = randmeanD(b,I)

Dslope                  = (log(I(1))-log(I(2)))/(b(2)-b(1));
val                     = 2 * Dslope * rand();