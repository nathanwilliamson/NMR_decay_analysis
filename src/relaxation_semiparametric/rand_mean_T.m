function val = rand_mean_T(t)

minT = t(1);
maxT = t(end);
val = minT + (maxT - minT) * rand();

end