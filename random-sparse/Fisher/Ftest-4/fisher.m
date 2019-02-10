function rhs=fisher(t,ut,dummy,k,x)
u=ifft(ut);
alpha=1;
rhs=-alpha*(k.^2).*ut+fft(u)-fft(u.^2);