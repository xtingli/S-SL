function rhs=fisher(t,ut,dummy,k,x)
u=ifft(ut);
alpha=10;
c_1=1;
c_2=10;
rhs=-alpha*(k.^2).*ut+c_1*fft(u)-c_2*fft(u.^2);