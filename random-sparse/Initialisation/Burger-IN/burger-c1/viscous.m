function rhs=viscous(t,ut,dummy,k,x)
u=ifft(ut);
alpha=0.1;
%mu=-(x'.^2+x'.^4);
rhs=-alpha*(k.^2).*ut-1i.*k.*fft(0.5.*u.^2);