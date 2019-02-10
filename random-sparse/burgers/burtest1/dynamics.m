function rhs=dynamics(t,ut,dummy,k,x)

% @author: Xiuting Li
%
%Inputs:
%      ut    : the Fourier transformation of the system state
%      x     : space
%      k     : the frequency vector
%
%Output:
%      rhs   : the system dynamics
%
u=ifft(ut);
alpha=0.1;
%mu=-(x'.^2+x'.^4);
rhs=-alpha*(k.^2).*ut-1i.*k.*fft(0.5.*u.^2);