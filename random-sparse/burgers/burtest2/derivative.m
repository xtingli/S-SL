function der=derivative(x,k)

if k<=1&&k>=1
    f_0=sech(x);
    f_1=-2*(exp(x)-exp(-x))./(exp(x)+exp(-x)).^2;
    f_2=-2./(exp(x)+exp(-x))+4*(exp(x)+exp(-x)).^2./(exp(x)+exp(-x)).^3;
    f_3=10*(exp(x)-exp(-x))./(exp(x)+exp(-x))-12*(exp(x)-exp(-x)).^3./(exp(x)+exp(-x)).^4;
    f_4=10./(exp(x)+exp(-x))-56*(exp(x)-exp(-x)).^2./(exp(x)+exp(-x)).^3+...
        48*(exp(x)-exp(-x)).^4./(exp(x)+exp(-x)).^5;
    der=[f_0.',f_1.',f_2.'];
elseif k<=2&&k>=2
        f_0=3*sech(x);
        f_1=3*(-2*(exp(x)-exp(-x))./(exp(x)+exp(-x)).^2);
        f_2=3*(-2./(exp(x)+exp(-x))+4*(exp(x)-exp(-x)).^2./(exp(x)+exp(-x)).^3);
        f_3=3*(2*(exp(x)-exp(-x))./(exp(x)+exp(-x))-12*(exp(x)-exp(-x)).^3./(exp(x)+exp(-x)).^4);
        f_4=3*(10./(exp(x)+exp(-x))-56*(exp(x)-exp(-x)).^2./(exp(x)+exp(-x)).^3+...
            48*(exp(x)-exp(-x)).^4./(exp(x)+exp(-x)).^5);
        der=[f_0.',f_1.',f_2.'];
elseif k<=3&&k>=3
            f_0=exp(-(x+2).^2);
            f_1=-2*(x+2).*exp(-(x+2).^2);
            f_2=-2*exp(-(x+2).^2)+4*(x+2).^2.*exp(-(x+2).^2);
            f_3=4*(x+2).*exp(-(x+2).^2)+8*(x+2).*exp(-(x+2).^2)-8*(x+2).^3.*exp(-(x+2).^2);
            f_4=4*exp(-(x+2).^2)-8*(x+2).^2.*exp(-(x+2).^2)+8*exp(-(x+2).^2)-16*(x+2).^2.*...
                exp(-(x+2).^2)-24*(x+2).^2.*exp(-(x+2).^2)+16*(x+2).^4.*exp(-(x+2).^2);
            der=[f_0.',f_1.',f_2.'];
elseif k<=4&&k>=4
                f_0=0.2*(exp(-(x+2).^2));
                f_1=0.2*(-2*(x+2).*exp(-(x+2).^2));
                f_2=0.2*(-2*exp(-(x+2).^2)+4*(x+2).^2.*exp(-(x+2).^2));
                f_3=0.2*(4*(x+2).*exp(-(x+2).^2)+8*(x+2).*exp(-(x+2).^2)-8*(x+2).^3.*exp(-(x+2).^2));
                f_4=0.2*(4*exp(-(x+2).^2)-8*(x+2).^2.*exp(-(x+2).^2)+8*exp(-(x+2).^2)-16*(x+2).^2.*...
                exp(-(x+2).^2)-24*(x+2).^2.*exp(-(x+2).^2)+16*(x+2).^4.*exp(-(x+2).^2));
               der=[f_0.',f_1.',f_2.'];
elseif k<=5&&k>=5
        f_0=2*sech(x);
        f_1=2*(-2*(exp(x)-exp(-x))./(exp(x)+exp(-x)).^2);
        f_2=2*(-2./(exp(x)+exp(-x))+4*(exp(x)-exp(-x)).^2./(exp(x)+exp(-x)).^3);
        f_3=2*(2*(exp(x)-exp(-x))./(exp(x)+exp(-x))-12*(exp(x)-exp(-x)).^3./(exp(x)+exp(-x)).^4);
        f_4=2*(10./(exp(x)+exp(-x))-56*(exp(x)-exp(-x)).^2./(exp(x)+exp(-x)).^3+...
            48*(exp(x)-exp(-x)).^4./(exp(x)+exp(-x)).^5);
        der=[f_0.',f_1.',f_2.'];
else  
                                f_0=sin(pi*x);
                    f_1=pi*cos(pi*x);
                    f_2=-pi.^2*sin(pi*x);
                    f_3=-pi.^3*cos(pi*x);
                    f_4=pi.^4*sin(pi*x);
                    der=[f_0.',f_1.',f_2.'];
end
end