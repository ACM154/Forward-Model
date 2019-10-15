function output = smooth_ddelta2(xval,xinterval,a)
% 4-point cosine smoothed discrete delta function
% Reference: Yang et al (2009) Journal of Computational Physics, 228(20), 7821-7836.
    h=xinterval(2)-xinterval(1);
    output=zeros(length(xval),1);
    region1=abs((xval-a)/h)<1.5;
    region2=(abs((xval-a)/h)>=1.5) & (abs((xval-a)/h)<=2.5);
    r1=(xval(region1)-a)/h;
    r2=(xval(region2)-a)/h;
    output(region1)=(1/(4*pi))/h*(pi+2*sin(pi/4*(2*r1+1))-2*sin(pi/4*(2*r1-1)));
    output(region2)=(-1/(8*pi))/h*(-5*pi+2*pi*abs(r2)+4*sin(pi/4*(2*abs(r2)-1)));
end