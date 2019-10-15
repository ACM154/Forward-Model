function output = smooth_ddelta1(xval,xinterval,a)
% 4-point cosine discrete delta function
% Reference: Yang et al (2009) Journal of Computational Physics, 228(20), 7821-7836.
    h=xinterval(2)-xinterval(1);
    output=zeros(length(xval),1);
    region=abs((xval-a)/h)<2;
    r=(xval(region)-a)/h;
    output(region)=(1/4)*(1+cos(pi*r/2))/h;
end