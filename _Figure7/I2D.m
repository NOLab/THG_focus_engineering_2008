function [res] = I2D(a)
% Appel des constantes
[c,e0,mu0,lambda_1200,lambda_400,omega_1200,omega_400,w0,NA,n1_400,n1_1200,n2_400,n2_1200,Chi3_1,Chi3_2,E0,f,f0]= constantes;
k=2*pi*n1_1200/lambda_1200;
%f0=1;
global rho ;
global zed;
global phi;
p=rho;
z=zed;
fi=phi;

res=exp(-n1_1200.*n1_1200.*sin(a).*sin(a)./(f0.*f0.*NA.*NA)).*sqrt(cos(a)).*besselj(2,(sin(a).*k.*p)).*exp(i*k*z*cos(a)).*sin(a).*sin(a).*sin(a).*(cos(a)+1);