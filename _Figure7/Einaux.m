function [res] =Einaux(theta,phi,option)
[c,e0,mu0,lambda_1200,lambda_400,omega_1200,omega_400,w0,NA,n1_400,n1_1200,n2_400,n2_1200,Chi3_1,Chi3_2,E0,f,f0]= constantes;
thetamax=asin(NA/n1_1200);
k=2*pi*n1_1200/lambda_1200;
res=exp(-n1_1200.*n1_1200.*sin(theta).*sin(theta)./(f0.*f0.*NA.*NA)).*(2*(cos(phi).*sin(theta)>=0)-1);