function [res]= integretoutmode2(longueur,largeur,b,bz)
'integration des matrices Aij'
[c,e0,mu0,lambda_1200,lambda_400,omega_1200,omega_400,w0,NA,n1_400,n1_1200,n2_400,n2_1200,Chi3_1,Chi3_2,E0,f,f0]= constantes;
alpha=asin(NA/n1_1200);
alphain=0.616 %NA=1.4
%alphain=0.521; %NA=1.14 dephasage de pi
%alphain=0.478; %dephasage de 0.8 pi 
inter=cos(alpha);
larr=(1.5*largeur)+1;

global rho ;
global zed;
tic;
A2A=sparse(larr,2*longueur+1);
A2B=sparse(larr,2*longueur+1);
A2C=sparse(larr,2*longueur+1);
A2D=sparse(larr,2*longueur+1);
A2E=sparse(larr,2*longueur+1);
A2F=sparse(larr,2*longueur+1);
A2G=sparse(larr,2*longueur+1);


for i=1:larr

    for j=0:longueur
        rho= b*(i-1);
        zed= bz*j;
A2A(i,longueur+1+j)=quadl(@I2A,0,alpha);
A2A(i,longueur+1-j)=conj(A2A(i,longueur+1+j));
A2B(i,longueur+1+j)=quadl(@I2B,0,alpha);
A2B(i,longueur+1-j)=conj(A2B(i,longueur+1+j));
A2C(i,longueur+1+j)=quadl(@I2C,0,alpha);
A2C(i,longueur+1-j)=conj(A2C(i,longueur+1+j));
A2D(i,longueur+1+j)=quadl(@I2D,0,alpha);
A2D(i,longueur+1-j)=conj(A2D(i,longueur+1+j));
A2E(i,longueur+1+j)=quadl(@I2E,0,alpha);
A2E(i,longueur+1-j)=conj(A2E(i,longueur+1+j));
A2F(i,longueur+1+j)=quadl(@I2F,0,alpha);
A2F(i,longueur+1-j)=conj(A2F(i,longueur+1+j));
A2G(i,longueur+1+j)=quadl(@I2G,0,alpha);
A2G(i,longueur+1-j)=conj(A2G(i,longueur+1+j));

end
end

save A2A.mat A2A;
save A2B.mat A2B;
save A2C.mat A2C;
save A2D.mat A2D;
save A2E.mat A2E;
save A2F.mat A2F;
save A2G.mat A2G;


res='meuh';

