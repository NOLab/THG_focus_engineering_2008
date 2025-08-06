function [res]= integretoutcyl(longueur,largeur,b,bz)
'integration des matrices Aij'
[c,e0,mu0,lambda_1200,lambda_400,omega_1200,omega_400,w0,NA,n1_400,n1_1200,n2_400,n2_1200,Chi3_1,Chi3_2,E0,f,f0]= constantes;
alpha = asin(NA/n1_1200);
%k=2*pi*n1/lambda;
%f0=1;
%ctte=0.5.*E0.*i.*k.*f.*(f./w0).*exp(-i.*k.*f);
global rho ;
global zed;
tic;
A00=sparse(floor(1.5*largeur),2*longueur+1);
A01=sparse(floor(1.5*largeur),2*longueur+1);
A02=sparse(floor(1.5*largeur),2*longueur+1);
A10=sparse(floor(1.5*largeur),2*longueur+1);
A11=sparse(floor(1.5*largeur),2*longueur+1);
A12=sparse(floor(1.5*largeur),2*longueur+1);
A13=sparse(floor(1.5*largeur),2*longueur+1);
A14=sparse(floor(1.5*largeur),2*longueur+1);

for i=1:floor(1.5*largeur)
    for j=0:longueur
        rho= b*(i-1);
        zed= bz*j;
A00(i,longueur+1+j)=quadl(@I00,0,alpha);
A00(i,longueur+1-j)=conj(A00(i,longueur+1+j));
A01(i,longueur+1+j)=quadl(@I01,0,alpha);
A01(i,longueur+1-j)=conj(A01(i,longueur+1+j));
A02(i,longueur+1+j)=quadl(@I02,0,alpha);
A02(i,longueur+1-j)=conj(A02(i,longueur+1+j));
A10(i,longueur+1+j)=quadl(@I10,0,alpha);
A10(i,longueur+1-j)=conj(A10(i,longueur+1+j));
A11(i,longueur+1+j)=quadl(@I11,0,alpha);
A11(i,longueur+1-j)=conj(A11(i,longueur+1+j));
A12(i,longueur+1+j)=quadl(@I12,0,alpha);
A12(i,longueur+1-j)=conj(A12(i,longueur+1+j));
A13(i,longueur+1+j)=quadl(@I13,0,alpha);
A13(i,longueur+1-j)=conj(A13(i,longueur+1+j));
A14(i,longueur+1+j)=quadl(@I14,0,alpha);
A14(i,longueur+1-j)=conj(A14(i,longueur+1+j));
end
end

save A00.mat A00;
save A01.mat A01;
save A02.mat A02;
save A10.mat A10;
save A11.mat A11;
save A12.mat A12;
save A13.mat A13;
save A14.mat A14;

res='meuh';

