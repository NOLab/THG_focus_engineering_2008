function [res]= integretoutBessel(longueur,largeur,b,bz)
'integration des matrices Aij d un faisceau de bessel'
[c,e0,mu0,lambda_1200,lambda_400,omega_1200,omega_400,w0,NA,n1_400,n1_1200,n2_400,n2_1200,Chi3_1,Chi3_2,E0,f,f0]= constantes;

alphamin = asin((NA*0.95)/n1_1200);
alphamax = asin(NA*1.05/n1_1200);

global rho ;
global zed;
tic;
larr=(1.5*largeur)+1;
A00B=sparse(larr,2*longueur+1);
A01B=sparse(larr,2*longueur+1);
A02B=sparse(larr,2*longueur+1);


for i=1:larr
    for j=0:longueur
        rho= b*(i-1);
        zed= bz*j;
        A00B(i,longueur+1+j)=quadl(@I00B,alphamin,alphamax);
        A00B(i,longueur+1-j)=conj(A00B(i,longueur+1+j));
        A01B(i,longueur+1+j)=quadl(@I01B,alphamin,alphamax);
        A01B(i,longueur+1-j)=conj(A01B(i,longueur+1+j));
        A02B(i,longueur+1+j)=quadl(@I02B,alphamin,alphamax);
        A02B(i,longueur+1-j)=conj(A02B(i,longueur+1+j));

    end
end

save A00B.mat A00B;
save A01B.mat A01B;
save A02B.mat A02B;


res='meuh';

