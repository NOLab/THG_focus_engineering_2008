function [res]= integretout3phase(longueur,largeur,b,bz,a1,a2)
'integration des matrices Aij d un faisceau 3phase'
[c,e0,mu0,lambda_1200,lambda_400,omega_1200,omega_400,w0,NA,n1_400,n1_1200,n2_400,n2_1200,Chi3_1,Chi3_2,E0,f,f0]= constantes;
alpha=asin(NA/n1_1200);
alphamax=alpha;
alphain=a1;
alphaout=a2;


larr=floor((1.5*largeur))+1;

global rho ;
global zed;
A00in=sparse(larr,2*longueur+1);
A01in=sparse(larr,2*longueur+1);
A02in=sparse(larr,2*longueur+1);

A00mid=sparse(larr,2*longueur+1);
A01mid=sparse(larr,2*longueur+1);
A02mid=sparse(larr,2*longueur+1);

A00out=sparse(larr,2*longueur+1);
A01out=sparse(larr,2*longueur+1);
A02out=sparse(larr,2*longueur+1)
for i=1:larr

    for j=0:longueur
        rho= b*(i-1);
        zed= bz*j;



        A00in(i,longueur+1+j)=quadgk(@I00,0,alphain);
        A00in(i,longueur+1-j)=conj(A00in(i,longueur+1+j));
        A01in(i,longueur+1+j)=quadgk(@I01,0,alphain);
        A01in(i,longueur+1-j)=conj(A01in(i,longueur+1+j));
        A02in(i,longueur+1+j)=quadgk(@I02,0,alphain);
        A02in(i,longueur+1-j)=conj(A02in(i,longueur+1+j));

        A00mid(i,longueur+1+j)=quadgk(@I00,alphain,alphaout);
        A00mid(i,longueur+1-j)=conj(A00mid(i,longueur+1+j));
        A01mid(i,longueur+1+j)=quadgk(@I01,alphain,alphaout);
        A01mid(i,longueur+1-j)=conj(A01mid(i,longueur+1+j));
        A02mid(i,longueur+1+j)=quadgk(@I02,alphain,alphaout);
        A02mid(i,longueur+1-j)=conj(A02mid(i,longueur+1+j));


        A00out(i,longueur+1+j)=quadgk(@I00,alphaout,alphamax);
        A00out(i,longueur+1-j)=conj(A00out(i,longueur+1+j));
        A01out(i,longueur+1+j)=quadgk(@I01,alphaout,alphamax);
        A01out(i,longueur+1-j)=conj(A01out(i,longueur+1+j));
        A02out(i,longueur+1+j)=quadgk(@I02,alphaout,alphamax);
        A02out(i,longueur+1-j)=conj(A02out(i,longueur+1+j));


    end
end
save A00in.mat A00in;
save A01in.mat A01in;
save A02in.mat A02in;

save A00mid.mat A00mid;
save A01mid.mat A01mid;
save A02mid.mat A02mid;

save A00out.mat A00out;
save A01out.mat A01out;
save A02out.mat A02out;

clear A00in
clear A01in
clear A02in

clear A00mid
clear A01mid
clear A02mid

clear A00out
clear A01out
clear A02out
res='meuh';