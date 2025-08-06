function [C] = interfaceangle2(alpha,longueur,largeur)
larg=2*largeur+1;
long=2*longueur+1;
C=ones(larg,larg,long);
for z=1:long
    C(:,:,z)=interfaceangle(alpha,largeur);
end

