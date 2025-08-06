function A= periode(freq,largeur,longueur)
larg=2*largeur+1;
long=2*longueur+1;
 A=ones(larg,larg,long);
 if (freq==0)
b='meuh'
 elseif (freq==0.5)
     for k=1:longueur
         A(:,:,k)=zeros(larg,larg);
     end
 else
 for c=1:long
A(:,:,c)=(1+cos(pi*(c-longueur-1)/freq))*ones(larg,larg);
end
end
