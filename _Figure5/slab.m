function A= slab(direction,epaisseur,largeur,longueur)
larg=2*largeur+1;
long=2*longueur+1;
 A=zeros(larg,larg,long);
Chi3_2=1;

if(direction=='z')
 for c=longueur+1-epaisseur:longueur+1+epaisseur
A(:,:,c)=ones(larg,larg);
end
end


if(direction=='x')
    A=zeros(larg,larg,long);
    for k=longueur+1-epaisseur:longueur+1+epaisseur
       for l=1:larg
          for c=1:long
              A(k,l,c)=Chi3_2;
          end
       end
   end
end

if(direction=='y')
   A=zeros(larg,larg,long);
    for k=1:larg
       for l=longueur+1-epaisseur:longueur+1+epaisseur
          for c=1:long
              A(k,l,c)=Chi3_2;
          end
       end
   end
end
