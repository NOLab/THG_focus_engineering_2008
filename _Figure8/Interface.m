function [A]=Interface(direction,position,largeur,longueur)
[c,e0,mu0,lambda_1200,lambda_400,omega_1200,omega_400,w0,NA,n1_400,n1_1200,n2_400,n2_1200,Chi3_1,Chi3_2,E0,f,f0]= constantes;
larg=(2*largeur)+1;
long=(2*longueur)+1;
if(direction=='x')
    A=zeros(larg,larg,long);
    for k=1:position
       for l=1:larg
          for c=1:long
              A(k,l,c)=1;
          end
       end
   end
end

if(direction=='y')
   A=zeros(larg,larg,long);
    for k=1:larg
       for l=1:position
          for c=1:long
              A(k,l,c)=1;
          end
       end
   end
end

if(direction=='z')
    A=zeros(larg,larg,long);
    for k=1:position
        A(:,:,k)=ones(larg,larg);
    end
end
