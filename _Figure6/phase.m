function A= phase(x,y,z,larg,long,b,bz)
lar=2*larg+1;
lon=2*long+1;
[c,e0,mu0,lambda_1200,lambda_400,omega_1200,omega_400,w0,NA,n1_400,n1_1200,n2_400,n2_1200,Chi3_1,Chi3_2,E0,f,f0]= constantes;
A=ones(lar,lar,lon);
%for g=1:q
%    for j=1:q
%       for k=1:q
%        A(g,j,k)=exp(i*(6*pi*n1/lambda)*b*(k-(m+1)+((g-(m+1))*x+(j-(m+1))*y)/z));
%       end
%   end
%end

for g=1:lar
   for j=1:lar
      for k=1:lon
       A(g,j,k)=exp(-i*(2*pi*n1_400/lambda_400)*(bz*z*(k-long-1) +(g-larg-1)*x*b +(j-larg-1)*y*b)/sqrt(abs(z)^2+abs(y)^2+abs(x)^2));
      end
    end
end