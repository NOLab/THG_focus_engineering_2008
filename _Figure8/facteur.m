function fa = facteur(x,y,z)
[c,e0,mu0,lambda_1200,lambda_400,omega_1200,omega_400,w0,NA,n1_400,n1_1200,n2_400,n2_1200,Chi3_1,Chi3_2,E0,f,f0]= constantes;
fa=(omega_400*omega_400/(c*c))*exp(i*2*pi*(n1_400/lambda_400)*(z+(x*x+y*y)/(2*z)))/z;