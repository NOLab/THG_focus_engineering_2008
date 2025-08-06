function res=phasee(x)
[c,e0,mu0,lambda_1200,lambda_400,omega_1200,omega_400,w0,NA,n1_400,n1_1200,n2_400,n2_1200,Chi3_1,Chi3_2,E0,f,f0]= constantes;
k=2*pi*n1_1200/lambda_1200;
res=exp(-i.*k*x);
