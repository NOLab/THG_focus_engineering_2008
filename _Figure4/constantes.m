function[c,e0,mu0,lambda_1200,lambda_400,omega_1200,omega_400,w0,NA,n1_400,n1_1200,n2_400,n2_1200,Chi3_1,Chi3_2,E0,f,f0]= constantes()


load ctmat.mat
c=299792458;
e0=8.8542*10^(-12); %epsilon 0
mu0=4*pi*10^(-7); %mu 0 

lambda_1200=ctmat(2);
lambda_400=lambda_1200/3;
omega_1200 = 2*pi*c/lambda_1200;
omega_400 = 3*omega_1200;
n1_400=ctmat(3); %bk7
n1_1200=ctmat(4);
n2_400=n1_400;
n2_1200=n1_1200;
NA=ctmat(1);
Chi3_1=0;
Chi3_2=1;
f=0.002; %focale de l obj
E0=1 ;% Amplitude en entree
f0=ctmat(5); %facteur de remplissage de l obsj
w0=(f0*f*NA/n1_1200); % waist en entree d obj
