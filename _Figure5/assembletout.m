function [res]= assembletout(longueur,largeur,b,bz,option);
'debut assemblages du maillage d excitatyion'
%debut du chrono et declaration des constantes

tic;
[c,e0,mu0,lambda_1200,lambda_400,omega_1200,omega_400,w0,NA,n1_400,n1_1200,n2_400,n2_1200,Chi3_1,Chi3_2,E0,f,f0]= constantes;

alpha = asin(NA/n1_1200);
k=2*pi*n1_1200/lambda_1200;

ctteH=0.5.*E0.*i.*k.*f.*(f./w0).*exp(-i.*k.*f);
ct=0.5.*i.*k.*f.*E0.*exp(-i.*k.*f);

Itot=1e9;


res='meuh';

%appel de la fonction integretout qui calcul les matrices Aij pour le
%maillage n*b donne
if (strcmp(option,'custom')==0)
    [res]= integretoutv2(longueur,largeur,b,bz);
end

if (strcmp(option,'Bessel'))
    [res]= integretoutbessel(longueur,largeur,b,bz);    
    load A00B.mat ;
    load A01B.mat;
    load A02B.mat;
end

if (strcmp(option,'Mode20'))
    [res]= integretoutmode2(longueur,largeur,b,bz);
    load A2A.mat ;
    load A2B.mat;
    load A2C.mat;
    load  A2D.mat;
    load  A2E.mat ;
    load  A2F.mat ;
    load  A2G.mat ;

end

if (strcmp(option,'3phase'))
    load alpha.mat
    a1=alpha(1);
    a2=alpha(2);
    [res]= integretout3phase(longueur,largeur,b,bz,a1,a2);
end


%recuperation des matrices calculees par integretout

load A00.mat ;
load A01.mat;
load A02.mat;
load  A10.mat;
load  A11.mat ;
load  A12.mat ;
load  A13.mat ;
load  A14.mat ;



if strcmp(option,'Hermit')
    'calcul de foyer d excitation avec mode Hermit Gaussien'

    %initialisation des matrices resultats

    Ex=zeros(2*largeur+1,2*largeur+1,2*longueur+1);
    Ey=zeros(2*largeur+1,2*largeur+1,2*longueur+1);
    Ez=zeros(2*largeur+1,2*largeur+1,2*longueur+1);
    for ix=1:2*largeur+1
        for ygrec=1:2*largeur+1
            for zed=1:2*longueur+1

                %calcul de \rho  et \phi en fonction de x et y

                %                 if (round(sqrt((ix-largeur-1)*(ix-largeur-1)+(ygrec-largeur-1)*(ygrec-largeur-1)))>largeur)
                %                     res='meuh';
                %                 else
                rau=1+round(sqrt((ix-largeur-1)*(ix-largeur-1)+(ygrec-largeur-1)*(ygrec-largeur-1)));
                if(ix==(largeur+1))
                    if(ygrec>=(largeur+1))
                        fi=pi/2;
                    else
                        fi=-pi/2;
                    end
                else
                    fi=atan2((ygrec-largeur-1),(ix-largeur-1));
                end

                %Remplissage de la la matrice

                Ex(ix,ygrec,zed)=ctteH*(i*A11(rau,zed)*cos(fi)+i*A14(rau,zed)*cos(3*fi)-(A11(rau,zed)+2*A12(rau,zed))*sin(fi)-A14(rau,zed)*sin(3*fi));
                Ey(ix,ygrec,zed)=ctteH*(-i*A12(rau,zed)*sin(fi)+i*A14(rau,zed)*sin(3*fi)+A12(rau,zed)*cos(fi)+A14(rau,zed)*cos(3*fi));
                Ez(ix,ygrec,zed)=ctteH*(-2*A10(rau,zed)+2*A13(rau,zed)*(cos(2*fi)+sin(2*fi)));

            end
        end
    end




end
if  strcmp(option,'Gauss')
    'calcul de foyer d excitation avec mode Gaussien'
    Ex=zeros(2*largeur+1,2*largeur+1,2*longueur+1);
    Ey=zeros(2*largeur+1,2*largeur+1,2*longueur+1);
    Ez=zeros(2*largeur+1,2*largeur+1,2*longueur+1);
    for ix=1:2*largeur+1
        for ygrec=1:2*largeur+1
            for zed=1:2*longueur+1

                %calcul de \rho  et \phi en fonction de x et y

                %                 if (round(sqrt((ix-largeur-1)*(ix-largeur-1)+(ygrec-largeur-1)*(ygrec-largeur-1)))>largeur)
                %                     res='meuh';
                %                 else
                rau=1+round(sqrt((ix-largeur-1)*(ix-largeur-1)+(ygrec-largeur-1)*(ygrec-largeur-1)));
                if(ix==(largeur+1))
                    if(ygrec>=(largeur+1))
                        fi=pi/2;
                    else
                        fi=-pi/2;
                    end
                else
                    fi=atan2((ygrec-largeur-1),(ix-largeur-1));
                end

                %Remplissage de la la matrice

                Ex(ix,ygrec,zed)=ct*(A00(rau,zed)+A02(rau,zed)*cos(2*fi));
                Ey(ix,ygrec,zed)=ct*(A02(rau,zed)*sin(2*fi));
                Ez(ix,ygrec,zed)=ct*(-2*i*A01(rau,zed)*cos(fi));

            end
        end
    end
end

if  strcmp(option,'Bessel')
    'calcul de foyer d excitation avec mode Bessel'
    Ex=zeros(2*largeur+1,2*largeur+1,2*longueur+1);
    Ey=zeros(2*largeur+1,2*largeur+1,2*longueur+1);
    Ez=zeros(2*largeur+1,2*largeur+1,2*longueur+1);
    for ix=1:2*largeur+1
        for ygrec=1:2*largeur+1
            for zed=1:2*longueur+1

                %calcul de \rho  et \phi en fonction de x et y

                %                 if (round(sqrt((ix-largeur-1)*(ix-largeur-1)+(ygrec-largeur-1)*(ygrec-largeur-1)))>largeur)
                %                     res='meuh';
                %                 else
                rau=1+round(sqrt((ix-largeur-1)*(ix-largeur-1)+(ygrec-largeur-1)*(ygrec-largeur-1)));
                if(ix==(largeur+1))
                    if(ygrec>=(largeur+1))
                        fi=pi/2;
                    else
                        fi=-pi/2;
                    end
                else
                    fi=atan2((ygrec-largeur-1),(ix-largeur-1));
                end

                %Remplissage de la la matrice

                Ex(ix,ygrec,zed)=ct*(A00B(rau,zed)+A02B(rau,zed)*cos(2*fi));
                Ey(ix,ygrec,zed)=ct*(A02B(rau,zed)*sin(2*fi));
                Ez(ix,ygrec,zed)=ct*(-2*i*A01B(rau,zed)*cos(fi));

            end
        end
    end
end


if  strcmp(option,'Mode01')
    Ex=zeros(2*largeur+1,2*largeur+1,2*longueur+1);
    Ey=zeros(2*largeur+1,2*largeur+1,2*longueur+1);
    Ez=zeros(2*largeur+1,2*largeur+1,2*longueur+1);

    ct=i;
    'calcul de foyer d excitation avec mode 01'

    for ix=1:2*largeur+1
        for ygrec=1:2*largeur+1
            for zed=1:2*longueur+1

                %calcul de \rho  et \phi en fonction de x et y

                %                 if (round(sqrt((ix-largeur-1)*(ix-largeur-1)+(ygrec-largeur-1)*(ygrec-largeur-1)))>largeur)
                %                     res='meuh';
                %                 else
                rau=1+round(sqrt((ix-largeur-1)*(ix-largeur-1)+(ygrec-largeur-1)*(ygrec-largeur-1)));
                if(ix==(largeur+1))
                    if(ygrec>=(largeur+1))
                        fi=pi/2;
                    else
                        fi=-pi/2;
                    end
                else
                    fi=atan2((ygrec-largeur-1),(ix-largeur-1));
                end

                %Remplissage de la la matrice

                Ex(ix,ygrec,zed)=ctteH*((A11(rau,zed)+2*A12(rau,zed))*i*sin(fi)+A14(rau,zed)*sin(3*fi));
                Ey(ix,ygrec,zed)=ctteH*(-i*A12(rau,zed)*cos(fi)-i*A14(rau,zed)*cos(3*fi));
                Ez(ix,ygrec,zed)=ctteH*(2*A13(rau,zed)*sin(2*fi));

            end
        end
    end

end

if strcmp(option,'Mode10')
    Ex=zeros(2*largeur+1,2*largeur+1,2*longueur+1);
    Ey=zeros(2*largeur+1,2*largeur+1,2*longueur+1);
    Ez=zeros(2*largeur+1,2*largeur+1,2*longueur+1);

    ct=i;
    'calcul de foyer d excitation avec mode 10'

    for ix=1:2*largeur+1
        for ygrec=1:2*largeur+1
            for zed=1:2*longueur+1

                %calcul de \rho  et \phi en fonction de x et y

                %                 if (round(sqrt((ix-largeur-1)*(ix-largeur-1)+(ygrec-largeur-1)*(ygrec-largeur-1)))>largeur)
                %                     res='meuh';
                %                 else
                rau=1+round(sqrt((ix-largeur-1)*(ix-largeur-1)+(ygrec-largeur-1)*(ygrec-largeur-1)));
                if(ix==(largeur+1))
                    if(ygrec>=(largeur+1))
                        fi=pi/2;
                    else
                        fi=-pi/2;
                    end
                else
                    fi=atan2((ygrec-largeur-1),(ix-largeur-1));
                end

                %Remplissage de la la matrice
                Ex(ix,ygrec,zed)=ctteH*((A11(rau,zed)*i*cos(fi)+A14(rau,zed)*i*cos(3*fi)));
                Ey(ix,ygrec,zed)=ctteH*(-i*A12(rau,zed)*sin(fi)+i*A14(rau,zed)*sin(3*fi));
                Ez(ix,ygrec,zed)=ctteH*(-2*A10(rau,zed)+2*i*A13(rau,zed)*cos(2*fi));

            end
        end
    end

end

if strcmp(option,'Mode10x')
    Ex=zeros(2*largeur+1,2*largeur+1,2*longueur+1);
    Ey=zeros(2*largeur+1,2*largeur+1,2*longueur+1);
    Ez=zeros(2*largeur+1,2*largeur+1,2*longueur+1);

    ct=i;
    'calcul de foyer d excitation avec mode 10'

    for ix=1:2*largeur+1
        for ygrec=1:2*largeur+1
            for zed=1:2*longueur+1

                %calcul de \rho  et \phi en fonction de x et y

                %                 if (round(sqrt((ix-largeur-1)*(ix-largeur-1)+(ygrec-largeur-1)*(ygrec-largeur-1)))>largeur)
                %                     res='meuh';
                %                 else
                rau=1+round(sqrt((ix-largeur-1)*(ix-largeur-1)+(ygrec-largeur-1)*(ygrec-largeur-1)));
                if(ix==(largeur+1))
                    if(ygrec>=(largeur+1))
                        fi=pi/2;
                    else
                        fi=-pi/2;
                    end
                else
                    fi=atan2((ygrec-largeur-1),(ix-largeur-1));
                end

                %Remplissage de la la matrice
                Ex(ix,ygrec,zed)=ctteH*((A11(rau,zed)*i*cos(fi)+A14(rau,zed)*i*cos(3*fi)));
                Ey(ix,ygrec,zed)=0;
                Ez(ix,ygrec,zed)=0;

            end
        end
    end

end

if strcmp(option,'RadDonut')
    Ex=zeros(2*largeur+1,2*largeur+1,2*longueur+1);
    Ey=zeros(2*largeur+1,2*largeur+1,2*longueur+1);
    Ez=zeros(2*largeur+1,2*largeur+1,2*longueur+1);

    ct=1;
    'calcul de foyer d excitation avec mode radonut'

    for ix=1:2*largeur+1
        for ygrec=1:2*largeur+1
            for zed=1:2*longueur+1

                %calcul de \rho  et \phi en fonction de x et y

                %                 if (round(sqrt((ix-largeur-1)*(ix-largeur-1)+(ygrec-largeur-1)*(ygrec-largeur-1)))>largeur)
                %                     res='meuh';
                %                 else
                rau=1+round(sqrt((ix-largeur-1)*(ix-largeur-1)+(ygrec-largeur-1)*(ygrec-largeur-1)));
                if(ix==(largeur+1))
                    if(ygrec>=(largeur+1))
                        fi=pi/2;
                    else
                        fi=-pi/2;
                    end
                else
                    fi=atan2((ygrec-largeur-1),(ix-largeur-1));
                    %fi=atan((ygrec-largeur-1)/(ix-largeur-1));

                end


                %Remplissage de la la matrice
                Ex(ix,ygrec,zed)=ctteH*((A11(rau,zed)-A12(rau,zed))*i*cos(fi));
                Ey(ix,ygrec,zed)=ctteH*(((A11(rau,zed)-A12(rau,zed))*i*sin(fi)));
                Ez(ix,ygrec,zed)=ctteH*(-4*A10(rau,zed));

            end
        end
    end

end

if strcmp(option,'AzDonut')
    Ex=zeros(2*largeur+1,2*largeur+1,2*longueur+1);
    Ey=zeros(2*largeur+1,2*largeur+1,2*longueur+1);
    Ez=zeros(2*largeur+1,2*largeur+1,2*longueur+1);

    ct=1;
    'calcul de foyer d excitation avec mode Azdonut'

    for ix=1:2*largeur+1
        for ygrec=1:2*largeur+1
            for zed=1:2*longueur+1

                %calcul de \rho  et \phi en fonction de x et y

                %                 if (round(sqrt((ix-largeur-1)*(ix-largeur-1)+(ygrec-largeur-1)*(ygrec-largeur-1)))>largeur)
                %                     res='meuh';
                %                 else
                rau=1+round(sqrt((ix-largeur-1)*(ix-largeur-1)+(ygrec-largeur-1)*(ygrec-largeur-1)));
                if(ix==(largeur+1))
                    if(ygrec>=(largeur+1))
                        fi=pi/2;
                    else
                        fi=-pi/2;
                    end
                else
                    fi=atan2((ygrec-largeur-1),(ix-largeur-1));
                    %fi=atan((ygrec-largeur-1)/(ix-largeur-1));
                end


                %Remplissage de la la matrice
                Ex(ix,ygrec,zed)=ctteH*(i*(A11(rau,zed)+3*A12(rau,zed))*sin(fi));
                Ey(ix,ygrec,zed)=ctteH*(((A11(rau,zed)+3*A12(rau,zed))*(-i)*cos(fi)));
                Ez(ix,ygrec,zed)=0;

            end
        end
    end

end
if strcmp(option,'2phase')

    Ex=zeros(2*largeur+1,2*largeur+1,2*longueur+1);
    Ey=zeros(2*largeur+1,2*largeur+1,2*longueur+1);
    Ez=zeros(2*largeur+1,2*largeur+1,2*longueur+1);

    load  A00in.mat;
    load  A01in.mat;
    load  A02in.mat;

    load  A00out.mat;
    load  A01out.mat;
    load  A02out.mat;


    ct=1;
    def=exp(i*pi);
    'calcul de foyer d excitation avec mode meux'

    for ix=1:2*largeur+1
        for ygrec=1:2*largeur+1
            for zed=1:2*longueur+1

                %calcul de \rho  et \phi en fonction de x et y
                rau=1+round(sqrt((ix-largeur-1)*(ix-largeur-1)+(ygrec-largeur-1)*(ygrec-largeur-1)));
                if(ix==(largeur+1))
                    if(ygrec>=(largeur+1))
                        fi=pi/2;
                    else
                        fi=-pi/2;
                    end
                else
                    fi=atan2((ygrec-largeur-1),(ix-largeur-1));
                end

                %Remplissage de la la matrice
                Ex(ix,ygrec,zed)=ct*A00out(rau,zed)+A02out(rau,zed)*cos(2*fi)+def*(A00in(rau,zed)+A02in(rau,zed)*cos(2*fi));
                Ey(ix,ygrec,zed)=ct*sin(2*fi)*(A02out(rau,zed)+def*A02in(rau,zed));
                Ez(ix,ygrec,zed)=-ct*2*i*cos(fi)*(A01out(rau,zed)+def*A01in(rau,zed));

            end
        end
    end

end
if strcmp(option,'3phase')

    Ex=zeros(2*largeur+1,2*largeur+1,2*longueur+1);
    Ey=zeros(2*largeur+1,2*largeur+1,2*longueur+1);
    Ez=zeros(2*largeur+1,2*largeur+1,2*longueur+1);

    load  A00in.mat;
    load  A01in.mat;
    load  A02in.mat;
    load  A00mid.mat;
    load  A01mid.mat;
    load  A02mid.mat;
    load  A00out.mat;
    load  A01out.mat;
    load  A02out.mat;


    load deph.mat
    defin=1;
    defmid=exp(i*pi*deph(1));
    defout=exp(i*pi*deph(2));
    'calcul de foyer d excitation avec masque de phase 3 elements'

    for ix=1:2*largeur+1
        for ygrec=1:2*largeur+1
            for zed=1:2*longueur+1

                %calcul de \rho  et \phi en fonction de x et y
                rau=1+round(sqrt((ix-largeur-1)*(ix-largeur-1)+(ygrec-largeur-1)*(ygrec-largeur-1)));
                if(ix==(largeur+1))
                    if(ygrec>=(largeur+1))
                        fi=pi/2;
                    else
                        fi=-pi/2;
                    end
                else
                    fi=atan2((ygrec-largeur-1),(ix-largeur-1));
                end

                %Remplissage de la la matrice
                Ex(ix,ygrec,zed)=(A00in(rau,zed)*defin+A00mid(rau,zed)*defmid+A00out(rau,zed)*defout)+(A02in(rau,zed)*defin+A02mid(rau,zed)*defmid+A02out(rau,zed)*defout)*cos(2*fi);
                Ey(ix,ygrec,zed)=sin(2*fi)*(A02in(rau,zed)*defin+A02mid(rau,zed)*defmid+A02out(rau,zed)*defout);
                Ez(ix,ygrec,zed)=-2*i*cos(fi)*(A01in(rau,zed)*defin+A01mid(rau,zed)*defmid+A01out(rau,zed)*defout);

            end
        end
    end

end

if strcmp(option,'custom')
    'on utilise les matrics  deja definies dans le workspace'
    load Ex.mat
    load Ey.mat
    load Ez.mat
end

if strcmp(option,'Mode20')
    'calcul de foyer d excitation avec mode20'

    %initialisation des matrices resultats

    Ex=zeros(2*largeur+1,2*largeur+1,2*longueur+1);
    Ey=zeros(2*largeur+1,2*largeur+1,2*longueur+1);
    Ez=zeros(2*largeur+1,2*largeur+1,2*longueur+1);
    for ix=1:2*largeur+1
        for ygrec=1:2*largeur+1
            for zed=1:2*longueur+1

                %calcul de \rho  et \phi en fonction de x et y
                %
                %                 if (round(sqrt((ix-largeur-1)*(ix-largeur-1)+(ygrec-largeur-1)*(ygrec-largeur-1)))>largeur)
                %                     res='meuh';
                %                 else
                rau=1+round(sqrt((ix-largeur-1)*(ix-largeur-1)+(ygrec-largeur-1)*(ygrec-largeur-1)));
                if(ix==(largeur+1))
                    if(ygrec>=(largeur+1))
                        fi=pi/2;
                    else
                        fi=-pi/2;
                    end
                else
                    fi=atan2((ygrec-largeur-1),(ix-largeur-1));
                end

                %Remplissage de la la matrice

                Ex(ix,ygrec,zed)=ctteH*(A2A(rau,zed)+2*A2B(rau,zed)-2*A00(rau,zed)-2*cos(2*fi)*(A2C(rau,zed)+A2D(rau,zed)+A02(rau,zed))+cos(4*fi)*A2E(rau,zed));
                Ey(ix,ygrec,zed)=ctteH*((A2C(rau,zed)+A02(rau,zed))*(-2)*sin(2*fi)+A2E(rau,zed)*sin(4*fi));
                Ez(ix,ygrec,zed)=ctteH*2*i*(A2G(rau,zed)*cos(3*fi)+cos(fi)*(4*A01(rau,zed)-3*A2F(rau,zed)));
                %
                %                     Ex(ix,ygrec,zed)=ctteH*(A2A(rau,zed)+2*A2B(rau,zed)-2*cos(2*fi)*(A2C(rau,zed)+A2D(rau,zed))+cos(4*fi)*A2E(rau,zed));
                %                     Ey(ix,ygrec,zed)=ctteH*((A2C(rau,zed))*(-2)*sin(2*fi)+A2E(rau,zed)*sin(4*fi));
                %                     Ez(ix,ygrec,zed)=ctteH*2*i*(A2G(rau,zed)*cos(3*fi)+cos(fi)*(-3*A2F(rau,zed)));
                %

            end
        end
    end
end

Itotint=sum(sum(sum(abs(Ex.*Ex)+abs(Ey.*Ey)+abs(Ez.*Ez))));
coef=sqrt(Itot/Itotint);
Ex=coef*Ex;
Ey=coef*Ey;
Ez=coef*Ez;

save Ex.mat Ex
save Ey.mat Ey
save Ez.mat Ez

Px=Ex.*Ex.*Ex+Ex.*Ey.*Ey+Ez.*Ez.*Ex;
save Px.mat Px
clear Px
Py=Ey.*Ey.*Ey+Ey.*Ex.*Ex+Ez.*Ez.*Ey;
save Py.mat Py
clear Py
Pz=Ez.*Ez.*Ez+Ey.*Ey.*Ez+Ex.*Ex.*Ez;
save Pz.mat Pz
clear Pz
clear Ex
clear Ey
clear Ez


'fin de l assemblage'
toc
