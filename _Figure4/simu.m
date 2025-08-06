% Programme simulant le signal THG cree une excitation specifique et une
% serie de geometries differentes.
%Argument: A (float vect),ecxitation (string)
%A(1) contient le nombre de mailles en longueur
%A(2) contient le nombre de mailles en largeur
%A(3) contient la largeur de la maille
%A(4) contient le longueur de la maille
%A(5) contient la largeur de la maille
%A(6) contient le pas d'integration chp lointain dx
%A(5) contient la nombre de pas d integration ne champ lointain
%Excitation decrit le mode d excitation

function [resultat] = simu(A,B,C,excitation,option,detection, geometry,nom)


%demarrage du chronometre
tic;
%initialisation des variables
resultat=[];
bigresx=[];
bigresy=[];
bigresz=[];
res=[];

%recuperations des entrees
longueur=A(1)
largeur=A(3)
b=A(4)*1e-9
bz=A(2)*1e-9
absz=0.01
dx=A(6)*absz
dy=dx;
pas=A(5)


%sauvegarde des parametres


NA=B(1);
lambda_1200=B(2);%1.2*10^(-6); % longueur d onde centrale d excitation
n1_400=B(3);%1.5308; %bk7
n1_1200=B(4);%1.5049;
f0=B(5);%2; %facteur de remplissage de l obsj
ctmat=[NA, lambda_1200, n1_400, n1_1200, f0];
save ctmat.mat ctmat


alpha=[C(1),C(3)];
save alpha.mat alpha
deph=[C(2), C(4)];
save deph.mat deph

'Simulation THG avec faisceau d excitation'
'Geometrie:'

%Calcul de l'intensite d'excitation au cube sur un maillage bm*bm*bm
%assembletout calcul et sauvegarde les matrices specifiees par option, en
%utilisant les matrices Aij calculees par le prog integretout

dirname=['./'  option '/' detection '/' excitation '/' nom]
mkdir(dirname)



assembletout(longueur,largeur,b,bz,excitation);
%fin de la premiere etape de calcul ( sim de l excitation)

toc

%les resultats de assembletout pertinents pour cette simulation sont charges

load Px.mat
load Py.mat
load Pz.mat


for oy=1:(geometry+1)
    'debut de l integration nouvelle geometrie'
    %initialisation des resultas intermediaires
    resx=[];
    resy=[];
    resz=[];

    if (strcmp(detection,'arriere'))
        z=-1*absz
    else
        z=absz
    end

    %deph=0;
    %Calcul et sauvegarde la geometrie de l echantillon
    if (strcmp(option,'homogene')==1)
        C=1;
        lim=1;
    end
    if (strcmp(option,'normalisation')==1)
        C=Interface(direction,longueur+1,largeur,longueur);
        lim=1;
    end

    if (strcmp(option,'custom')==1)
        if (oy==1)
            C=Interface('z',longueur+1,largeur,longueur);
        end
        if (oy==2)
            C=slab('z',9,largeur,longueur);
        end
        if (oy==3)
            C=slab('z',12,largeur,longueur);
        end
        if (oy==4)
            C=slab('z',15,largeur,longueur);
        end
        if (oy==5)
            C=slab('z',19,largeur,longueur);
        end
        if( oy>5)
            C=0;
        end
        lim=1;
    end
    if (strcmp(option,'interfacex')==1)
        direction='x';
        lim=(oy-1)*floor((2*largeur+1)/geometry);
        C=Interface(direction,lim,largeur,longueur);
    end

    if (strcmp(option,'interfacey')==1)
        direction='y';
        lim=(oy-1)*floor((2*largeur+1)/geometry);
        C=Interface(direction,lim,largeur,longueur);
    end

    if (strcmp(option,'interfacez')==1)
        direction='z';
        lim=(oy-1)*floor((longueur)/geometry);
        C=Interface(direction,longueur+1+lim,largeur,longueur);
    end

    if (strcmp(option,'slabx')==1)
        direction='x';
        lim=(oy-1)*floor(longueur/geometry);
        C=slab(direction,lim,largeur,longueur);
    end

    if (strcmp(option,'slaby')==1)
        direction='y';
        lim=(oy-1)*floor(longueur/geometry);
        C=slab(direction,lim,largeur,longueur);
    end

    if (strcmp(option,'slabz')==1)
        direction='z';
        lim=(oy-1)*floor(longueur/geometry);
        C=slab(direction,lim,largeur,longueur);
    end

    if (strcmp(option,'periode')==1)
        if(strcmp(detection,'arriere'))
            freq=2.33+0.6*(oy-1);
            lim=freq;
        end
        if(strcmp(detection,'avant'))
            freq=10+oy*4;
            lim=freq;
        end
        C=periode(lim,largeur,longueur);
    end

    if (strcmp(option,'angle')==1)
        alpha=(oy-1)*floor(180/geometry)*pi/180;
        lim=alpha;
        C=interfaceangle2(alpha,longueur,largeur);
    end


    savefile=[dirname '/param.mat']
    A=[A lim];
    save(savefile,'A')


    %Boucle de calcul du champ lointain
    for f=1:2*pas+1
        for g=1:2*pas+1

            %Valeurs arbitraires, a optimiser apres avoir reflechi un peu.

            x=dx*(g-pas-1);
            y=dy*(f-pas-1);


            er=sqrt(z*z+x*x+y*y);
            er2=er*er;
            M11=1-(x*x)/er2;
            M22=1-(y*y)/er2;
            M12=-x*y/er2;
            M13=-x*z/er2;
            M23=-y*z/er2;
            M33=1-(z*z)/er2;


            %Calcul du terme de propagation proche a lointain (matrice (2m+1)^3 )


            T = phase(x,y,z,largeur,longueur,b,bz);

            %integration
            T=T.*C;
            fa=facteur(x,y,z);

            S0=T.*(M11.*Px+M12*Py+M13.*Pz);
            res1=sum(sum(sum(S0)));
            res1=dx*dy*b*b*bz*res1*fa;
            res11=(abs(res1))^2;
            resx=[resx res11];

            S0=T.*(M22.*Py+M12.*Px+M23.*Pz);
            res1=sum(sum(sum(S0)));
            res1=dx*dy*b*b*bz*res1*fa;
            res11=(abs(res1))^2;
            resy=[resy res11];

            S0=T.*(M33.*Pz+M13.*Px+M23.*Py);
            res1=sum(sum(sum(S0)));
            res1=dx*dy*b*b*bz*res1*fa;
            res11=(abs(res1))^2;
            resz=[resz res11];

        end
        f
        toc
    end
    resx=reshape(resx,2*pas+1,2*pas+1);
    resy=reshape(resy,2*pas+1,2*pas+1);
    resz=reshape(resz,2*pas+1,2*pas+1);
    bigresx=[bigresx resx];
    bigresy=[bigresy resy];
    bigresz=[bigresz resz];
    integralex=dx*dy*sum(sum(resx));
    integraley=dx*dy*sum(sum(resy));
    integralez=dx*dy*sum(sum(resz));

    resultat=[resultat [integralex;integraley;integralez]]
    save resultat.mat resultat
    save bigresx.mat bigresx
    save bigresy.mat bigresy
    save bigresz.mat bigresz
    toc;
end

%Sauvegarde des donnees

mkdir(dirname)

savefile=[dirname '/res.mat']
save(savefile,'resultat')

savefile=[dirname '/resx.mat'];
save(savefile,'bigresx')

savefile=[dirname '/resy.mat'];
save(savefile,'bigresy')

savefile=[dirname '/resz.mat'];
save(savefile,'bigresz')

tarfilename=[ excitation '-' option '-' detection]
tar(tarfilename,dirname);

tar('conditions.tar','*.m');
copyfile('conditions.tar',dirname)

end
