function res= traitementdatafinal(date,geo,dir,opt)
directory=['./' date '/' geo '/' dir ];
copyfile('decoupage.m',directory)
cd(directory)
N=ls;
s=size(N);
m=s(1)-3;
for k=1:m
    name=N(2+k,:);
    director=['./' name];
    copyfile('decoupage.m',director)
    cd(director)
    director=['./' opt];
    copyfile('decoupage.m',director)
    cd(director)
    load res.mat;
    Rx=resultat(1,:);
    Ry=resultat(2,:);
    Rz=resultat(3,:);
    Rs=resultat(1,:)+resultat(2,:)+resultat(3,:);
    if (k==1)
        load param.mat
        sze=size(A)-6;
        load res.mat
        s=size(resultat)
        sze=s(2);
        vect=[1:sze];
		if (strcmp(geo,'interfacez'))
           % vect=2*(A(4)*A(1)/A(5))*[-(A(5)/2):(A(5)/2)];
		  vect= [1:sze];
        end

        if (strcmp(geo,'interfacex'))
          %vect=2*(A(3)*A(2)/A(5))*[-(A(5)/2):(A(5)/2)];
          vect=[1:sze];
        end
        if (strcmp(geo,'slabx'))
            %vect=(A(3)*A(2)/A(5))*[0:A(5)];
            vect=[1:sze];
        end
        if (strcmp(geo,'slaby'))
            %vect=(A(3)*A(2)/A(5))*[0:A(5)];
            vect=[1:sze];
        end
        if (strcmp(geo,'slabz'))
            vect=[1:sze];
        end
        if (strcmp(geo,'periode'))
            vect=[1:sze];
        end
        if (strcmp(geo,'angle'))
            vect=(180/(sze-1))*[0:sze-1];
        end
        res=[vect];
        resx=[vect];
        resy=[vect];
        resz=[vect];
    end
    %% 
    
    res=[res;Rs]
    resx=[resx; Rx]
    resy=[resy; Ry]
    resz=[resz; Rz]
    load resx
    bigresx=16000*bigresx/(max(max(bigresx)));
    decoupage(bigresx,'em');
    cd ..
    cd ..
end

cd ..
cd ..
cd ..

res=res.';
resx=resx.';
resy=resy.';
resz=resz.';

savefile=['./' date '_' geo '-' dir '.mat'];
save(savefile,'res', '-ASCII');
savefile=['./' date '_' geo '-' dir 'x.mat'];
save(savefile,'resx', '-ASCII');
savefile=['./' date '_' geo '-' dir 'y.mat'];
save(savefile,'resy', '-ASCII');
savefile=['./' date '_' geo '-' dir 'z.mat'];
save(savefile,'resz', '-ASCII');

figure(1);
plot(res(:,1),res(:,2:end))
figure(2);
plot(res(:,1),res(:,2:end)./max(res(:,2:end)))