function [D] = interfaceangle(alpha,largeur)
larg=2*largeur+1;
D=ones(larg,larg);
for x=1:larg
    for y=1:larg
        if ((x-largeur-1)<=(y-largeur-1)*tan(alpha))
            D(x,y)=0;
        end
    end
end
