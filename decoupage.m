function []=decoupage(matx,nom)
s=size(matx)
s2=s(2)/s(1);
s1=s(1);
for k=1:s2
    mat=matx(:,(k-1)*s1+1:k*s1);
    num=int2str(k);
    name=[nom,'_' num '.tif'];
    imwrite(uint16(mat),name,'tif','Compression','none');
end