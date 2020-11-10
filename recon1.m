%% 
name='dao7'
t=Tiff([name,'.tif'],'r'); % target
im=read(t);
im=double(im);
im(1373,2318)=0;
for iter=1:2
    iter
    im_f=(fftshift((fft2(ifftshift(im)))));
    off=0;
    im_f(2561-off:2561+off,2561-off:2561+off)=0;
    im2=(abs(ifftshift(ifft2(fftshift(im_f)))));
%     mse(im,im2)
    im=im2;
end
figure
imagesc(im)
daspect([1 1 1])
im=im-mean(im(:))-std(im(:));
im(im<0)=0;
%% 
for idx=1:10
    hpoly=drawpolygon(gca);
    maskp=createMask(hpoly);
    im(maskp)=0;
end
%% 10/21 10/25
load imdao3
%% 
psfd=zeros(7);
psfd(4,4)=1;
psfd=imgaussfilt(psfd,1.5)
J=deconvlucy(im,psfd,10);
% J=deconvwnr(im,psf,10);
figure
imagesc(J)
daspect([1 1 1])
im=J;
save(['im',name],'im')
%% 
load cursor_info_ethan1
mask=zeros(5120,5120,3);
for idx=1:length(cursor_info_ethan1)
    mask(cursor_info_ethan1(idx).Position(2),cursor_info_ethan1(idx).Position(1),1)=1;
end
mask2=mask;
mask2=padarray(mask2,[100 100]);
offx=50;offy=-50;
mask2=mask2(round(size(mask2,1)/2)-2559-offx:round(size(mask2,1)/2)+2560-offx,round(size(mask2,2)/2)-2559-offy:round(size(mask2,2)/2)+2560-offy);
psf=mask2;
figure
imagesc(imgaussfilt(psf,5))
daspect([1 1 1])
%% 10/21 10/25
%for point source reconstruction, the rotation angle is very important and needs to be very precise!
sup1_2;
mask=zeros(5120);
a=a+11.52;
a=round(a*222.2)+1;
for idx=1:length(a)
    mask(a(idx,1),a(idx,2))=1;
end
figure
% mask=imgaussfilt(mask,5);
roa=238.63;
roar=roa*pi/180;
mask2=zeros(size(mask));
[cx0,cy0]=find(mask~=0);
x=cx0-2560;
y=cy0-2560;
cx=zeros(size(cx0));
cy=zeros(size(cx0));
for idx3=1:length(cx0)
    cx(idx3)=cos(roar)*x(idx3)-sin(roar)*y(idx3)+2560;
    cy(idx3)=sin(roar)*x(idx3)+cos(roar)*y(idx3)+2560;
    if (cx(idx3)>0 & cx(idx3)<5120 & cy(idx3)>0 & cy(idx3)<5120)
        mask2(round(cx(idx3)),round(cy(idx3)))=1;
    end
end

% mask2=imrotate(mask,-52.1);
% mask2=mask2+rand(size(mask2))*1e-4;
% mask2m=double(imregionalmax(mask2));
mask2=padarray(mask2,[500 500]);
% mask2=mask2(round(size(mask2,1)/2)-2559:round(size(mask2,1)/2)+2560,round(size(mask2,2)/2)-2559:round(size(mask2,2)/2)+2560);
mask2=mask2(round(size(mask2,1)/2)-2559+50:round(size(mask2,1)/2)+2560+50,round(size(mask2,2)/2)-2559-50:round(size(mask2,2)/2)+2560-50);
[cx0,cy0]=find(mask2~=0);
psf=mask2;
imagesc(imgaussfilt(psf,5))
daspect([1 1 1])
%% 
orx=2560;
ory=2560;
scale=1.1295;
[cx0,cy0]=find(psf~=0);
for idx3=1:length(cx0)
    cx=round((orx-cx0)*(scale-1)+cx0);
    cy=round((ory-cy0)*(scale-1)+cy0);
end
psf2=zeros(size(psf));
for idx3=1:length(cx0)
    psf2(cx(idx3),cy(idx3))=1;
end
psf=psf2;
figure
imagesc(imgaussfilt(psf,5))
daspect([1 1 1])
%% 10/21
% close all
S=1024;
fov=S;
recon_rec=zeros(S,S,10);
So=5;
ridx=1
figure(1)
imagesc(im)
orx=2560;
ory=2560;
% t1=1;
% t2=1;
scale0=1.1295;
mag0=scale0-1;
for disd=32:0.4:36.6
    mag=4.6/disd;
    scale=1+mag;
% for idx=-2:2:2
%     scale=scale0+0.005*idx;
%     mag=scale-1;
    smin=1;%round(S/2-fov/2*mag0/mag);
    smax=1024;%round(S/2+fov/2*mag0/mag);
    recon=zeros(smax-smin+1);
    [cx0,cy0]=find(psf~=0);
	for idx1=smin:smax
        iidx1=idx1-smin+1;
        cx=round((cx0-idx1*So)*scale+idx1*So);
        for idx2=smin:smax
            iidx2=idx2-smin+1;
%             figure(1)
%             hold on
            cy=round((cy0-idx2*So)*scale+idx2*So);
            list=zeros(1,length(cx));
            for idx3=1:length(cx)
                sx=(cx(idx3)-(idx1*So-orx)*mag)-idx1*So; % shift at 5120 coordinates, lens center to reconstruction point distance lateral
                sy=(cy(idx3)-(idx2*So-ory)*mag)-idx2*So;
                dis=sqrt(sx^2+sy^2); % distance at 5120 coordinates
                if (cy(idx3)>1) && (cx(idx3)>1) && (cx(idx3)<=5120) && (cy(idx3)<=5120) &&(dis<1500)
%                     list(idx3)=im(floor(cx(idx3)),floor(cy(idx3)));
                    recon(iidx1,iidx2)=recon(iidx1,iidx2)+im(floor(cx(idx3)),floor(cy(idx3)));
%                     if im(round(idx1*So+sx),round(idx2*So+sy))>20
%                         plot(cy(idx3),cx(idx3),'go')
%                     else
%                         plot(cy(idx3),cx(idx3),'ro')
%                     end
%                 else
%                     plot(cy(idx3),cx(idx3),'mo')
                end
            end
%             lista=list(list~=0);
%             listamean=mean(lista);
%             listastd=std(lista);
%             listc=lista((lista>listamean-2.5*listastd)&(lista<listamean+2.5*listastd));
%             recon(iidx1,iidx2)=recon(iidx1,iidx2)+mean(listc);
%             if (idx1==t1) && (idx2==t2)
%                 figure(1)
%                 hold on
%                 cy=round((cy0-idx2*So)*scale+idx2*So);
%                 for idx3=1:length(cx)
%                     sx=(cx(idx3)-(idx1*So-orx)*mag)-idx1*So; % shift at 5120 coordinates, lens center to reconstruction point distance lateral
%                     sy=(cy(idx3)-(idx2*So-ory)*mag)-idx2*So;
%                     dis=sqrt(sx^2+sy^2); % distance at 5120 coordinates
%                     if (cy(idx3)>1) && (cx(idx3)>1) && (cx(idx3)<=5120) && (cy(idx3)<=5120) &&(dis<1300)
% %                         recon(iidx1,iidx2)=recon(iidx1,iidx2)+im(floor(cx(idx3)),floor(cy(idx3)));
%                         if im(round(idx1*So+sx),round(idx2*So+sy))>20
%                             plot(cy(idx3),cx(idx3),'go')
%                         else
%                             plot(cy(idx3),cx(idx3),'ro')
%                         end
%                     else
%                         plot(cy(idx3),cx(idx3),'mo')
%                     end
%                 end
%             end
        end
    end
%     recon=imrotate(recon,180);
    recon_rec(:,:,ridx)=imresize(recon,[S,S]);
    ridx=ridx+1
end
figure
imagesc(recon)
% hold on
% plot(t2,t1,'o')
daspect([1 1 1])
% caxis([1000 3000])
%% 
% c=zeros(S);
[x,y]=meshgrid(linspace(-1,1,S));
c=x.^2+y.^2+5;
c=c./max(c(:));
for idx=1:size(recon_rec3,3)
    recon_rec3(:,:,idx)=recon_rec3(:,:,idx).*c;
end
%% 
recon_rec3=recon_rec;
%% 
% recon_rec(isnan(recon_rec))=0;
% recon_rec3=recon_rec;
recon_rec=recon_rec3;
recon_rec(recon_rec<400)=0;
% recon_rec=recon_rec3;
% recon_rec=imgaussfilt(recon_rec,3);
% CCmap=zeros(size(recon_rec));
% for idx=1:size(recon_rec,3)
%     CCidx=bwconncomp(recon_rec(:,:,idx));
%     temp=CCmap(:,:,idx);
%     for cidx=1:CCidx.NumObjects
%         temp(CCidx.PixelIdxList{cidx})=1;
%     end
%     CCmap(:,:,idx)=temp;

CC = bwconncomp(recon_rec,26);
% recon_rec2=zeros(size(recon_rec));
% for idx=1:CC.NumObjects
%     temp=CC.PixelIdxList{idx};
%     a=find(recon_rec==max(recon_rec(temp)));
%     recon_rec2(a)=recon_rec(a);
% end
recon_rec5=zeros(size(recon_rec));
for th=0.05:0.05:0.2
    for idx=1:CC.NumObjects
        temp=CC.PixelIdxList{idx};
        recon_rec(temp)=recon_rec(temp)/max(recon_rec(temp));
    end
    recon_rec(recon_rec<th)=0;

    CC = bwconncomp(recon_rec,26)
    recon_rec4=zeros(size(recon_rec));
    for idx=1:CC.NumObjects
        temp=CC.PixelIdxList{idx};
        recon_rec(temp)=recon_rec(temp)+rand(size(temp))*0.005;
        a=find(recon_rec==max(recon_rec(temp)));
        recon_rec4(a)=recon_rec3(a);
        recon_rec5(a)=recon_rec3(a);
    end
end
%% 
recon=sum(recon_rec5(:,:,:),3);
%%
close all

for idx=1:8
    figure
    imagesc(recon_rec3(:,:,idx))
    daspect([1 1 1])
    axis off
end

%% 
% recon2=flip(recon,2);
figure
imagesc(imgaussfilt(recon,2))
% imagesc(recon2)
daspect([1 1 1])
caxis([0 100])
%% 10/21
close all
for idx=1:8
    figure
    imagesc(imgaussfilt(recon_rec4(:,:,idx),4))
    daspect([1 1 1])
    axis off
    caxis([0 100])
end
%%
maxi=max(recon_rec4(:));
maxi=0.1;
for idx=1:10
    temp=imgaussfilt(recon_rec4(:,:,idx),3);
    temp=temp/maxi;
    imwrite(temp,[num2str(idx),'.jpg'])
end
