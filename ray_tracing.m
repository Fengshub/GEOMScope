%% read measurement image
name='Newmeasurement';
t=Tiff([name,'.tif'],'r'); % target image
im=read(t);
im=double(im);
% Frequenc domain filtering
for iter=1:2 % number of iterations for 
    im_f=(fftshift((fft2(ifftshift(im)))));
    off=0;
    im_f(2561-off:2561+off,2561-off:2561+off)=0;
    im2=(abs(ifftshift(ifft2(fftshift(im_f)))));
    im=im2;
end
figure
imagesc(im)
daspect([1 1 1])
%% load processed image
load rawdata.mat
%% R-L deconvolution to deblur (test version)
psfd=zeros(7);
psfd(4,4)=1;
psfd=imgaussfilt(psfd,1.5);
J=deconvlucy(im,psfd,10);
% J=deconvwnr(im,psf,10);
figure
imagesc(J)
daspect([1 1 1])
im=J;
save(['im',name],'im')
%% generate psf (for sensor pixel 5120 x 5120)
%for point source reconstruction, the rotation angle is very important and needs to be very precise!
acoordinate; % load the designed lens unit coordinates (in mm)
mask=zeros(5120); % psf pattern (ideal assumption)
a=a+11.52; % assign metric coordinate value to pixel coordinate
a=round(a*222.2)+1;
for idx=1:length(a)
    mask(a(idx,1),a(idx,2))=1; % generate psf
end
figure
% mask=imgaussfilt(mask,5);
roa=238.63; % rotation angle after calibration of lens array tilting to sensor
roar=roa*pi/180; % rotation angle in degrees
mask2=zeros(size(mask));
[cx0,cy0]=find(mask~=0);
x=cx0-2560;
y=cy0-2560;
cx=zeros(size(cx0));
cy=zeros(size(cx0));
for idx3=1:length(cx0) % rotation of psf peaks (can be replaced by rotational matrix / imrotate with proper interpolation
    cx(idx3)=cos(roar)*x(idx3)-sin(roar)*y(idx3)+2560;
    cy(idx3)=sin(roar)*x(idx3)+cos(roar)*y(idx3)+2560;
    if (cx(idx3)>0 & cx(idx3)<5120 & cy(idx3)>0 & cy(idx3)<5120)
        mask2(round(cx(idx3)),round(cy(idx3)))=1;
    end
end

mask2=padarray(mask2,[500 500]); % adjust lateral shift of lens array to image sensor
mask2=mask2(round(size(mask2,1)/2)-2559+50:round(size(mask2,1)/2)+2560+50,round(size(mask2,2)/2)-2559-50:round(size(mask2,2)/2)+2560-50);
[cx0,cy0]=find(mask2~=0);
psf=mask2;
imagesc(imgaussfilt(psf,5))
daspect([1 1 1])
%% calculate psf peak coordinate scaled at base distance for object reconstruction
orx=2560; % origin of pixel coordinate (5120 in both x and y)
ory=2560;
scale=1.1295; % calculated ratio between object-sensor distance to object-MLA distance
[cx0,cy0]=find(psf~=0);
for idx3=1:length(cx0) % calculate scaled PSF peak coordinates
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
%% 3D reconstruction by ray tracing
% close all
S=1024; % reconstruction voxel in x and y (same as sensor size)
fov=S; % reconstruction field of view in voxels (can be larger than sensor size)
recon_rec=zeros(S,S,10); % 3D reconstruction stack
So=5; % ration between voxel size and pixel size
ridx=1 % reconstruction distance index
figure(1)
imagesc(im)
orx=2560;
ory=2560;
% t1=1;
% t2=1;
scale0=1.1295; % base obj plane scale
mag0=scale0-1; % base obj plane magnification 
for disd=32:0.4:36.6 % reconstruction distance
    mag=4.6/disd; % calculate magnification
    scale=1+mag; % calculate scale
% for idx=-2:2:2
%     scale=scale0+0.005*idx;
%     mag=scale-1;
    smin=1;%round(S/2-fov/2*mag0/mag);  % min coordinate in voxel (equal to 1 when reconstruction fov = sensor size)
    smax=1024;%round(S/2+fov/2*mag0/mag); % max coordinate in voxel
    recon=zeros(smax-smin+1); % initialize current plan voxels
    [cx0,cy0]=find(psf~=0);
	for idx1=smin:smax % ray tracing scan voxels along x
        iidx1=idx1-smin+1;
        cx=round((cx0-idx1*So)*scale+idx1*So);
        for idx2=smin:smax % ray tracing scan voxels along y
            iidx2=idx2-smin+1;
%             figure(1)
%             hold on
            cy=round((cy0-idx2*So)*scale+idx2*So);
            list=zeros(1,length(cx));
            for idx3=1:length(cx) % % ray tracing scan effective imaging lenses
                sx=(cx(idx3)-(idx1*So-orx)*mag)-idx1*So; % shift at 5120 coordinates, lens center to reconstruction point distance lateral
                sy=(cy(idx3)-(idx2*So-ory)*mag)-idx2*So;
                dis=sqrt(sx^2+sy^2); % distance at 5120 coordinates
                if (cy(idx3)>1) && (cx(idx3)>1) && (cx(idx3)<=5120) && (cy(idx3)<=5120) &&(dis<1500) % identify valid back projections
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
            % ------------------for debug sampling voxel
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
    recon_rec(:,:,ridx)=imresize(recon,[S,S]); % record current depth reconstruction
    ridx=ridx+1 % move to next depth
end
figure
imagesc(recon)
% hold on
% plot(t2,t1,'o')
daspect([1 1 1])
% caxis([1000 3000])
%% clustering algorithm backup
recon_rec3=recon_rec;
%% multi-level custering algorithm (for sparse objects)
recon_rec=recon_rec3;

CC = bwconncomp(recon_rec,26); % calculate connected map

recon_rec5=zeros(size(recon_rec));
for th=0.05:0.05:0.2 % thresholding / separating beads in same clusters by relative intensity level
    for idx=1:CC.NumObjects
        temp=CC.PixelIdxList{idx};
        recon_rec(temp)=recon_rec(temp)/max(recon_rec(temp));
    end
    recon_rec(recon_rec<th)=0;

    CC = bwconncomp(recon_rec,26) % custering the thresholded reconstruction stack
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
%% plot clustering 3D stack
close all
for idx=1:8
    figure
    imagesc(recon_rec3(:,:,idx))
    daspect([1 1 1])
    axis off
end
