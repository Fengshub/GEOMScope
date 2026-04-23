%%
clear;clc;close all;
%%
load('fluorescent_particles_measurement_uint8.mat'); % load measurement data
im=double(im)/255;
figure
imagesc(imgaussfilt(im,3))
daspect([1 1 1])
title('measurement')
%% for point source reconstruction, the rotation angle is very important and needs to be very precise!
% load coordinates of lens units in (-10mm,10mm)
lenscoordinates;
mask=zeros(5120);
% convert to measurement matrix coordinate
a=a+11.52;
a=round(a*222.2)+1;
for idx=1:length(a)
    mask(a(idx,1),a(idx,2))=1;
end
% calibrate rotation angles
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
mask2=padarray(mask2,[500 500]);
mask2=mask2(round(size(mask2,1)/2)-2559+50:round(size(mask2,1)/2)+2560+50,round(size(mask2,2)/2)-2559-50:round(size(mask2,2)/2)+2560-50);
psf=mask2;
figure
imagesc(imgaussfilt(psf,5))
title('psf')
daspect([1 1 1])
%% 10/21 11/20
% close all
S=1024; % reconstruction grid size
fov=S;
recon_rec=zeros(S,S,10); % allocate reconstruction volume
So=5; % downsampling rate
ridx=1
orx=2560;
ory=2560;
disi=4;% image plane distance
[cx0,cy0]=find(psf~=0);% lens unit coordinates
for disd=25:0.5:28
    disp(['reconstructing plane at distance ',num2str(disd),' mm'])
    mag=disi/disd;% image/object magnification
    scale=1+mag;
    smin=1;
    smax=1024;
    recon=zeros(smax-smin+1);
    % scan voxel coordinates and back propagate image pixels
	for idx1=smin:smax
        iidx1=idx1-smin+1;
        cx=round((cx0-idx1*So)*scale+idx1*So);
        for idx2=smin:smax
            iidx2=idx2-smin+1;
            cy=round((cy0-idx2*So)*scale+idx2*So);% find contributing image pixels for current voxel as (cx, cy)
            list=zeros(1,length(cx));
            for idx3=1:length(cx)
                sx=(cx(idx3)-idx1*So)*1;
                sy=(cy(idx3)-idx2*So)*1;
                dis=sqrt(sx^2+sy^2); % lens center to reconstruction point lateral distance
                if (cy(idx3)>1) && (cx(idx3)>1) && (cx(idx3)<=5120) && (cy(idx3)<=5120) &&(dis<1500) % selecting effective image pixels
                    recon(iidx1,iidx2)=recon(iidx1,iidx2)+im(floor(cx(idx3)),floor(cy(idx3)));% pixel back propagation
                end
            end
        end
    end
    recon_rec(:,:,ridx)=imresize(recon,[S,S]);
    ridx=ridx+1;
end
imagesc(recon)
%%
close all
for idx=1:ridx-1 % plot all reconstructed matrix slices
    figure
    imagesc(recon_rec(:,:,idx))
    daspect([1 1 1])
    title(['reconstruction plane ',num2str(25+0.5*(idx-1)),' mm'])
end
