%%
clear;clc;close all;
%%
load('fig1_measurement_uint8.mat'); % load measurement data
im=double(im)/255;
figure
imagesc(imgaussfilt(im,3))
daspect([1 1 1])
title('measurement')
%% 
% load coordinates of lens units in (-10mm,10mm)
lenscoordinates;
mask=zeros(3072);
% convert to measurement matrix coordinate
a=a+10;
a=round(a*153.6)+1; 
for idx=1:length(a)
    mask(a(idx,1),a(idx,2))=1;
end
psf=mask;
figure
imagesc(imgaussfilt(psf,5))
title('psf')
daspect([1 1 1])
%% 
close all
S=1024; % reconstruction grid size
recon_rec=zeros(S,S,21); % allocate reconstruction volume
So=3; % downsampling rate
ridx=1;
disi=4; % image plane distance
[cx0,cy0]=find(psf~=0); % lens unit coordinates
for disd=17:0.5:23 % scan object distance
    disp(['reconstructing plane at distance ',num2str(disd),' mm'])
    recon=zeros(S);
    mag=disi/disd; % image/object magnification
    scale=1+mag;
    % scan voxel coordinates and back propagate image pixels
    for idx1=1:S
        cx=round((cx0-idx1*So)*scale+idx1*So); 
        for idx2=1:S
            cy=round((cy0-idx2*So)*scale+idx2*So); % find contributing image pixels for current voxel as (cx, cy)
            for idx3=1:length(cx)
                sx=(cx(idx3)-idx1*So)*1;
                sy=(cy(idx3)-idx2*So)*1;
                dis=sqrt(sx^2+sy^2); % lens center to reconstruction point lateral distance
                if (dis<1300) && (idx2*So+sy>1) && (idx1*So+sx>1) && (idx1*So+sx<S*So) && (idx2*So+sy<S*So) % selecting effective image pixels
                    recon(idx1,idx2)=recon(idx1,idx2)+im(round(idx1*So+sx),round(idx2*So+sy)); % pixel back propagation
                end
            end
        end
    end
    recon_rec(:,:,ridx)=recon;
    ridx=ridx+1;
end
%%
for idx=1:ridx-1 % plot all reconstructed matrix slices
    figure
    imagesc(recon_rec(:,:,idx))
    daspect([1 1 1])
    title(['reconstruction plane ',num2str(17+0.5*(idx-1)),' mm'])
end
