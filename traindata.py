import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
from skimage.draw import polygon
import numpy as np
import random
S=512
N=8 # number of sides of polygon
R=8 # max raidus of polygon
nt=100
nd=10
nump=25
numpb=6
zx=np.zeros((nt,S,S,nd))
zy=np.zeros((nt,S,S,nd))
for idx1 in range(0,nt):
    bsigma=random.random()*1.5
    for idx2 in range(0,nd):
        zy_temp=np.zeros((S,S))
        n=np.random.randint(nump)+numpb
        for idx3 in range(0,n):
            z_temp=np.zeros((R*2,R*2),'float32')
            #: fill polygon
            a=np.sort(np.random.rand(N))*2*np.pi
            r=np.random.randint(4,R-1,(N,1))
            x=np.round(np.cos(a)*r+R)
            y=np.round(np.sin(a)*r+R)
            rr,cc=polygon(x,y,z_temp.shape)
            z_temp[rr,cc]=1
            posx=np.random.randint(S-R*2-1)+R
            posy=np.random.randint(S-R*2-1)+R
            zy_temp[posx-R:posx+R,posy-R:posy+R]=zy_temp[posx-R:posx+R,posy-R:posy+R]+z_temp
        zy_temp[zy_temp>0]=1
        zy[idx1,:,:,idx2]=ndimage.gaussian_filter(zy_temp,(bsigma,bsigma))
        zx[idx1,:,:,idx2]=ndimage.gaussian_filter(zy_temp,(bsigma,bsigma))
for idx1 in range(0,nt):
    for idx2 in range(0,nd):
        for idx3 in range(0,nd):
            if (idx2!=idx3):
                dis=np.min([np.abs(idx2-idx3),np.abs(idx2-idx3+10),np.abs(idx2-idx3-10)])
                sigma=dis*bsigma+bsigma
                temp=ndimage.gaussian_filter(zy[idx1,:,:,idx3],(sigma,sigma))
                zx[idx1,:,:,idx2]=zx[idx1,:,:,idx2]+temp
        zx[idx1,:,:,idx2]=(zx[idx1,:,:,idx2]-np.mean(zx[idx1,:,:,idx2]))/np.std(zx[idx1,:,:,idx2])
zx=np.swapaxes(zx,0,2)
zy=np.swapaxes(zy,0,2)
zy=np.reshape(zy,(512,512,1000),order='F')
zx=np.reshape(zx,(512,512,1000),order='F')
zy=np.swapaxes(zy,0,2)
zx=np.swapaxes(zx,0,2)
plt.figure()
plt.imshow(zy[:,:,3])
plt.show()
plt.figure()
plt.imshow(zx[:,:,3])
plt.show()
