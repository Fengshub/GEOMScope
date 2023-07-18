# GEOMScope: Large Field-of-View 3D Lensless Microscopy with Low Computational Complexity
This is a temporary repository for work <br />
Tian, F., Hu, J., Yang, W., GEOMScope: Large Field-of-View 3D Lensless Microscopy with Low Computational Complexity. Laser & Photonics Reviews 2021, 15, 2100072.
## [Paper](https://onlinelibrary.wiley.com/doi/full/10.1002/lpor.202100072)
![Diagram of geomscope](https://github.com/Fengshub/GEOMScope/blob/master/imgs/imager_schematic2.png)
### Clone this repository:
```
git clone https://github.com/Yang-Research-Laboratory/Learned-Lensless-3D-Camera.git
```

## Ray tracing algorithm and clustering algorithm for sparse objects
We provided the [**code**](https://github.com/Fengshub/GEOMScope/blob/master/ray_tracing.m) for 3D object reconstruction from measured image with ray tracing. The ray tracing algorithm is to generate for intermediate 3D reconstruction results before post processing with CNN or clustering algorithms. 
The raw measurement data example is provided in [**folder**](https://github.com/Fengshub/GEOMScope/blob/master/data/rawdata.mat), with camera TIF image in [**drive folder**](https://drive.google.com/file/d/1g8AMWo8S9XkrMhYhekwzK4__OOIc8-9L/view?usp=drive_link). <br />
The algorithms includes data loading, PSF initialization, ray tracing and multi-threshold clustering. <br />
Our random microlens array has each lens unit [**coordinates**](https://github.com/Fengshub/GEOMScope/blob/master/acoordinate.m) recorded in mm units, the coordinates are used to initialize point spread function (PSF).

## CAD files and macros
We provide CAD files of designed microlens array [**mold**](https://drive.google.com/file/d/1RoEhVsZKAH5rszAVYK7zxNvhVvGT_JkV/view?usp=sharing) for 3D printing.
We also provide lens array & emission filter holder we used, and several [**CAD**](https://github.com/Fengshub/GEOMScope/tree/master/CAD) molds for imaging targets.

![postprocessing algorithm](https://github.com/Fengshub/GEOMScope/blob/master/imgs/postprocessing%20algorithm.PNG)
## UNet model for post processing
We provide [**code**](https://github.com/Fengshub/GEOMScope/blob/master/unet.py) to build convolutional network modules for post processing, and [**training**](https://github.com/Fengshub/GEOMScope/blob/master/traindata.py).

## Zemax Ray tracing simulation
We provide ZemaxOpticStudio lens [**zmx**](https://github.com/Fengshub/GEOMScope/blob/master/Zemax/lensarray_nsq.zmx) file in non sequential mode on the full microlens array simulation with our [**CAD target**](https://github.com/Fengshub/GEOMScope/blob/master/CAD/letter_v4.step).
