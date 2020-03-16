# BING++
[Sequential Optimization for Efficient High-Quality Object Proposal Generation](https://zimingzhang.wordpress.com/publications/)

This source code is free for academic usage under LICENSE.

Tested on Ubuntu 14.04 

Need to pass OPENCV_PATH which includes lib and include folders to cmake as instructed below.

To build:
```
mkdir build
cd build
cmake .. -DOPENCV_PATH=/path/to/opencv/ (e.g. /home/user/opencv3.1.0/)
make
```

To run:
```
./BINGpp/BING++ /path/to/data/ (e.g. /datasets/VOC2007/)
```

Notes:

Annotations.tar.gz has to be extracted into VOC2007 folder.

JPEGImages folder from VOC2007 dataset is not included due to its size.

Included vlfeat for 64bit linux in ext, if you are using another arch please add appropriate libvl or edit CMakelists.txt in Objectness.

