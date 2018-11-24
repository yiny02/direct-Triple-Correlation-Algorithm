# direct Triple-Correlation Algorithm
1.	System Requirements: 
    All matlab scripts (.m) were developed and tested on Matlab 2016b (Windows 10 64-bit and Windows 7 64-bit);
    The smTriCorrCPU.mexw64 is a Matlab executable function and requires the OS to be win-64;
    The smTriCorrGPU.mexw64 is a Matlab executable function. It requires 1) the OS to be win-64 and 2) an NVDIA Graphic Card that supports CUDA 8.0 computing capability.
    Both .mexw64 files were tested on Matlab 2016b or later on Windows 7/10 64-bit OS.
2.	Installation guide: 
    Please follow the general installation guide of Matlab on https://www.mathworks.com/;
    Please follow the general installation guide of CUDA 8.0 on https://developer.nvidia.com/cuda-80-ga2-download-archive;
    Please note that the latest version of CUDA is CUDA 9.0 to-date, but make sure if the matlab you are about to install supports CUDA 9.0 capability.
3.	Demo: 
    A demo file TripleDemo.m file was provided for test and guidance of using other functions.
    A quick test of all functions could be to run the TripleDemo.m file as the main function in Matlab with given (default) parameters. 
    Please note that instead of providing small data set for test, we provided a simulation script for user to generate testing data set with different parameters. The default simulation, which is called in the TripleDemo.m main function, generates 50% R-G-B and 50% R-B-G patterns in a same canvas with randomized orientations. The output of this test, a variable named triple_trans.mat will generate a triple-correlation cube that being able to resolve such heterogenous patterns.
    The time cost of the core functions was given in Figure 2, with the RAM of 64GB DDR4, the CPU of Core i7 7800X, and a GPU of NVDIA GTX1060 (6GB) if applied. We note that when the image is extremely dense, RAM cost should be concerned.
4.	Introduction of use: 
    Please 1) prepare the reconstructed coordinates (of each of the three color-channel) into two columns with x in col1 and y in col2; 2) submit the three channel coordinates as well as the parameters to the functions following the instruction of Triple-Correlation computing section in the TripleDemo.m file. The default parameters were the ones used in a manuscript. Reproduction of the raw data given in that manuscript is available upon request to Eli.Rothenberg@nyumc.org, or Yandong.Yin@nyumc.org.
