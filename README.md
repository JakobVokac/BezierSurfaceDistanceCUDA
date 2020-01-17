# BezierSurfaceDistanceCUDA

This is the CUDA translation of the distance computation algorithm for parametric and bezier surfaces (see [BezierSurfaceDistance](https://github.com/JakobVokac/BezierSurfaceDistance)), mean't to calculate distances for the aortic valve model.

The project is setup in Nsight Eclipse Edition. Unfortunately, the IDE does not provide an option to export the project as a CMake project, so the IDE is required to compile the code.

In order to compile the program for a specific GPU, the user must know the compute capability of that GPU. Using Nsight Eclipse Edition, the build can then be set for the specific compute capability under: 

**Project -> Properties -> Build -> Settings**

The IDE builds an executable which can then be run on other devices with the help of **nvcc**.

The current built executable in the Debug folder of the project can be run with the **nvcc** compiler or the **nvprof** profiler on a GPU with the compute architecture 5.0 (see [CUDA Programming Guide](https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#compute-capabilities) and [NVCC compiler](https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#introduction)).



The executable of the project is setup to take input from the **input.txt** file, which should contain numbers in triples (3 in a line, for each input point) and nothing else.
