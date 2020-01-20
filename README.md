# BezierSurfaceDistanceCUDA

This is the CUDA translation of the distance computation algorithm for parametric and bezier surfaces (see [BezierSurfaceDistance](https://github.com/JakobVokac/BezierSurfaceDistance)), mean't to calculate distances for the aortic valve model.

The project is setup in Nsight Eclipse Edition. Unfortunately, the IDE does not provide an option to export the project as a CMake project, so the IDE is required to compile the code.

To import the project, follow the steps below, taken (and slightly modified) from [Nsight Eclipse Plugins Edition section 2.11](https://docs.nvidia.com/cuda/nsight-eclipse-plugins-guide/index.html#import-nsight-project):

1.  Create a new C++ Project in your workspace in Nsight Eclipse edition.
2.  Specify the project name and choose Empty project type with CUDA toolchains.
3.  Right click on the project to import the source files. **Import > General > File System >(From directory)** or copy the source files from the existing project, this should include only the src files, not the Debug folder with the given executable.
4.  Import the project settings like include paths and symbols using the following right click menu **Import > C/C++ > C/C++ Project Settings >Next...** and import the **projSettings.xml** file.
5.  Select the location of the project settigns file and select the project and configuration on the next wizard page.
6.  Complete the wizard. The project settings will be imported from the file exported from Nsight Eclipse Edition.
7.  Build the project by clicking on the hammer button on the main toolbar.
8.  (Optional) You might have to create a new configuration for the executable:
  - Once the project is built, check the Debug folder for the executable, should have the same name as your project.
  - Go under **Run > Run Configurations...** and select the executable for your configuration.
  - Select any options you might need and press **Apply**. You should now be able to run the program.

In order to compile the program for a specific GPU, the user must know the compute capability of that GPU. Using Nsight Eclipse Edition, the build can then be set for the specific compute capability under: 

**Project -> Properties -> Build -> Settings** or **File -> Properties -> Build -> Settings**

The IDE builds an executable which can then be run on other devices with the help of **nvcc**.

The current built executable in the Debug folder of the project can be run with the **nvcc** compiler or the **nvprof** profiler on a GPU with the compute architecture 5.0 (see [CUDA Programming Guide](https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#compute-capabilities) and [NVCC compiler](https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#introduction)).



The executable of the project is setup to take input from the **input.txt** file, which should contain numbers in triples (3 in a line, for each input point) and nothing else.
