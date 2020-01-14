################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../src/geometry/surface/BottomParametric.cu \
../src/geometry/surface/TopParametric.cu \
../src/geometry/surface/bicubicsrf.cu \
../src/geometry/surface/compositeBicubicsrf.cu \
../src/geometry/surface/mySurface.cu 

OBJS += \
./src/geometry/surface/BottomParametric.o \
./src/geometry/surface/TopParametric.o \
./src/geometry/surface/bicubicsrf.o \
./src/geometry/surface/compositeBicubicsrf.o \
./src/geometry/surface/mySurface.o 

CU_DEPS += \
./src/geometry/surface/BottomParametric.d \
./src/geometry/surface/TopParametric.d \
./src/geometry/surface/bicubicsrf.d \
./src/geometry/surface/compositeBicubicsrf.d \
./src/geometry/surface/mySurface.d 


# Each subdirectory must supply rules for building sources it contributes
src/geometry/surface/%.o: ../src/geometry/surface/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-9.2/bin/nvcc -G -g -O0 -Xcompiler -dc -std=c++11 -gencode arch=compute_50,code=sm_50  -odir "src/geometry/surface" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-9.2/bin/nvcc  -G -g -O0 -Xcompiler -dc -std=c++11 --compile --relocatable-device-code=true -gencode arch=compute_50,code=compute_50 -gencode arch=compute_50,code=sm_50  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


