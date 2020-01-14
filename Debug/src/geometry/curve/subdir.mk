################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../src/geometry/curve/cubiccrv.cu 

OBJS += \
./src/geometry/curve/cubiccrv.o 

CU_DEPS += \
./src/geometry/curve/cubiccrv.d 


# Each subdirectory must supply rules for building sources it contributes
src/geometry/curve/%.o: ../src/geometry/curve/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-9.2/bin/nvcc -G -g -O0 -Xcompiler -dc -std=c++11 -gencode arch=compute_50,code=sm_50  -odir "src/geometry/curve" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-9.2/bin/nvcc  -G -g -O0 -Xcompiler -dc -std=c++11 --compile --relocatable-device-code=true -gencode arch=compute_50,code=compute_50 -gencode arch=compute_50,code=sm_50  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


