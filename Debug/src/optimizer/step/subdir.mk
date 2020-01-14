################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../src/optimizer/step/Newton.cu 

OBJS += \
./src/optimizer/step/Newton.o 

CU_DEPS += \
./src/optimizer/step/Newton.d 


# Each subdirectory must supply rules for building sources it contributes
src/optimizer/step/%.o: ../src/optimizer/step/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-9.2/bin/nvcc -G -g -O0 -Xcompiler -dc -std=c++11 -gencode arch=compute_50,code=sm_50  -odir "src/optimizer/step" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-9.2/bin/nvcc  -G -g -O0 -Xcompiler -dc -std=c++11 --compile --relocatable-device-code=true -gencode arch=compute_50,code=compute_50 -gencode arch=compute_50,code=sm_50  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


