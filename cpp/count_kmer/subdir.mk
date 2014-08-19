
# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
FreqKmer.o \
classData.o \
classPattern.o \
count_kmer.o \
test.o \
testPart2.o \
test_intramacro.o

CPP_SRCS += \
FreqKmer.cpp \
classData.cpp \
classPattern.cpp \
count_kmer.cpp \
test.cpp \
testPart2.cpp \
test_intramacro.cpp

OBJS += \
FreqKmer.o \
classData.o \
classPattern.o \
count_kmer.o \
test.o \
testPart2.o \
test_intramacro.o

CPP_DEPS += \
FreqKmer.d \
classData.d \
classPattern.d \
count_kmer.d \
test.d \
testPart2.d \
test_intramacro.d


# Each subdirectory must supply rules for building sources it contributes
%.o: %.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++0x -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


