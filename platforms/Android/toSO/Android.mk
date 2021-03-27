include $(CLEAR_VARS)

# override strip command to strip all symbols from output library; no need to ship with those..
# cmd-strip = $(TOOLCHAIN_PREFIX)strip $1 

LOCAL_ARM_MODE  := arm
LOCAL_PATH      := $(NDK_PROJECT_PATH)
LOCAL_MODULE    := libNeural_Inference
LOCAL_C_INCLUDES +=  C:\android-ndk-r13b\sources\cxx-stl\llvm-libc++\include
LOCAL_CFLAGS    := -Werror -std=c++11 -Werror=write-strings -fpermissive -fexceptions
LOCAL_SRC_FILES := $(LOCAL_PATH)/../android_source.cpp \
$(LOCAL_PATH)/../../../src/riffio/Resample.cpp \
$(LOCAL_PATH)/../../../src/mfcc/MFCC.cpp \

#LOCAL_LDFLAGS := -v
LOCAL_LDLIBS    := -llog

include $(BUILD_SHARED_LIBRARY)
