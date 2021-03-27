include $(CLEAR_VARS)

# override strip command to strip all symbols from output library; no need to ship with those..
# cmd-strip = $(TOOLCHAIN_PREFIX)strip $1 

LOCAL_ARM_MODE  := arm
LOCAL_PATH      := $(NDK_PROJECT_PATH)
LOCAL_MODULE    := lib_Recogniser
#LOCAL_C_INCLUDES += C:\android-ndk-r10e\sources\cxx-stl\stlport\stlport
LOCAL_C_INCLUDES +=C:\android-ndk-r10e\sources\cxx-stl\gnu-libstdc++\4.9\include
LOCAL_CFLAGS    := -Werror -std=c++11 -Werror=write-strings -fpermissive -fexceptions
LOCAL_SRC_FILES := $(LOCAL_PATH)/../PhoneScoreProject.cpp \
$(LOCAL_PATH)/../common/UCLVisemeRec.cpp \
$(LOCAL_PATH)/../common/UCLPhoneRec.cpp \
$(LOCAL_PATH)/../common/RIFFIO.cpp \
$(LOCAL_PATH)/../common/Resample.cpp \
$(LOCAL_PATH)/../common/MFCC.cpp \
#$(LOCAL_PATH)/../common/LPCAnal.cpp \
#$(LOCAL_PATH)/../common/GetOpion.cpp \
#$(LOCAL_PATH)/../common/GenderClassify.cpp \
$(LOCAL_PATH)/../common/complex.cpp

LOCAL_LDLIBS    := -llog

include $(BUILD_SHARED_LIBRARY)