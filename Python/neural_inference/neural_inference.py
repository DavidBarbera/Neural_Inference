from cffi import FFI
import os

ffi = FFI()

ffi.cdef(r''' 
    float c_cffi_dtw_wav_model(const char* wav1, const char* wav2, const char* model);
    int nni_wav_file_metadata(const char* wav, int* nsamples, int* channels, int* samplerate, int* nbits, int* format);
    int nni_read_file(const char* wav, float* data, int* nsamples, int* channels, int* samplerate);
    int nni_signal_to_mfcc(float* signal, float* mfcc, int nsamp, int srate
                                                                            , float twindow_size, float tstep_size
                                                                            //, int mfcc_coeff
                                                                            //, int mfcc_lofreq, int mfcc_hifreq
                                                                            , int* nframes, int* nmfcc);
                                                                        
    void nni_process_signal(float* signal_in, float* signal_out, int nsamp_, int srate_, int* newsamp_, int* nsrate_, bool resample, bool norm);    
    int nni_wav_to_mfcc(const char* wav_f, const char* mfcc_f);
    ''')

_nni = ffi.dlopen("Neural_Inference.dll")


def dtw_wavs( wav1, wav2 ):
    
    wavfile1=ffi.new("char[]", wav1.encode('utf-8'))
    wavfile2=ffi.new("char[]", wav2.encode('utf-8'))
    
    model=ffi.new("char[]", b"model_weights.bin") 
    
    return _nni.c_cffi_dtw_wav_model( wavfile1, wavfile2, model)
    