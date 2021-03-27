from cffi import FFI
from ctypes.util import find_library as _find_library
import os as _os
import sys as _sys
import numpy as np






_ffi = FFI()

_ffi.cdef(r''' 
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
    int nni_mfcc_to_pp(float* pp_array, float* mfcc, int timesteps, const char* model);
    int nni_samples_to_pp(float* pp_array, float* samples, int nsamp, int srate, const char* model);
    int nni_get_timesteps(int nsamp, int srate);
    int nni_wav_to_pp(float* pp_array, char* wav, const char* model);
    float nni_dtw(float* pp1, int T1, float* pp2, int T2);
    ''')

#_nni = _ffi.dlopen("Neural_Inference.dll")

try:
    _libname = _find_library('Neural_Inference.dll')
    if _libname is None:
        raise OSError('Neural_Inference library not found')
    nni = _ffi.dlopen(_libname)
except OSError:
    if _sys.platform == 'darwin':
        _libname = 'libNeural_Inference.dylib'
    elif _sys.platform == 'win32':
        from platform import architecture as _architecture
        #_libname = 'Neural_Inference' + _architecture()[0] + '.dll'
        _libname = 'Neural_Inference' + '.dll'

    else:
        raise

    # hack for packaging tools like cx_Freeze, which
    # compress all scripts into a zip file
    # which causes __file__ to be inside this zip file

    _path = _os.path.dirname(_os.path.abspath(__file__))

    while not _os.path.isdir(_path):
        _path = _os.path.abspath(_os.path.join(_path, '..'))

    _nni = _ffi.dlopen(_os.path.join(
        _path, '', _libname))

# __libsndfile_version__ = _ffi.string(_snd.sf_version_string()).decode('utf-8', 'replace')
# if __libsndfile_version__.startswith('libsndfile-'):
#     __libsndfile_version__ = __libsndfile_version__[len('libsndfile-'):]

def dtw_wavs( wav1, wav2 ):

        # print(wav1)
        # print(wav2)
    wavfile1=_ffi.new("char[]", wav1.encode('utf-8'))
    wavfile2=_ffi.new("char[]", wav2.encode('utf-8'))
    
    _model=_os.path.join( _path, '', "model_weights.bin")

    model=_ffi.new("char[]", _model.encode('utf-8')) 
    
    return _nni.c_cffi_dtw_wav_model( wavfile1, wavfile2, model)


_ffi_types = {
    'float64': 'double',
    'float32': 'float',
    'int32': 'int',
    'int16': 'short'
}

def create_empty_array_mfcc(nframes, nmfcc=26, dtype='float32'):
    import numpy as np
    shape=(nframes,nmfcc)
    
    return np.empty(shape, dtype, order='C')

def create_empty_array(nsamples, dtype='float32'):
    """Create an empty array with appropriate shape."""
    import numpy as np
    shape = nsamples,
    
    return np.empty(shape, dtype, order='C')

def check_dtype(dtype):
    """Check if dtype string is valid and return ctype string."""
    try:
        return _ffi_types[dtype]
    except KeyError:
        raise ValueError("dtype must be one of {0!r} and not {1!r}".format(
            sorted(_ffi_types.keys()), dtype))

def mh_read_wav(wavfile):
    nsamp=_ffi.new("int *")
    srate=_ffi.new("int *")
    channels=_ffi.new("int *")
    nbits=_ffi.new("int *")
    wavformat=_ffi.new("int *")

    
    bwavfile=bytes(wavfile,'utf-8')
    cwavfile=_ffi.new("char[]", bwavfile )
    nsamples=_nni.nni_wav_file_metadata(cwavfile,nsamp,channels,srate,nbits,wavformat)
    nchan=int(channels[0])
    #print(f"result: {nsamples}\n samples: {int(nsamp[0])}\nchannels: {int(channels[0])}")
    #print(f"srate: {int(srate[0])}\n nbits: {int(nbits[0])}\nformat: {int(wavformat[0])}")
    data=create_empty_array(nsamples,'float32')
    """Check array and call low-level IO function."""
    if (data.ndim not in (1, 2) or
            data.ndim == 1 and nchan != 1 or
            data.ndim == 2 and data.shape[1] != nchan):
        raise ValueError("Invalid shape: {0!r}".format(data.shape))
    if not data.flags.c_contiguous:
        raise ValueError("Data must be C-contiguous")

    ctype = check_dtype(data.dtype.name)
    #print(ctype)

    assert data.dtype.itemsize == _ffi.sizeof(ctype)
    #print(data.dtype.itemsize)
    #print(ffi.sizeof(ctype))
    cdata = _ffi.cast(f"{ctype}*", data.__array_interface__['data'][0])
    _nni.nni_read_file(cwavfile, cdata, nsamp, channels, srate)
    
    return data, int(srate[0])


def embeddings(samples, sr=16000):
    nsamples=samples.shape[0]
    timesteps = _nni.nni_get_timesteps(nsamples, sr)

    #Samples
    csamples = _ffi.cast("float *", samples.ctypes.data)

    #pp
    pp = np.zeros((timesteps,45),dtype=np.float32)
    c_pp = _ffi.cast("float *", pp.ctypes.data)

    #Model
    _model=_os.path.join( _path, '', "model_weights.bin")
    model=_ffi.new("char[]", _model.encode('utf-8')) 

    nsamp=_ffi.new("int *")
    srate=_ffi.new("int *")

    #print(sr)
    nsamp[0]=nsamples
    srate[0]=sr
    #print(int(nsamp[0]),int(srate[0]))
    ns=_nni.nni_samples_to_pp(c_pp, csamples, int(nsamp[0]), int(srate[0]), model)#, resample, norm)
    #print(nsamp, srate)
    #print(int(nsamp[0]),int(srate[0]))

    #pp = np.frombuffer(_ffi.buffer(cpp, timesteps*45*4), dtype=np.float32).reshape((-1,45)) #float32 is 4 bytes

    return pp #, timesteps


def wav_to_pp(wavfile):
    #Model
    _model=_os.path.join( _path, '', "model_weights.bin")
    model=_ffi.new("char[]", _model.encode('utf-8'))

    nsamp=_ffi.new("int *")
    srate=_ffi.new("int *")
    channels=_ffi.new("int *")
    nbits=_ffi.new("int *")
    wavformat=_ffi.new("int *")

    
    bwavfile=bytes(wavfile,'utf-8')
    cwavfile=_ffi.new("char[]", bwavfile )
    nsamples=_nni.nni_wav_file_metadata(cwavfile,nsamp,channels,srate,nbits,wavformat)
    nchan=int(channels[0])

    timesteps = _nni.nni_get_timesteps(int(nsamp[0]), int(srate[0]))

    pp = np.zeros((timesteps,45),dtype=np.float32)
    c_pp = _ffi.cast("float *", pp.ctypes.data)

    _nni.nni_wav_to_pp(c_pp, cwavfile, model)

    return pp


def dtw(pp1, pp2):

    T1=pp1.shape[0]
    T2=pp2.shape[0]

    assert pp1.shape[1]==45
    assert pp2.shape[1]==45
    assert pp1.dtype==np.float32
    assert pp2.dtype==np.float32

    c_pp1 = _ffi.cast("float *", pp1.ctypes.data)
    c_pp2 = _ffi.cast("float *", pp2.ctypes.data)

    dtw_distance = _nni.nni_dtw(c_pp1, T1, c_pp2, T2)

    return dtw_distance