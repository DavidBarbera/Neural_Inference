//////////////////////////////////////////////////////
///  Python API via CFFI
//   Windows  ---> .dll
/////////////////////////////////////////////////////

#include <nni.h>

#define UNITS 128
#define IN 26

using namespace std::chrono;
using namespace nni;



#define EXPORT_API __declspec(dllexport)

extern "C" {

	EXPORT_API float c_cffi_dtw_wav_model(const char* wav1, const char* wav2, const char* model)
	{
		return dtw_wav(wav1, wav2, model, true, true, true, true);
	}

	EXPORT_API int nni_wav_file_metadata(const char* wav, int* nsamples, int* channels, int* samplerate, int* nbits, int* format) {

		return Wav_metadata(wav, nsamples, channels, samplerate, nbits, format);
	}

	EXPORT_API int nni_read_file(const char* wav, float* data, int* nsamples, int* channels, int* samplerate) {

		//return wav_to_samples(wav, data, nsamples[0], samplerate[0]);
		return extract_audio_data(wav, nsamples, channels, samplerate, data);
	}

	EXPORT_API int nni_signal_to_mfcc(float* signal, float* mfcc, int nsamp, int srate
																			, float twindow_size, float tstep_size
																			//, int mfcc_coeff
																			//, int mfcc_lofreq, int mfcc_hifreq
																			, int* nframes, int* nmfcc) {

		return signal_to_mfcc(signal, mfcc, nsamp, srate, twindow_size, tstep_size, nframes, nmfcc);
	}

	EXPORT_API void nni_process_signal(float* signal_in, float* signal_out, int nsamp_, int srate_, int* newsamp_, int* nsrate_, bool resample, bool norm) {
		preprocess_signal(signal_in, signal_out, nsamp_, srate_, newsamp_, nsrate_, resample, norm);
	}

	EXPORT_API int nni_wav_to_mfcc(const char* wav_f, const char* mfcc_f) {

		return wav_to_mfcc(wav_f, mfcc_f);
	}

	EXPORT_API int nni_mfcc_to_pp(float* pp_array, float* mfcc, int timesteps, const char* model) {

		return mfcc_to_pp_array(pp_array, mfcc, timesteps, model);
	}

	EXPORT_API int nni_samples_to_pp(float* pp_array, float* samples, int nsamp, int srate, const char* model){//, bool resample=true, bool norm=true) {

		return samples_to_pp_array(pp_array, samples, nsamp, srate, model, true, true);// , resample, norm);
	}

	EXPORT_API int nni_get_timesteps(int nsamp, int srate) {

		return get_timesteps(nsamp, srate);
	}

	EXPORT_API int nni_wav_to_pp(float* pp_array, char* wav, const char* model) {

		return wav_to_pp_array(pp_array, wav, model, true, true); //, bool resample = true, bool norm = true)
	}

	EXPORT_API float nni_dtw(float* pp1, int T1, float* pp2, int T2) {

		return dtw_pp(pp1, T1, pp2, T2);
	}

}//extern "C"
