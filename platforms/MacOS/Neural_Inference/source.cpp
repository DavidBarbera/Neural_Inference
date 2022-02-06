//////////////////////////
////  UNITY API
////  ANDROID ---> .so\
////   MACOS ----> .dylib/bundle
////////////////////////


#include "../../../src/nni.h"

#define UNITS 128
#define IN 26

using namespace std::chrono;
using namespace nni;



#ifdef __cplusplus
extern "C" {

	float c_dtw_samples_samples_model(float* samples1, int nsamples1, int srate1, float* samples2, int nsamples2, int srate2, const char* model, bool resamp1 = true, bool resamp2 = true, bool norm1 = true, bool norm2 = true)
	{
		return dtw_samples_samples(samples1, nsamples1, srate1, samples2, nsamples2, srate2, model, resamp1, resamp2, norm1, norm2);
	}

	float c_dtw_samples_wav_model(float* samples, int nsamp, int srate, const char* wav, const char* model, bool resamps = true, bool resampw = true, bool norms = true, bool normw = true)
	{
		return dtw_samples_wav(samples, nsamp, srate, wav, model, resamps, resampw, norms, normw);
	}

	float c_dtw_wav_model(const char* wav1, const char* wav2, const char* model, bool resamp1 = true, bool resamp2 = true, bool norm1 = true, bool norm2 = true)
	{
		return dtw_wav(wav1, wav2, model, resamp1, resamp2, norm1, norm2);
	}

	//float c_dtw_wav(const char* wav1, const char* wav2)
	//{
	//	return dtw_wav(wav1, wav2);
	//}

	//float c_dtw_Andrew1_vs_Andrew2()
	//{
	//	return dtw_Andrew1_vs_Andrew2();
	//}

	//int c_addition(int a, int b)
	//{
	//	return a + b;
	//}

	//int c_can_load_wav_file(const char* wav)
	//{
	//	return can_load_wav_file(wav);
	//}

	float c_samples_vs_pp_e(float* samples, int nsamp, int srate, const char* ppfile, const char* model, int& err, bool resamp = true, bool norm = true) {
		return c_dtw_samples_vs_pp(samples, nsamp, srate, ppfile, model, err, resamp, norm);
	}

	float c_samples_vs_pp_pp_e(float* samples, int nsamp, int srate, const char* ppfile1, const char* ppfile2, const char* model, int& err, bool resamp = true, bool norm = true) {
		return c_dtw_samples_vs_pp_pp(samples, nsamp, srate, ppfile1, ppfile2, model, err, resamp, norm);
	}

	float c_pp_vs_pp_e(const char* ppfile1, const char* ppfile2, int& err) {
		return c_dtw_pp_vs_pp(ppfile1, ppfile2, err);
	}

	float c_debug_samples_vs_pp_e(float* samples, int nsamp, int srate, const char* ppfile, const char* model, const char* bindir, const char* word, int& err, bool resample = true, bool norm = true) {
		return c_debug_dtw_samples_vs_pp(samples, nsamp, srate, ppfile, model, bindir, word, err, resample, norm);
	}

	float c_debug_samples_vs_pp_store_samples_to_bin_only(float* samples, int nsamp, int srate, const char* ppfile, const char* model, const char* bindir, const char* word,const char* gender, int& err, bool resample = true, bool norm = true) {
		return c_debug_dtw_samples_vs_pp_store_samples_to_bin_only(samples, nsamp, srate, ppfile, model, bindir, word, gender, err, resample, norm);
	}


}//extern "C"
#endif //__cplusplus


