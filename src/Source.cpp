#include "nni.h"

#define UNITS 128
#define IN 26

using namespace std::chrono;
using namespace nni;


int main(int argc, char **argv) {
	//Running tests:
	//high_resolution_clock::time_point t1 = high_resolution_clock::now();
	//test_samples_to_mfcc(); 
	//test_from_wav_to_samples();
	//test_mfcc_to_pt();
	//test_dtw_on_pts();

	//test_parts();


	////samples to pp-------------------------------
	//// prepare data:
	////Load Samples
	//std::ifstream infile;
	//infile.open("tests/test-data/duck_samples.bin", std::ios::binary | std::ios::in);
	//infile.seekg(0, infile.end);
	//int infile_length = infile.tellg();
	//infile.seekg(0, infile.beg);
	//char *buffer = new char[infile_length];
	//infile.read(buffer, infile_length);
	//infile.close();

	////Inputs:
	//float *samples = (float*)buffer;
	//int nsamp = infile_length / sizeof(float);
	//int srate = 22050;
	//printf("samples: %i", nsamp);

	////Outputs:
	//Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pp;

	//pp = samples_to_pp((float*&)samples, nsamp, srate);

	//delete[] buffer;

	////Outputting the result into a text file:
	////Text precision:
	//Eigen::IOFormat precision(9);

	//std::ofstream file;
	//std::printf("\n pp: (%i,%i)\n", pp.rows(), pp.cols());
	//file.open("pp.txt", std::ios::out);
	//file << pp.format(precision) << "\n";
	//file.close();
	//
	//char fullwav1[256];
	//char fullwav2[256];
	//const char *wavdirname = "C:/Users/dbarbera/Documents/Repositories/Neural_Inference/";
	//const char *wav1 = "tests/test-data/Henry-1-Andrew.wav";
	//const char *wav2 = "tests/test-data/Henry-2-Andrew.wav";

	//sprintf_s(fullwav1, "%s%s", wavdirname, wav1);

	//float dtw_distance = dtw_wav(wav1, wav2);
	//printf("\nDTW distance: %.9f", dtw_distance);

	//printf("Calculating dtw distance...\n");
	//printf("\nDTW distance: %.9f", dtw_Andrew1_vs_Andrew2());

	////wav vs wav:
	//printf("Calculating dtw distance...\n");
	//printf("\nDTW distance - wav vs wav: %.9f\n", dtw_Andrew1_vs_Andrew2_model());

	////samples vs wav:
	//printf("Calculating dtw distance...\n");
	//printf("\nDTW distance - samples vs wav: %.9f\n", dtw_Andrew1_vs_Andrew2_model_samples_wav());

	////samples vs samples:
	//printf("Calculating dtw distance...\n");
	//printf("\nDTW distance - samples vs samples: %.9f\n", dtw_Andrew1_vs_Andrew2_model_samples_samples());


	//samples vs pp
	//dtw_AmbulanceWav_vs_AmbulancePP();


	//dtw_AmbulanceWav_vs_AmbulancePP_AmbulancePP();
	
	//pp vs pp
	//high_resolution_clock::time_point start;
	//high_resolution_clock::time_point end;

	//start=high_resolution_clock::now();
	//dtw_AmbulanceWav_vs_AmbulancePP_AmbulancePP();
	////samples vs (pp,pp) and returns the minimum
	////dtw_AmbulancePP_vs_AmbulancePP();
	//end=high_resolution_clock::now();

	//auto dur_us = duration<double, std::micro>(end - start).count();
	//auto dur_ms = duration<double, std::milli>(end - start).count();
	//printf("Time: %lfus %lfms", dur_us, dur_ms);

	//printf("Exe direcory:\n%s\n", ExePath().data());
	////printf(" wavs: \n%s\n", wav_path.data());

	//wavs_to_mfccs_files();
	//wavs_to_pp_files();


	//test_align();
	//test_data_storage();
	
	//high_resolution_clock::time_point start;
	//high_resolution_clock::time_point end;

	//start=high_resolution_clock::now();
	////test_pp_align_on_dictionary();
	//wavs_to_pp_files();
	//wavs_to_mfccs_files();
	//end= high_resolution_clock::now();

	//auto dur_us = duration<double, std::micro>(end - start).count();
	//auto dur_ms = duration<double, std::milli>(end - start).count();
	//printf("\n\nTime: %lfus %lfms\n", dur_us, dur_ms);
	

	//high_resolution_clock::time_point start;
	//high_resolution_clock::time_point end;

	//start = high_resolution_clock::now();
	//test_SLTfeatures_v1();
	//end = high_resolution_clock::now();

	//auto dur_us = duration<double, std::micro>(end - start).count();
	//auto dur_ms = duration<double, std::milli>(end - start).count();
	//printf("\n\nTime: %lfus %lfms\n", dur_us, dur_ms);


	//dtw_at_scale();

	//test_to_wav_file();

	//timestamp();
	//test_c_debug_dtw_samples_vs_pp();
    //dtw_wav_vs_pp_at_scale();
	//wavs_to_pp_files_resample_norm();

	////---- Exec self-contained programs -----
	//dtw_at_scale();
	//dtw_at_scale_reproducibility();
	//dtw_wav_vs_pp_at_scale_reproducibility();
	//dtw_wav_vs_pp_at_scale();
	//wavs_to_pp_files();
	//wavs_to_mfccs_files();
	////---------------------------------------
	int srate = 16000;
	float twindow_size = 0.03;
	float tstep_size = 0.01;
	int nsamp = 7339;

	int _window_size = (int)(0.5f + ((float)srate) * ((float)twindow_size));
	int _step_size = (int)(0.5f + ((float)srate) * ((float)tstep_size));
	printf("_window_size: %i\n_step_size: %i",_window_size, _step_size);
	int _nsamp = nsamp;
	int _srate = srate;
	//int _mfcc_lofreq = mfcc_lofreq;
	//int _mfcc_hifreq = mfcc_hifreq;
	//int _mfcc_coeff = mfcc_coeff;
	int _nframes = (int)((_nsamp - _window_size + _step_size) / _step_size);

	printf("\nwindow: %i\ninstep: %i\nsrate: %i\nnsamp: %i\n\nnframes: %i\n", _window_size, _step_size, srate, nsamp, _nframes);

	test_signal_to_mfcc();


	//--------------------------------------
	printf("\n\n Press a key to exit.\n");
	getchar();
	return 0;
}


//////////////////////////
////  UNITY API
////  Windows  ---> .dll
////////////////////////
//
//#define EXPORT_API __declspec(dllexport)
//
//extern "C" {
//
//	EXPORT_API float c_dtw_samples_samples_model(float* samples1, int nsamples1, int srate1, float* samples2, int nsamples2, int srate2, const char* model, bool resamp1=true, bool resamp2=true, bool norm1=true, bool norm2=true)
//	{
//		return dtw_samples_samples(samples1, nsamples1, srate1, samples2, nsamples2, srate2, model, resamp1, resamp2, norm1, norm2);
//	}
//
//	EXPORT_API float c_dtw_samples_wav_model(float* samples, int nsamp, int srate, const char* wav, const char* model, bool resamps=true, bool resampw=true, bool norms=true, bool normw = true)
//	{
//		return dtw_samples_wav(samples, nsamp, srate, wav, model, resamps, resampw, norms, normw);
//	}
//
//	EXPORT_API float c_dtw_wav_model(const char* wav1, const char* wav2, const char* model, bool resamp1=true, bool resamp2=true, bool norm1=true, bool norm2=true)
//	{
//		return dtw_wav(wav1, wav2, model, resamp1, resamp2, norm1, norm2);
//	}
//
//	EXPORT_API float c_cffi_dtw_wav_model(const char* wav1, const char* wav2, const char* model)
//	{
//		return dtw_wav(wav1, wav2, model, true, true, true, true);
//	}
//	
//	//EXPORT_API float c_dtw_wav(const char* wav1, const char* wav2)
//	//{
//	//	return dtw_wav(wav1, wav2);
//	//}
//
//	//EXPORT_API float c_dtw_Andrew1_vs_Andrew2()
//	//{
//	//	return dtw_Andrew1_vs_Andrew2();
//	//}
//
//	//EXPORT_API int c_addition(int a, int b)
//	//{
//	//	return a+b;
//	//}
//
//	//EXPORT_API int c_can_load_wav_file(const char *wav)
//	//{
//	//	return can_load_wav_file(wav);
//	//}
//
//	EXPORT_API float c_samples_vs_pp_e(float* samples, int nsamp, int srate, const char* ppfile, const char* model, int& err, bool resamp = true, bool norm = true) {
//		return c_dtw_samples_vs_pp(samples, nsamp, srate, ppfile, model, err, resamp, norm);
//	}
//
//	EXPORT_API float c_samples_vs_pp_pp_e(float* samples, int nsamp, int srate, const char* ppfile1, const char* ppfile2, const char* model, int& err, bool resamp = true, bool norm = true) {
//		return c_dtw_samples_vs_pp_pp(samples, nsamp, srate, ppfile1, ppfile2, model, err, resamp, norm);
//	}
//
//	EXPORT_API float c_pp_vs_pp_e(const char* ppfile1, const char* ppfile2, int& err) {
//		return c_dtw_pp_vs_pp(ppfile1, ppfile2, err);
//	}
//
//	EXPORT_API float c_debug_samples_vs_pp_e(float* samples, int nsamp, int srate, const char* ppfile, const char* model, const char* bindir, const char* word, int& err, bool resample = true, bool norm = true) {
//		return c_debug_dtw_samples_vs_pp(samples, nsamp, srate, ppfile, model, bindir, word, err, resample, norm);
//	}
//
//	EXPORT_API float c_debug_samples_vs_pp_store_samples_to_bin_only(float* samples, int nsamp, int srate, const char* ppfile, const char* model, const char* bindir, const char* word, const char* gender,	 int& err, bool resample = true, bool norm = true) {
//		return c_debug_dtw_samples_vs_pp_store_samples_to_bin_only(samples, nsamp, srate, ppfile, model, bindir, word, gender, err, resample, norm);
//	}
//
//}//extern "C"
//
//
//////////////////////////
////  UNITY API
////  ANDROID ---> .so
////////////////////////
//#ifdef __cplusplus
//extern "C" {
//
//	float c_dtw_samples_samples_model(float* samples1, int nsamples1, int srate1, float* samples2, int nsamples2, int srate2, const char* model, bool resamp1 = true, bool resamp2 = true, bool norm1 = true, bool norm2 = true)
//	{
//		return dtw_samples_samples(samples1, nsamples1, srate1, samples2, nsamples2, srate2, model, resamp1, resamp2, norm1, norm2);
//	}
//
//	float c_dtw_samples_wav_model(float* samples, int nsamp, int srate, const char* wav, const char* model, bool resamps = true, bool resampw = true, bool norms = true, bool normw = true)
//	{
//		return dtw_samples_wav(samples, nsamp, srate, wav, model, resamps, resampw, norms, normw);
//	}
//
//	float c_dtw_wav_model(const char* wav1, const char* wav2, const char* model, bool resamp1 = true, bool resamp2 = true, bool norm1 = true, bool norm2 = true)
//	{
//		return dtw_wav(wav1, wav2, model, resamp1, resamp2, norm1, norm2);
//	}
//
//	//float c_dtw_wav(const char* wav1, const char* wav2)
//	//{
//	//	return dtw_wav(wav1, wav2);
//	//}
//
//	//float c_dtw_Andrew1_vs_Andrew2()
//	//{
//	//	return dtw_Andrew1_vs_Andrew2();
//	//}
//
//	//int c_addition(int a, int b)
//	//{
//	//	return a + b;
//	//}
//
//	//int c_can_load_wav_file(const char* wav)
//	//{
//	//	return can_load_wav_file(wav);
//	//}
//
//	float c_samples_vs_pp_e(float* samples, int nsamp, int srate, const char* ppfile, const char* model, int& err, bool resamp = true, bool norm = true) {
//		return c_dtw_samples_vs_pp(samples, nsamp, srate, ppfile, model, err, resamp, norm);
//	}
//
//	float c_samples_vs_pp_pp_e(float* samples, int nsamp, int srate, const char* ppfile1, const char* ppfile2, const char* model, int& err, bool resamp = true, bool norm = true) {
//		return c_dtw_samples_vs_pp_pp(samples, nsamp, srate, ppfile1, ppfile2, model, err, resamp, norm);
//	}
//
//	float c_pp_vs_pp_e(const char* ppfile1, const char* ppfile2, int& err) {
//		return c_dtw_pp_vs_pp(ppfile1, ppfile2, err);
//	}
//
//	float c_debug_samples_vs_pp_e(float* samples, int nsamp, int srate, const char* ppfile, const char* model, const char* bindir, const char* word, int& err, bool resample = true, bool norm = true) {
//		return c_debug_dtw_samples_vs_pp(samples, nsamp, srate, ppfile, model, bindir, word, err, resample, norm);
//	}
//
//	float c_debug_samples_vs_pp_store_samples_to_bin_only(float* samples, int nsamp, int srate, const char* ppfile, const char* model, const char* bindir, const char* word, int& err, bool resample = true, bool norm = true) {
//		return c_debug_dtw_samples_vs_pp_store_samples_to_bin_only(samples, nsamp, srate, ppfile, model, bindir, word, err, resample, norm);
//	}
//
//
//}//extern "C"
//#endif //__cplusplus




//-----------------------------------------------------------------------------------------------------------------------------------------------
/////////////////////////////
///Python API via CFFI
//////////////////////////////

//
////////////////////////
//  Windows  ---> .dll
////////////////////////

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

/////////////////////////////
////  ANDROID/Linux ---> .so
/////////////////////////////
//#ifdef __cplusplus
//extern "C" {
//
//	
//float c_cffi_dtw_wav_model(const char* wav1, const char* wav2, const char* model)
//{
//	return dtw_wav(wav1, wav2, model, true, true, true, true);
//}
//}//extern "C"
//#endif //__cplusplus