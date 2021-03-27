#include <nni.h>

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

