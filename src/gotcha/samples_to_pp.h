
//#include <assert.h>

namespace nni {
	namespace gotcha {

		void saveData(std::string fileName, Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> matrix) {
			//https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
			const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");

			std::ofstream file(fileName);
			if (file.is_open())
			{
				file << matrix.format(CSVFormat);
				file.close();
			}
		}


		//template <typename Derived>
		Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> mfcc_to_pp(float*& rmfcc, int timesteps) {//float*& rmfcc //Eigen::MatrixBase<Derived> e_input
			//Dynamic parameters
			//This should come from a flatbuffer, at the moment is just hardcoded
			int units = 128;
			int in = 26;
			int out = 45;
			int T = timesteps; // this parameter comes on runtime
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> e_input(rmfcc, T, in);

			//std::ofstream outfile;
			//outfile.open("tests/test-data/duck_mfcc_INTEGRATION.bin", std::ios::binary | std::ios::out);
			//outfile.write((char*)rmfcc, (T - 1)*(in) * sizeof(float));
			//outfile.close();

			//Loading model: reading weigths
			std::ifstream infile;
			const char* model_weights_path = "C:/Users/dbarbera/Documents/Repositories/Neural_Inference/tests/test-data/model_weights.bin";
			infile.open(model_weights_path, std::ios::binary | std::ios::in);
			infile.seekg(0, infile.end);
			int infile_length = infile.tellg();
			infile.seekg(0, infile.beg);
			//bin_size = ((26 * 384)+(128*384)+(1*384)) * sizeof(float);//just forward
			//int bin_size = (2 * ((in * 3 * units) + (units * 3 * units) + (1 * 3 * units))) * sizeof(float); //forward and backward
			char* buffer_l1 = new char[infile_length];
			infile.read(buffer_l1, infile_length * sizeof(char));
			infile.close();
			float* weights = (float*)buffer_l1;
			//std::printf("bin_size: %i\n", bin_size);
			//std::printf("file_length: %i\n", infile_length);


			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l1(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l2(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l3(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l4(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l5(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l6(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l7(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l8(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l9(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l10(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l11(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l12(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l13(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l14(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l15(T, out);


			bidirectional_gru(output_l1, e_input, weights, units, in, T);
			batch_normalization(output_l2, output_l1, weights, 2 * units, T);
			bidirectional_gru(output_l3, output_l2, weights, units, 2 * units, T);
			batch_normalization(output_l4, output_l3, weights, 2 * units, T);
			bidirectional_gru(output_l5, output_l4, weights, units, 2 * units, T);
			batch_normalization(output_l6, output_l5, weights, 2 * units, T);
			bidirectional_gru(output_l7, output_l6, weights, units, 2 * units, T);
			batch_normalization(output_l8, output_l7, weights, 2 * units, T);
			bidirectional_gru(output_l9, output_l8, weights, units, 2 * units, T);
			batch_normalization(output_l10, output_l9, weights, 2 * units, T);
			bidirectional_gru(output_l11, output_l10, weights, units, 2 * units, T);
			batch_normalization(output_l12, output_l11, weights, 2 * units, T);
			bidirectional_gru(output_l13, output_l12, weights, units, 2 * units, T);
			batch_normalization(output_l14, output_l13, weights, 2 * units, T);
			timedistributed_dense(output_l15, output_l14, weights, 2 * units, out, T); //performs softmax as well


			delete[] buffer_l1;

			return output_l15;
		}

		Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> mfcc_to_pp(float*& rmfcc, int timesteps, const char* model_weights_path) {//float*& rmfcc //Eigen::MatrixBase<Derived> e_input
			//Dynamic parameters
			//This should come from a flatbuffer, at the moment is just hardcoded
			int units = 128;
			int in = 26;
			int out = 45;
			int T = timesteps; // this parameter comes on runtime
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> e_input(rmfcc, T, in);

			//std::ofstream outfile;
			//outfile.open("tests/test-data/duck_mfcc_INTEGRATION.bin", std::ios::binary | std::ios::out);
			//outfile.write((char*)rmfcc, (T - 1)*(in) * sizeof(float));
			//outfile.close();

			//Loading model: reading weigths
			std::ifstream infile;
			//const char *model_weights_path = "C:/Users/dbarbera/Documents/Repositories/Neural_Inference/tests/test-data/model_weights.bin";
			infile.open(model_weights_path, std::ios::binary | std::ios::in);
			infile.seekg(0, infile.end);
			int infile_length = infile.tellg();
			infile.seekg(0, infile.beg);
			//bin_size = ((26 * 384)+(128*384)+(1*384)) * sizeof(float);//just forward
			//int bin_size = (2 * ((in * 3 * units) + (units * 3 * units) + (1 * 3 * units))) * sizeof(float); //forward and backward
			char* buffer_l1 = new char[infile_length];
			infile.read(buffer_l1, infile_length * sizeof(char));
			infile.close();
			float* weights = (float*)buffer_l1;
			//std::printf("bin_size: %i\n", bin_size);
			//std::printf("file_length: %i\n", infile_length);


			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l1(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l2(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l3(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l4(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l5(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l6(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l7(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l8(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l9(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l10(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l11(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l12(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l13(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l14(T, 2 * units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_l15(T, out);


			bidirectional_gru(output_l1, e_input, weights, units, in, T);
			batch_normalization(output_l2, output_l1, weights, 2 * units, T);
			bidirectional_gru(output_l3, output_l2, weights, units, 2 * units, T);
			batch_normalization(output_l4, output_l3, weights, 2 * units, T);
			bidirectional_gru(output_l5, output_l4, weights, units, 2 * units, T);
			batch_normalization(output_l6, output_l5, weights, 2 * units, T);
			bidirectional_gru(output_l7, output_l6, weights, units, 2 * units, T);
			batch_normalization(output_l8, output_l7, weights, 2 * units, T);
			bidirectional_gru(output_l9, output_l8, weights, units, 2 * units, T);
			batch_normalization(output_l10, output_l9, weights, 2 * units, T);
			bidirectional_gru(output_l11, output_l10, weights, units, 2 * units, T);
			batch_normalization(output_l12, output_l11, weights, 2 * units, T);
			bidirectional_gru(output_l13, output_l12, weights, units, 2 * units, T);
			batch_normalization(output_l14, output_l13, weights, 2 * units, T);
			timedistributed_dense(output_l15, output_l14, weights, 2 * units, out, T); //performs softmax as well


			delete[] buffer_l1;

			//saveData(std::string("matrix.csv"), output_l15);

			return output_l15;
		}



		void preprocess_signal(float* signal_in, float* signal_out, int nsamp_, int srate_, int* newsamp_, int* nsrate_, bool resample = false, bool norm = false) {
			int nchan, orate, rnsamp;// , nsamp;

			//resampling (to keep compatibility with previous training methods where audio is 16PCM@16kHz and is converted to 16PCM@22050Hz
			orate = srate_;
			if (resample == true && orate != 22050) {
				CResample resamp;
				resamp.init(22050.0 / (srate_));
				int newsamp = (int)(22050.0 * (nsamp_) / (srate_));
				assert(newsamp == *newsamp_);
				int used;
				rnsamp = resamp.process(signal_in, nsamp_, 1, &used, signal_out, newsamp);
				resamp.close();

				*newsamp_ = rnsamp;
				*nsrate_ = 22050;
			}
			else {
				assert(nsamp_ == *newsamp_);
				std::memcpy(signal_out, signal_in, sizeof(float) * (nsamp_));
			}

			// normalise to unit variance
			if (norm == true) {
				double sumsq = 0;
				for (int i = 0; i < nsamp_; i++) sumsq += ((double)((signal_out)[i])) * ((double)((signal_out)[i]));
				double factor = std::sqrt(sumsq / (nsamp_));
				for (int i = 0; i < nsamp_; i++) ((signal_out)[i]) = (float)(((signal_out)[i]) / factor);
			}

		}//preprocess_signal

		void preprocess_samples(float*& rfsp, int& nsamp_, int& srate_, bool resample = true, bool norm = true) {
			int nchan, orate, rnsamp;// , nsamp;
			float* fsp;

			//resampling (to keep compatibility with previous training methods where audio is 16PCM@16kHz and is converted to 16PCM@22050Hz
			orate = srate_;
			if (resample == true && orate != 22050) {
				fsp = (float*)calloc(nsamp_, sizeof(float));
				//assert(fsp == 0);
				std::memcpy(fsp, rfsp, sizeof(float) * nsamp_);

				CResample resamp;
				resamp.init(22050.0 / srate_);
				int newsamp = (int)(22050.0 * nsamp_ / srate_);
				rfsp = (float*)calloc(newsamp, sizeof(float));
				int used;
				rnsamp = resamp.process(fsp, nsamp_, 1, &used, rfsp, newsamp);
				free(fsp);
				resamp.close();

				nsamp_ = rnsamp;
				srate_ = 22050;

				//printf("\nresampling to 22050..");
			}

			// normalise to unit variance
			if (norm == true) {
				double sumsq = 0;
				for (int i = 0; i < nsamp_; i++) sumsq += ((double)((rfsp)[i])) * ((double)((rfsp)[i]));
				double factor = std::sqrt(sumsq / nsamp_);
				for (int i = 0; i < nsamp_; i++) ((rfsp)[i]) = (float)(((rfsp)[i]) / factor);
			}
		}

		inline Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> samples_to_pp(float*& samples, int& nsamp, int& srate, bool resample = true, bool norm = true) {

			float* rmfcc = NULL;
			int nframes = 0;
			int nmfcc = 0;
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pp;

			preprocess_samples(samples, nsamp, srate, resample, norm);

			MFCC c;

			int timesteps = c.get_mfcc((float*&)samples, nsamp, srate, (float*&)rmfcc, &nframes, &nmfcc);
			int in = 26;
			//Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mfcc_(rmfcc,timesteps,in);
			pp = mfcc_to_pp((float*&)rmfcc, timesteps);
			//pp = mfcc_to_pp(mfcc_, timesteps);

			return pp;
		}

		inline Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> samples_to_pp(float*& samples, int& nsamp, int& srate, const char* model) {

			float* rmfcc = NULL;
			int nframes = 0;
			int nmfcc = 0;
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pp;


			MFCC c;

			int timesteps = c.get_mfcc((float*&)samples, nsamp, srate, (float*&)rmfcc, &nframes, &nmfcc);
			int in = 26;
			//Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mfcc_(rmfcc,timesteps,in);
			//pp = mfcc_to_pp((float*&)rmfcc, timesteps);
			pp = mfcc_to_pp((float*&)rmfcc, timesteps, model);
			//pp = mfcc_to_pp(mfcc_, timesteps);

			return pp;
		}

		inline Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> samples_to_pp(float*& samples, int& nsamp, int& srate, const char* model, bool resample = true, bool norm = true) {

			float* rmfcc = NULL;
			int nframes = 0;
			int nmfcc = 0;
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pp;

			preprocess_samples(samples, nsamp, srate, resample, norm);

			MFCC c;

			int timesteps = c.get_mfcc((float*&)samples, nsamp, srate, (float*&)rmfcc, &nframes, &nmfcc);
			int in = 26;
			//Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mfcc_(rmfcc,timesteps,in);
			//pp = mfcc_to_pp((float*&)rmfcc, timesteps);
			pp = mfcc_to_pp((float*&)rmfcc, timesteps, model);
			//pp = mfcc_to_pp(mfcc_, timesteps);

			return pp;
		}

		int get_timesteps(int nsamp, int srate) {

			int		winsize = (int)(0.5 + srate * MFCC_WINSIZE);
			int		stpsize = (int)(0.5 + srate * MFCC_STEPSIZE);
			CMFCC	mfcc(MFCC_COEFF, 1, winsize, MFCC_LOFREQ, MFCC_HIFREQ, srate);
			int		nframe = ((nsamp - winsize + stpsize) / stpsize);

			return nframe;
		}

		int samples_to_pp_array(float*& pp_array, float*& samples, int nsamp, int srate, const char* model, bool resample = true, bool norm = true) {

			//float* rmfcc = NULL;
			//int nframes = 0;
			//int nmfcc = 0;
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pp;

			//preprocess_samples(samples, nsamp, srate, resample, norm);

			//MFCC c;

			//int timesteps = c.get_mfcc((float*&)samples, nsamp, srate, (float*&)rmfcc, &nframes, &nmfcc);
			//int in = 26;
			////Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mfcc_(rmfcc,timesteps,in);
			////pp = mfcc_to_pp((float*&)rmfcc, timesteps);
			//pp = mfcc_to_pp((float*&)rmfcc, timesteps, model);
			////pp = mfcc_to_pp(mfcc_, timesteps);

			//Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> signal((float*& )samples, nsamp, 1);
			//saveData(std::string("signal_samples_to_pp.csv"), signal);

			pp = samples_to_pp((float*&)samples, nsamp, srate, model, resample, norm);

			std::copy(pp.data(), pp.data() + pp.size(), pp_array);

			return nsamp;
		}

		int mfcc_to_pp_array(float*& pp_array, float*& rmfcc, int timesteps, const char* model, bool resample = true, bool norm = true) {

			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pp;

			pp = mfcc_to_pp((float*&)rmfcc, timesteps, model);
			//pp = mfcc_to_pp(mfcc_, timesteps);

			std::copy(pp.data(), pp.data() + pp.size(), pp_array);

			return timesteps;
		}
		

		int wav_to_samples(const char* wav, float*& rfsp, int& nsamp_, int& srate_) {
			//float *samples = NULL;
			//int nsamp, srate;

			//wavfile_to_samples("tests/test-data/duck.wav", samples, nsamp, srate);
			int	nsamp, srate;
			//float *rfsp;

			int nchan, orate, rnsamp;
			float* fsp, * fsp2;

			//const char wavfile[] = "tests/test-data/duck.wav";
			if (LoadWav(wav, &nsamp_, &nchan, &srate_, &fsp, &fsp2) == 0) {
				fprintf(stderr, "failed to read WAV ond test \"test_from_wav_to_samples\"");
				//getchar();
				return -1;
			}
			if (nchan == 2) free(fsp2);
			rfsp = fsp;

			return 0;
		}

		//int wav_to_samples(const char* wav, float* &rfsp, int &nsamp_, int &srate_, bool resample=true, bool norm=true) {
		//	//float *samples = NULL;
		//	//int nsamp, srate;

		//	//wavfile_to_samples("tests/test-data/duck.wav", samples, nsamp, srate);
		//	int	nsamp, srate;
		//	//float *rfsp;

		//	int nchan, orate, rnsamp;
		//	float	*fsp, *fsp2;

		//	//const char wavfile[] = "tests/test-data/duck.wav";
		//	if (LoadWav(wav, &nsamp_, &nchan, &srate_, &fsp, &fsp2) == 0) {
		//		fprintf(stderr, "failed to read WAV ond test \"test_from_wav_to_samples\"");
		//		//getchar();
		//		return -1;
		//	}
		//	if (nchan == 2) free(fsp2);
		//	rfsp = fsp;
		//	preprocess_samples(rfsp, nsamp_, srate_, resample, norm);

		//	////resampling:
		//	//orate = srate;
		//	//if (resample==true && orate != 22050) {
		//	//	CResample resamp;
		//	//	resamp.init(22050.0 / srate);
		//	//	int newsamp = (int)(22050.0*nsamp / srate);
		//	//	rfsp = (float *)calloc(newsamp, sizeof(float));
		//	//	int used;
		//	//	nsamp = rnsamp = resamp.process(fsp, nsamp, 1, &used, rfsp, newsamp);
		//	//	free(fsp);
		//	//	resamp.close();

		//	//	nsamp_ = nsamp;
		//	//	srate_ = 22050;

		//	//	printf("\nresampling to 22050..");
		//	//} else {
		//	//	rfsp = fsp;
		//	//	rnsamp = nsamp;
		//	//	nsamp_ = nsamp;
		//	//	srate_ = srate;
		//	//}

		//	// normalise to unit
		//	//if (norm == true) {
		//	//	double sumsq = 0;
		//	//	double component;
		//	//	for (int i = 0; i < nsamp; i++) {
		//	//		component = (rfsp)[i];
		//	//		sumsq += component*component;
		//	//	}
		//	//	double factor = std::sqrt(sumsq / nsamp); //android
		//	//	for (int i = 0; i < nsamp; i++) {
		//	//		(rfsp)[i] = (float)((rfsp)[i] / factor);
		//	//	}
		//	//}


		//	//// normalise to unit variance
		//	//if (norm == true) {
		//	//	double sumsq = 0;
		//	//	for (int i = 0; i < nsamp; i++) sumsq += (rfsp)[i] * (rfsp)[i];
		//	//	double factor = std::sqrt(sumsq / nsamp);
		//	//	for (int i = 0; i < nsamp; i++) (rfsp)[i] = (float)((rfsp)[i] / factor);
		//	//}

		//	//for testing
		//	//std::ofstream outfile;
		//	//outfile.open("tests/test-data/duck_samples_test_out.bin", std::ios::binary | std::ios::out);
		//	//outfile.write((const char*)rfsp, sizeof(float)*(nsamp));
		//	//outfile.close();
		//	//printf(" *rnsamp: %i \n newsamp: %i", rnsamp, newsamp);
		//	//free(rfsp);
		//	return 0;
		//}



		Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> wav_to_pp(const char* wav, bool resample=true) {
			float *samples = NULL;
			int nsamp = 0, srate = 0;
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pp;

			wav_to_samples(wav, (float*&)samples, nsamp, srate);
			pp = samples_to_pp((float*&)samples, nsamp, srate);


			return pp;
		}

		Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> wav_to_pp(const char* wav, const char* model, bool resample=true, bool norm=true) {
			float *samples = NULL;
			int nsamp = 0, srate = 0;
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pp;

			wav_to_samples(wav, (float*&)samples, nsamp, srate);
			pp = samples_to_pp((float*&)samples, nsamp, srate, model, resample, norm);


			return pp;
		}

		int wav_to_pp_array(float*& pp_array, const char* wav, const char* model, bool resample = true, bool norm = true) {
			float* samples = NULL;
			int nsamp = 0, srate = 0;
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pp;

			wav_to_samples(wav, (float*&)samples, nsamp, srate);
			//Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> signal((float*&)samples, nsamp, 1);
			//saveData(std::string("signal_wav_to_pp.csv"), signal);

			pp = samples_to_pp((float*&)samples, nsamp, srate, model, resample, norm);

			//saveData(std::string("matrix_wav_to_pp.csv"), pp);
			std::copy(pp.data(), pp.data() + pp.size(), pp_array);

			return nsamp;
		}

		float dtw_wav(const char* wav1, const char* wav2, bool resample1 = true, bool resample2 = true) {
			std::ifstream infile;
			int infile_length = 0;
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pp1;
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pp2;

			pp1 = wav_to_pp(wav1,resample1);
			pp2 = wav_to_pp(wav2,resample2);

			float dtw_distance = dtw(pp1, pp2);


			//float dtw_distance = 0.99999f;

			return dtw_distance;		
		}

		float dtw_wav(const char* wav1, const char* wav2, const char* model, bool resample1=true, bool resample2=true, bool norm1=true, bool norm2=true) {
			//std::ifstream infile;
			//int infile_length = 0;
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pp1;
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pp2;

			pp1 = wav_to_pp(wav1, model, resample1, norm1);
			pp2 = wav_to_pp(wav2, model, resample2, norm2);

			float dtw_distance = dtw(pp1, pp2);


			//float dtw_distance = 0.99999f;

			return dtw_distance;
		}

		float dtw_pp(float* pp1_array, int T1, float* pp2_array, int T2) {

			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> pp1(pp1_array, T1, 45);
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> pp2(pp2_array, T2, 45);

			float dtw_distance = dtw(pp1, pp2);

			return dtw_distance;
		}

		float dtw_samples_wav(float*& samples, int& nsamp, int& srate, const char* wav, const char* model, bool resamps=true, bool resampw=true, bool norms=true, bool normw=true) {
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pp1;
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pp2;

			pp1 = samples_to_pp(samples, nsamp, srate, model, resamps, norms);
			pp2 = wav_to_pp(wav, model, resampw, normw);

			float dtw_distance = dtw(pp1, pp2);

			//float dtw_distance = 0.99999f;

			return dtw_distance;
		}


		float dtw_samples_samples(float*& samples1, int& nsamp1, int& srate1, float*& samples2, int& nsamp2, int& srate2, const char* model, bool resamp1=true, bool resamp2=true, bool norm1=true, bool norm2=true) {
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pp1;
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pp2;

			pp1 = samples_to_pp(samples1, nsamp1, srate1, model, resamp1, norm1);
			pp2 = samples_to_pp(samples2, nsamp2, srate2, model, resamp2, norm2);

			float dtw_distance = dtw(pp1, pp2);


			//float dtw_distance = 0.99999f;

			return dtw_distance;
		}

		float c_dtw_samples_vs_pp(float*& samples, int& nsamp, int& srate, const char* ppfile, const char* model, bool resamp=true, bool norm = true) {
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pp1;
			//Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pp2;
			
			pp1 = samples_to_pp(samples, nsamp, srate, model, resamp, norm);

			//Loading posteriors from binary file
			std::ifstream infile;
			infile.open(ppfile, std::ios::binary | std::ios::in);
			infile.seekg(0, infile.end);
			int infile_length = infile.tellg();
			infile.seekg(0, infile.beg);
			char *buffer = new char[infile_length];
			infile.read(buffer, infile_length * sizeof(char));
			infile.close();
			float *posteriors = (float*)buffer;

			int out = 45;
			int T = (int)infile_length/(sizeof(float)*out);

			printf("pp.shape: (%i,%i)\n", T, out);

			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> pp2(posteriors, T, out);

			float dtw_distance = dtw(pp1, pp2);
			printf("\nDTW distance: %.9f", dtw_distance);

			delete[] buffer;

			return dtw_distance;
		}

		float c_dtw_wav_vs_pp(const char* wav, const char* ppfile, const char* model, bool resample=true, bool norm = true) {
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pp1;
			//Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pp2;

			pp1 = wav_to_pp(wav, model, resample, norm);

			//Loading posteriors from binary file
			std::ifstream infile;
			infile.open(ppfile, std::ios::binary | std::ios::in);
			infile.seekg(0, infile.end);
			int infile_length = infile.tellg();
			infile.seekg(0, infile.beg);
			char* buffer = new char[infile_length];
			infile.read(buffer, infile_length * sizeof(char));
			infile.close();
			float* posteriors = (float*)buffer;

			int out = 45;
			int T = (int)infile_length / (sizeof(float) * out);

			//printf("pp.shape: (%i,%i)\n", T, out);

			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> pp2(posteriors, T, out);

			float dtw_distance = dtw(pp1, pp2);
			//printf("\nDTW distance: %.9f", dtw_distance);

			delete[] buffer;

			return dtw_distance;
		}
		
		void save_samples(float*& samples, int& nsamp, const char* filename) {
			std::ofstream outfile;
			//const char *model_weights_path = "C:/Users/dbarbera/Documents/Repositories/Neural_Inference/tests/test-data/model_weights.bin";
			outfile.open(filename, std::ios::binary | std::ios::out);
			//char* buffer = new char[nsamp*sizeof(float)];
			char* buffer = (char*)samples;
			outfile.write(buffer, nsamp * sizeof(float));
			outfile.close();
	
		}

		int to_binary_file(const char* filename, float*& samples, int& nsamp) {
			std::ofstream outfile;
			//printf("filename: %s", filename);
			outfile.open(filename, std::ios::binary | std::ios::out);
			if (!outfile.is_open()) {
				return 1;
			}
			//char* buffer = new char[nsamp * sizeof(float)];
			//buffer = (char*)samples;
			char* buffer = (char*)samples;
			outfile.write(buffer, nsamp * sizeof(float));
			outfile.close();
			//delete[] buffer;
			return 0;
		}

		void to_wav_file(const char* filename, float*& samples, int& nsamp, int& srate) {
			SaveWav(filename, nsamp, 1, srate, samples, NULL);
		}

		std::string timestamp() {
			std::tm time;// = std::chrono::system_clock::now();
			std::time_t nowt = std::time(0);// = std::chrono::system_clock::to_time_t(time);

			std::chrono::high_resolution_clock::time_point p = std::chrono::high_resolution_clock::now();

			std::chrono::milliseconds ms = std::chrono::duration_cast<std::chrono::milliseconds>(p.time_since_epoch());

			std::chrono::seconds ss = std::chrono::duration_cast<std::chrono::
				seconds>(ms);
			std::time_t t = ss.count();
			std::size_t fractional_seconds = ms.count() % 1000;

//This is platform dependent to make it thread-safe
#if defined(__unix__)
			localtime_r(&nowt, &time);
#elif defined(_MSC_VER)
			localtime_s(&time, &nowt);
#else
			static std::mutex mtx;
			std::lock_guard<std::mutex> lock(mtx);
			time = *std::localtime(&nowt);
#endif

			std::string s(24,'\0');// , ' ');
			std::strftime(&s[0], s.size(), "%Y-%m-%d-%H-%M-%S", &time);
			//printf("%s", s.data());
			//std::string name(s.data());
			//s += ".wav";
			std::string filename(24, '\0');
			//Only for windows:
			//sprintf_s(&filename[0], filename.size(),"%s-%3d", s.data(), fractional_seconds);
			std::snprintf(&filename[0], filename.size(), "%s-%3zu", s.data(), fractional_seconds);

			//printf("\n %s.%i \n", s.data(),fractional_seconds);
			//std::cout << filename << std::endl;
			//std::cout << time.tm_year+1900;

			return filename;
		}

		//void test_to_wav_file() {
		//	std::string wavfilein("tests/test-data/ambulance_F.wav");
		//	std::string wavfileout("tests/test-data/ambulance_F_out.wav");
		//	float* samples = NULL;
		//	int nsamp = 0, srate = 0;
		//	wav_to_samples(wavfilein.data(), (float*&)samples, nsamp, srate, false);

		//	to_wav_file(wavfileout.data(), samples, nsamp, srate);

		//}

		//version with error handling
		float c_dtw_samples_vs_pp(float*& samples, int& nsamp, int& srate, const char* ppfile, const char* model, int& err, bool resamp = true, bool norm = true) {
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pp1;

			pp1 = samples_to_pp(samples, nsamp, srate, model, resamp, norm);

			//Loading posteriors from binary file
			std::ifstream infile;
			infile.open(ppfile, std::ios::binary | std::ios::in);
			err = 0;
			if (!infile.is_open()) {
				err = 1;
				return 1111.f;
			}

			infile.seekg(0, infile.end);
			int infile_length = infile.tellg();
			infile.seekg(0, infile.beg);
			char* buffer = new char[infile_length];
			infile.read(buffer, infile_length * sizeof(char));
			infile.close();
			float* posteriors = (float*)buffer;

			int out = 45;
			int T = (int)infile_length / (sizeof(float) * out);

			printf("pp.shape: (%i,%i)\n", T, out);

			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> pp2(posteriors, T, out);

			float dtw_distance = dtw(pp1, pp2);
			printf("\nDTW distance: %.9f", dtw_distance);

			delete[] buffer;

			return dtw_distance;
		}


		float c_dtw_samples_vs_pp_pp(float*& samples, int& nsamp, int& srate, const char* ppfile1, const char* ppfile2, const char* model, bool resamp = true, bool norm = true) {
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pps;

			pps = samples_to_pp(samples, nsamp, srate, model, resamp, norm);

			//Posteriors on phonetis system of 45 phones (ARPAbet+sil)
			int out = 45;
			//Loading posteriors from binary file1
			std::ifstream infile1;
			infile1.open(ppfile1, std::ios::binary | std::ios::in);
			infile1.seekg(0, infile1.end);
			int infile1_length = infile1.tellg();
			infile1.seekg(0, infile1.beg);
			char *buffer1 = new char[infile1_length];
			infile1.read(buffer1, infile1_length * sizeof(char));
			infile1.close();
			float *posteriors1 = (float*)buffer1;

			
			int T1 = (int)infile1_length / (sizeof(float)*out);

			printf("pp.shape: (%i,%i)\n", T1, out);

			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> pp1(posteriors1, T1, out);


			//Loading posteriors from binary file2
			std::ifstream infile2;
			infile2.open(ppfile2, std::ios::binary | std::ios::in);
			infile2.seekg(0, infile2.end);
			int infile2_length = infile2.tellg();
			infile2.seekg(0, infile2.beg);
			char *buffer2 = new char[infile2_length];
			infile2.read(buffer2, infile2_length * sizeof(char));
			infile2.close();
			float *posteriors2 = (float*)buffer2;

			int T2 = (int)infile2_length / (sizeof(float)*out);

			printf("pp.shape: (%i,%i)\n", T2, out);

			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> pp2(posteriors2, T2, out);

			float dtw_distance1 = dtw(pps, pp1);
			float dtw_distance2 = dtw(pps, pp2);
			printf("\nDTW distance1: %.9f\nDTW distance2: %.9f\n", dtw_distance1, dtw_distance2);
			float dtw_distance = fmin(dtw_distance1, dtw_distance2);
			printf("\n min DTW distance: %.9f\n", dtw_distance);

			delete[] buffer2;
			delete[] buffer1;

			return dtw_distance;
		}

		//Version with minimal error handling
		float c_dtw_samples_vs_pp_pp(float*& samples, int& nsamp, int& srate, const char* ppfile1, const char* ppfile2, const char* model, int& err, bool resamp = true, bool norm = true) {
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pps;

			pps = samples_to_pp(samples, nsamp, srate, model, resamp, norm);
			err = 0;
			//Posteriors on phonetis system of 45 phones (ARPAbet+sil)
			int out = 45;
			//Loading posteriors from binary file1
			
			std::ifstream infile1;

			infile1.open(ppfile1, std::ios::binary | std::ios::in);

			if (!infile1.is_open()) {
				err = 1;
				return 1111.f;
			}

			infile1.seekg(0, infile1.end);
			int infile1_length = infile1.tellg();
			infile1.seekg(0, infile1.beg);
			char* buffer1 = new char[infile1_length];
			infile1.read(buffer1, infile1_length * sizeof(char));
			infile1.close();
			float* posteriors1 = (float*)buffer1;


			int T1 = (int)infile1_length / (sizeof(float) * out);

			printf("pp.shape: (%i,%i)\n", T1, out);

			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> pp1(posteriors1, T1, out);


			//Loading posteriors from binary file2
			std::ifstream infile2;

			infile2.open(ppfile2, std::ios::binary | std::ios::in);

			if (!infile2.is_open()) {
				err = 1;
				return 1111.f;
			}

			infile2.seekg(0, infile2.end);
			int infile2_length = infile2.tellg();
			infile2.seekg(0, infile2.beg);
			char* buffer2 = new char[infile2_length];
			infile2.read(buffer2, infile2_length * sizeof(char));
			infile2.close();
			float* posteriors2 = (float*)buffer2;

			int T2 = (int)infile2_length / (sizeof(float) * out);

			printf("pp.shape: (%i,%i)\n", T2, out);

			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> pp2(posteriors2, T2, out);

			float dtw_distance1 = dtw(pps, pp1);
			float dtw_distance2 = dtw(pps, pp2);
			printf("\nDTW distance1: %.9f\nDTW distance2: %.9f\n", dtw_distance1, dtw_distance2);
			float dtw_distance = fmin(dtw_distance1, dtw_distance2);
			printf("\n min DTW distance: %.9f\n", dtw_distance);

			delete[] buffer2;
			delete[] buffer1;

			return dtw_distance;
		}

		//------------- debugging the plugin ---------------
		float c_debug_dtw_samples_vs_pp(float*& samples, int& nsamp, int& srate, const char* ppfile, const char* model, const char* bindir, const char* word, int& err, bool resample = true, bool norm = true) {
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pp1;
			std::string dir(bindir);
			std::string sword(word);
			std::string before_norm(dir + "07-cc_plugin-" + timestamp().data() + "-" + sword);
			to_binary_file((before_norm+".bin").data(), samples, nsamp);

			pp1 = samples_to_pp(samples, nsamp, srate, model, resample, norm);

			std::ifstream infile;
			infile.open(ppfile, std::ios::binary | std::ios::in);
			err = 0;
			if (!infile.is_open()) {
				err = 1;
				return 1111.f;
			}

			infile.seekg(0, infile.end);
			int infile_length = infile.tellg();
			infile.seekg(0, infile.beg);
			char* buffer = new char[infile_length];
			infile.read(buffer, infile_length * sizeof(char));
			infile.close();
			float* posteriors = (float*)buffer;

			int out = 45;
			int T = (int)infile_length / (sizeof(float) * out);

			printf("pp.shape: (%i,%i)\n", T, out);

			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> pp2(posteriors, T, out);

			float dtw_distance = dtw(pp1, pp2);
			printf("\nDTW distance: %.9f", dtw_distance);

			delete[] buffer;

			return dtw_distance;
		}

		float c_debug_dtw_samples_vs_pp_store_samples_to_bin_only(float*& samples, int& nsamp, int& srate, const char* ppfile, const char* model, const char* bindir, const char* word, const char* gender, int& err, bool resample = true, bool norm = true) {
			
			std::string dir(bindir);
			std::string sword(word);
			std::string sgender(gender);
			std::string before_norm(dir + "05-cc_plugin-"+timestamp().data() + "-" + sword +"-"+sgender);

			err = 0;
			err=to_binary_file((before_norm + ".bin").data(), samples, nsamp);
			if (err != 0) {
				return 11111.f;
			}

			return 22222.f;
		}

		float c_debug_dtw_samples_vs_pp_pp(float*& samples, int& nsamp, int& srate, const char* ppfile1, const char* ppfile2, const char* model, bool norm = true) {
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pps;

			pps = samples_to_pp(samples, nsamp, srate, model, true, norm);

			//Posteriors on phonetis system of 45 phones (ARPAbet+sil)
			int out = 45;
			//Loading posteriors from binary file1
			std::ifstream infile1;
			infile1.open(ppfile1, std::ios::binary | std::ios::in);
			infile1.seekg(0, infile1.end);
			int infile1_length = infile1.tellg();
			infile1.seekg(0, infile1.beg);
			char* buffer1 = new char[infile1_length];
			infile1.read(buffer1, infile1_length * sizeof(char));
			infile1.close();
			float* posteriors1 = (float*)buffer1;


			int T1 = (int)infile1_length / (sizeof(float) * out);

			printf("pp.shape: (%i,%i)\n", T1, out);

			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> pp1(posteriors1, T1, out);


			//Loading posteriors from binary file2
			std::ifstream infile2;
			infile2.open(ppfile2, std::ios::binary | std::ios::in);
			infile2.seekg(0, infile2.end);
			int infile2_length = infile2.tellg();
			infile2.seekg(0, infile2.beg);
			char* buffer2 = new char[infile2_length];
			infile2.read(buffer2, infile2_length * sizeof(char));
			infile2.close();
			float* posteriors2 = (float*)buffer2;

			int T2 = (int)infile2_length / (sizeof(float) * out);

			printf("pp.shape: (%i,%i)\n", T2, out);

			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> pp2(posteriors2, T2, out);

			float dtw_distance1 = dtw(pps, pp1);
			float dtw_distance2 = dtw(pps, pp2);
			printf("\nDTW distance1: %.9f\nDTW distance2: %.9f\n", dtw_distance1, dtw_distance2);
			float dtw_distance = fmin(dtw_distance1, dtw_distance2);
			printf("\n min DTW distance: %.9f\n", dtw_distance);

			delete[] buffer2;
			delete[] buffer1;

			return dtw_distance;
		}

		float c_dtw_pp_vs_pp(const char* ppfile1, const char* ppfile2, int& err) {
			//Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pp1;
			//Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pp2;

			//Posteriors on phonetis system of 45 phones (ARPAbet+sil)
			int out = 45;
			//Loading posteriors from binary file1
			std::ifstream infile1;
			infile1.open(ppfile1, std::ios::binary | std::ios::in);
			if (!infile1.is_open()) {
				err = 1;
				return 1111.f;
			}
			infile1.seekg(0, infile1.end);
			int infile1_length = infile1.tellg();
			infile1.seekg(0, infile1.beg);
			char *buffer1 = new char[infile1_length];
			infile1.read(buffer1, infile1_length * sizeof(char));
			infile1.close();
			float *posteriors1 = (float*)buffer1;


			int T1 = (int)infile1_length / (sizeof(float)*out);

			printf("pp.shape: (%i,%i)\n", T1, out);

			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> pp1(posteriors1, T1, out);


			//Loading posteriors from binary file2
			std::ifstream infile2;
			infile2.open(ppfile2, std::ios::binary | std::ios::in);
			if (!infile2.is_open()) {
				err = 1;
				return 1111.f;
			}
			infile2.seekg(0, infile2.end);
			int infile2_length = infile2.tellg();
			infile2.seekg(0, infile2.beg);
			char *buffer2 = new char[infile2_length];
			infile2.read(buffer2, infile2_length * sizeof(char));
			infile2.close();
			float *posteriors2 = (float*)buffer2;

			int T2 = (int)infile2_length / (sizeof(float)*out);

			printf("pp.shape: (%i,%i)\n", T2, out);

			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> pp2(posteriors2, T2, out);

			float dtw_distance = dtw(pp1, pp2);
			printf("\nDTW distance (pp vs pp): %.9f\n", dtw_distance);

			delete[] buffer2;
			delete[] buffer1;

			return dtw_distance;
		}

		//cases:
		float dtw_AmbulancePP_vs_AmbulancePP() {
			float *samples = NULL;
			int nsamp = 0, srate = 0;
			int err;

			std::string ppfile1("C:/Users/dbarbera/Documents/Repositories/Neural_Inference/tests/test-data/ambulance_M.pp");
			std::string ppfile2("C:/Users/dbarbera/Documents/Repositories/Neural_Inference/tests/test-data/ambulance_F.pp");

			float distance = c_dtw_pp_vs_pp(ppfile1.data(), ppfile2.data(), err);

			return distance;
		}

		float dtw_AmbulanceWav_vs_AmbulancePP_AmbulancePP() {
			float *samples = NULL;
			int nsamp = 0, srate = 0;

			std::string wavfile("C:/Users/dbarbera/Documents/Repositories/Neural_Inference/tests/test-data/ambulance_F.wav");
			std::string ppfile1("C:/Users/dbarbera/Documents/Repositories/Neural_Inference/tests/test-data/ambulance_M.pp");
			std::string ppfile2("C:/Users/dbarbera/Documents/Repositories/Neural_Inference/tests/test-data/ambulance_F.pp");
			std::string modelfile("C:/Users/dbarbera/Documents/Repositories/Neural_Inference/tests/test-data/model_weights.bin");

			wav_to_samples(wavfile.data(), (float*&)samples, nsamp, srate);

			float distance = c_dtw_samples_vs_pp_pp(samples, nsamp, srate, ppfile1.data(), ppfile2.data(), modelfile.data());

			return distance;
		}

		//cases:
		float dtw_AmbulanceWav_vs_AmbulancePP() {
			float *samples = NULL;
			int nsamp = 0, srate = 0;
			
			std::string wavfile("C:/Users/dbarbera/Documents/Repositories/Neural_Inference/tests/test-data/ambulance_F.wav");
			std::string ppfile("C:/Users/dbarbera/Documents/Repositories/Neural_Inference/tests/test-data/ambulance_M.pp");
			std::string modelfile("C:/Users/dbarbera/Documents/Repositories/Neural_Inference/tests/test-data/model_weights.bin");

			wav_to_samples(wavfile.data(), (float*&)samples, nsamp, srate);
			
			float distance = c_dtw_samples_vs_pp(samples, nsamp, srate, ppfile.data(), modelfile.data());

			return distance;
		}

		float dtw_Andrew1_vs_Andrew2() {
			std::string fullwav1;
			std::string fullwav2;
			std::string wavdirname("C:/Users/dbarbera/Documents/Repositories/Neural_Inference/");
			std::string wav1("tests/test-data/Henry-1-Andrew.wav");
			std::string wav2("tests/test-data/Henry-2-Andrew.wav");

			fullwav1 = wavdirname + wav1;
			fullwav2 = wavdirname + wav2;
			

			float dtw_distance = dtw_wav(fullwav1.data(), fullwav2.data());
			//printf("\nDTW distance: %.9f", dtw_distance);
			return dtw_distance;
		}

		//cases: passing the model path
		float dtw_Andrew1_vs_Andrew2_model() {
			std::string fullwav1;
			std::string fullwav2;
			std::string fullmodel;

			std::string wavdirname("C:/Users/dbarbera/Documents/Repositories/Neural_Inference/");
			std::string wav1("tests/test-data/Henry-1-Andrew.wav");
			std::string wav2("tests/test-data/Henry-2-Andrew.wav");
			std::string model("tests/test-data/model_weights.bin");

			fullwav1 = wavdirname + wav1;
			fullwav2 = wavdirname + wav2;
			fullmodel = wavdirname + model;


			float dtw_distance = dtw_wav(fullwav1.data(), fullwav2.data(), fullmodel.data());
			//printf("\nDTW distance: %.9f", dtw_distance);
			return dtw_distance;
		}

		float dtw_Andrew1_vs_Andrew2_model_samples_wav() {
			std::string fullsamples1;
			std::string fullwav2;
			std::string fullmodel;

			std::string wavdirname("C:/Users/dbarbera/Documents/Repositories/Neural_Inference/");
			std::string samples1("tests/test-data/22050_andrew1samples.bin");
			std::string wav2("tests/test-data/Henry-2-Andrew.wav");
			std::string model("tests/test-data/model_weights.bin");

			fullsamples1 = wavdirname + samples1;
			fullwav2 = wavdirname + wav2;
			fullmodel = wavdirname + model;

			//Loading PCM data from binary
			std::ifstream infile;
			infile.open(fullsamples1.data(), std::ios::binary | std::ios::in);
			infile.seekg(0, infile.end);
			int infile_length = infile.tellg();
			infile.seekg(0, infile.beg);
			char *buffer_samples1 = new char[infile_length];
			infile.read(buffer_samples1, infile_length * sizeof(char));
			infile.close();
			float* samples = (float*)buffer_samples1;
			int nsamples = infile_length/sizeof(float);
			int srate = 22050;
			bool norm = true;
			float dtw_distance = dtw_samples_wav(samples, nsamples, srate, fullwav2.data(), fullmodel.data(), norm);

			samples = NULL;
			delete[] buffer_samples1;

			return dtw_distance;

		}

		float dtw_Andrew1_vs_Andrew2_model_samples_samples() {
			std::string fullsamples1;
			std::string fullsamples2;
			std::string fullmodel;
			std::string wavdirname("C:/Users/dbarbera/Documents/Repositories/Neural_Inference/");
			std::string path_samples1("tests/test-data/22050_andrew1samples.bin");
			std::string path_samples2("tests/test-data/22050_andrew2samples.bin");
			std::string model("tests/test-data/model_weights.bin");


			fullsamples1 = wavdirname + path_samples1;
			fullsamples2 = wavdirname + path_samples2;
			fullmodel = wavdirname + model;

			//Loading PCM data from binary1
			std::ifstream infile1;
			infile1.open(fullsamples1.data(), std::ios::binary | std::ios::in);
			infile1.seekg(0, infile1.end);
			int infile1_length = infile1.tellg();
			infile1.seekg(0, infile1.beg);
			char *buffer_samples1 = new char[infile1_length];
			infile1.read(buffer_samples1, infile1_length * sizeof(char));
			infile1.close();
			float* samples1 = (float*)buffer_samples1;
			int nsamples1 = infile1_length / sizeof(float);
			int srate1 = 22050;
			bool norm1 = true;

			//Loading PCM data from binary2
			std::ifstream infile2;
			infile2.open(fullsamples2, std::ios::binary | std::ios::in);
			infile2.seekg(0, infile2.end);
			int infile2_length = infile2.tellg();
			infile2.seekg(0, infile2.beg);
			char *buffer_samples2 = new char[infile2_length];
			infile2.read(buffer_samples2, infile2_length * sizeof(char));
			infile2.close();
			float* samples2 = (float*)buffer_samples2;
			int nsamples2 = infile2_length / sizeof(float);
			int srate2 = 22050;
			bool norm2 = true;

			float dtw_distance = dtw_samples_samples(samples1, nsamples1, srate1,  samples2, nsamples2, srate2, fullmodel.data(), norm1, norm2);

			samples1 = NULL;
			samples2 = NULL;
			delete[] buffer_samples1;
			delete[] buffer_samples2;


			return dtw_distance;

		}

	}//namespace gotcha
}//namespace nni