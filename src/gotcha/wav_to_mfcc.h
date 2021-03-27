
namespace nni {
	namespace gotcha {

		int signal_to_mfcc(float* signal, float* mfcc, int nsamp, int srate
			, float twindow_size, float tstep_size
			//, int mfcc_coeff
			//, int mfcc_lofreq, int mfcc_hifreq
			, int* nframes, int* nmfcc) {

			MFCCpy c(nsamp, srate, twindow_size, tstep_size);
			int timesteps = c.get_mfcc(signal, mfcc, nframes, nmfcc);

			//std::ofstream outfile;
			//outfile.open(mfcc_f, std::ios::binary | std::ios::out);
			//outfile.write((char*)rmfcc, (nframes) * (nmfcc) * sizeof(float));
			//outfile.close();

			return timesteps;
		}

		int test_signal_to_mfcc() {

			int nsamp = 0, srate = 0;
			float twindow_size = 0.03, tstep_size = 0.01;
			int _nmfcc = 26;
			//std::string wav_f("duck.wav");
			std::string wav_f("accept-ww-F.wav");

			float* signal;
			float* mfcc;

			wav_to_samples(wav_f.data(), (float*&)signal, nsamp, srate);
			printf("\nnsamples: %i\nsrate: %i\n", nsamp, srate);

			int _window_size = (int)(0.5f + ((float)srate) * ((float)twindow_size))
				;
			int  _step_size = (int)(0.5f + ((float)srate) * ((float)tstep_size));

			assert(_step_size != 0);
			int _nframes = (int)((nsamp - _window_size + _step_size) / _step_size);
			printf("\nwindow: %f seconds -> %i samples\nstep: %f seconds -> %i samples\n", twindow_size, _window_size, tstep_size, _step_size);
			printf("\n%i samples @ %i Hz -> %i frames\n", nsamp, srate, _nframes);

			mfcc = new float[(_nframes)*_nmfcc];

			int nframes = 0;
			int nmfcc = 0;

			int timesteps = signal_to_mfcc(signal, mfcc, nsamp, srate, twindow_size, tstep_size, &nframes, &nmfcc);

			printf("\ntimesteps: %i\nframes: %i\nnmfcc: %i\n", timesteps, nframes, nmfcc);

			std::ofstream outfile;
			outfile.open("duck_mfcc_test.bin", std::ios::binary | std::ios::out);
			outfile.write((char*)mfcc, (nframes) * (nmfcc) * sizeof(float));
			outfile.close();

			delete[] signal;
			delete[] mfcc;

			return 0;
		}

		int wav_to_mfcc(const char* wav_f, const char* mfcc_f) {
			float* samples = NULL;
			int nsamp = 0, srate = 0;

			wav_to_samples(wav_f, (float*&)samples, nsamp, srate);

			float* rmfcc = NULL;
			int nframes = 0;
			int nmfcc = 0;

			MFCC c;
			int timesteps = c.get_mfcc((float*&)samples, nsamp, srate, (float*&)rmfcc, &nframes, &nmfcc);

			std::ofstream outfile;
			outfile.open(mfcc_f, std::ios::binary | std::ios::out);
			outfile.write((char*)rmfcc, (nframes) * (nmfcc) * sizeof(float));
			outfile.close();


			delete[] samples;

			return 0;

		}

	}//namespace gotcha
}//namespace nni