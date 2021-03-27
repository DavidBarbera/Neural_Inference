
namespace nni {

	namespace plugin {
		//Tests Unity plugin usage

		int can_load_wav_file(const char* wav) {
			//const char *wav = "tests/test-data/Henry-1-Andrew.wav";
			//const char *abswav = "C:/Users/dbarbera/AppData/LocalLow/DefaultCompany/Gotcha Demo/tests/tests/test-data/Henry-1-Andrew.wav";
			int	nsamp, srate;
			//float *rfsp;

			int nchan, orate, rnsamp;
			float	*fsp, *fsp2;

			//const char wavfile[] = "tests/test-data/duck.wav";
			//const char * wavfile = "C:/Users/dbarbera/test.wav";
			int error = LoadWav(wav, &nsamp, &nchan, &srate, &fsp, &fsp2);
			if (error < 0) {
				//fprintf(stderr, "failed to read WAV ond test \"test_from_wav_to_samples\"");
				//getchar();
				return error;
			}

			if (nchan == 2) free(fsp2);
			if (nchan == 1) free(fsp);

			return 0;
		}
	}//namespace plugin
}//namespace nni