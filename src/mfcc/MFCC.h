/// (c) Mark Huckvale, adapted to crossplatform code by David Barbera
#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <assert.h>

#ifdef _WIN32
#ifdef _WIN64
#define M_PI   3.14159265358979323846264338327950288
#else
#define M_PI   3.14159265358979323846264338327950288
#endif //_WIN64
#endif //_WIN32

#define MFCC_PREEMP	0.95
#define MFCC_LIFTER	23

// MFCC records
#define MFCC_COEFF	12
#define MFCC_SIZE	26
#define MFCC_LOFREQ	100
#define MFCC_HIFREQ	6000
#define MFCC_WINSIZE	0.03 //0.03
#define MFCC_STEPSIZE	0.01
#define VAD_THRESH	1.0			// i.e. 10dB

namespace nni {
	namespace mfcc {
		/* filter bank parameters */
		struct filter_rec {
			int	lidx;	/* low index */
			int	cidx;	/* centre index */
			int	hidx;	/* high index */
			double	*win;	/* FFT weighting window */
		};

		class CMFCC
		{
			struct filter_rec	*m_filtab;
			int					m_nfilter;
			int					m_ncoeff;
			int					m_energy;
			int					m_winsize;
			int					m_fftsize;
			float				*m_wsp;
			float				*m_fsp;
			float				*m_chan;
			float				*m_ac;
			float				*m_kwin;

		protected:
			int designfbank(double lf, double hf, int nfft, int srate);
			double besselI0(double p);
			void KaiserWindow(float *win, int len, double alpha);
			void autocorrelation(float *x, float *rr, int n);

		public:
			CMFCC(int ncoeff, int energy, int winsize, double lofreq, double hifreq, double srate);
			~CMFCC(void);
			void Calc(float *sig, float *mfcc);
			void CalcAMFCC(float *sig, float *mfcc);
		};


		class MFCC {
		public:
			float	expavg[MFCC_SIZE];
			float	gmean[MFCC_SIZE];
			float	gsdev[MFCC_SIZE];
			float	gPCA[MFCC_SIZE * 2];
			// global flags
			int		do_amfcc = 0;			// calculate MFCCs using high order autocorrelations
			int		do_zscore = 0;		// normalise data with z-scores
			int		do_addPCA = 0;		// convert z-scores into a few PCA coefficients on each vector
			int		do_zread = 0;			// read zscore mean and stddev from file rather than calculate them
			//Data structures:
			struct mfcc_rec {
				int		ssamp;				// start sample of analysis window
				int		esamp;				// end sample of analysis window
				char	 phone[4];			// assigned phone
				float	vad;				// voice activity
				float	data[MFCC_SIZE];	// MFCC coeffs + deltas
			};
			struct mfcc_rec	*mdata = NULL;
			int	mcount = 0;
			int nb_frames = 0;
			float *mfcc = NULL;

			// calculate the MFCC coefficients
			int CalcMFCC(float *sig, int nsamp, int srate)
			{
				int		winsize = (int)(0.5 + srate * MFCC_WINSIZE);
				int		stpsize = (int)(0.5 + srate * MFCC_STEPSIZE);
				CMFCC	mfcc(MFCC_COEFF, 1, winsize, MFCC_LOFREQ, MFCC_HIFREQ, srate);
				int		nframe = ((nsamp - winsize + stpsize) / stpsize);
				nb_frames = nframe;
				// get memory table
				mdata = (struct mfcc_rec *)calloc(nframe + 1, sizeof(struct mfcc_rec));

				// get mfcc coefficients
				mcount = 0;
				for (int i = 0; i < nsamp - winsize; i += stpsize) {
					mdata[mcount].ssamp = i;
					mdata[mcount].esamp = i + winsize - 1;
					if (do_amfcc)
						mfcc.CalcAMFCC(sig + i, mdata[mcount].data);
					else
						mfcc.Calc(sig + i, mdata[mcount].data);
					mcount++;
				}

				// calculate deltas
				if (mcount > 2) {
					for (int i = 1; i < mcount - 1; i++) {
						for (int j = 0; j <= MFCC_COEFF; j++)
							mdata[i].data[j + MFCC_COEFF + 1] = (mdata[i + 1].data[j] - mdata[i - 1].data[j]) / 2;
					}
					for (int j = 0; j <= MFCC_COEFF; j++) {
						mdata[0].data[j + MFCC_COEFF + 1] = mdata[1].data[j + MFCC_COEFF + 1];
						mdata[mcount - 1].data[j + MFCC_COEFF + 1] = mdata[mcount - 2].data[j + MFCC_COEFF + 1];
					}
				}

				return(mcount);
			}

			float	minlevel, curlevel;

			void SpeechActivity()
			{
				int	i;
				float	alpha;

				for (i = 0; i < mcount; i++) {
					float level = mdata[i].data[MFCC_COEFF];
					// update minimum level
					if (level < minlevel) {
						// follow rapidly going down
						alpha = 2 / (1 + 0.05 / MFCC_STEPSIZE);
					}
					else {
						// follow slowly going up
						alpha = 2 / (1 + 10.0 / MFCC_STEPSIZE);
					}
					minlevel = (1 - alpha)*minlevel + alpha * level;
					// update current level
					if (level > curlevel) {
						// follow rapidly going up
						alpha = 2 / (1 + 0.05 / MFCC_STEPSIZE);
					}
					else {
						// follow slowly going down
						alpha = 2 / (1 + 0.1 / MFCC_STEPSIZE);
					}
					curlevel = (1 - alpha)*curlevel + alpha * level;
					// make decision
					mdata[i].vad = (curlevel - minlevel);
					//		printf("lev=%g min=%g cur=%g vad=%g\n",level,minlevel,curlevel,curlevel-minlevel);
				}
			}

			// Normalise MFCC vectors
			void Normalise()
			{
				int	i, j;
				float	alpha = 2.0 / (1 + 10 / MFCC_STEPSIZE);			// exponential average over 10ms

				if (do_zscore) {
					if (do_zread) {
						for (i = 0; i < mcount; i++) {
							for (j = 0; j < MFCC_SIZE; j++)
								mdata[i].data[j] = (float)((mdata[i].data[j] - gmean[j]) / gsdev[j]);
						}
					}
					else {
						for (j = 0; j < MFCC_SIZE; j++) {
							double sum = 0;
							double sumsq = 0;
							int cnt = 0;
							for (i = 0; i < mcount; i++) if (mdata[i].vad > VAD_THRESH) {
								sum += mdata[i].data[j];
								sumsq += mdata[i].data[j] * mdata[i].data[j];
								cnt++;
							}
							if (cnt > 0) {
								double mean = sum / cnt;
								double sdev = sqrt(cnt*sumsq - sum * sum) / cnt; //android/windows issue std::
								for (i = 0; i < mcount; i++)
									mdata[i].data[j] = (float)((mdata[i].data[j] - mean) / sdev);
								gmean[j] = mean;
								gsdev[j] = sdev;
							}
						}
					}
					if (do_addPCA) {
						/*for (i = 0; i < do_addPCA; i++) {
							double pca = 0;
							for (j = 0; j < MFCC_SIZE; j++) {
								pca += (gmean[j] - PCAzscore[0][2 * j])*PCAzscore[i + 1][2 * j];
								pca += (gsdev[j] - PCAzscore[0][2 * j + 1])*PCAzscore[i + 1][2 * j + 1];
							}
							gPCA[i] = (float)pca;
						}*/
					}
				}
				else {
					for (i = 0; i < mcount; i++) {
						if (mdata[i].vad > VAD_THRESH) {
							for (j = 0; j < MFCC_SIZE; j++)
								expavg[j] = (1 - alpha)*expavg[j] + alpha * mdata[i].data[j];
						}
						for (j = 0; j < MFCC_SIZE; j++) {
							mdata[i].data[j] -= expavg[j];
						}
					}
				}
			}

			//------- Signal Processing part -------------------------------------
			// compute MFCC coefficients
			int get_mfcc(float *fsp, int nsamp, int srate, float* &mfcc, int *nframes, int *nmfcc) {
				CalcMFCC(fsp, nsamp, srate);

				// Compute speech activity
				minlevel = curlevel = mdata[0].data[MFCC_COEFF];
				SpeechActivity();

				// Normalise
				for (int j = 0; j < MFCC_SIZE; j++) expavg[j] = 0;
				Normalise();

				mfcc = new float[(nb_frames) * MFCC_SIZE];
				//printf("mfcc nb_frames: %i\n", nb_frames);
				char* char_mfcc = new char[(nb_frames)* MFCC_SIZE*sizeof(float)];
				char* char_num = new char[sizeof(float)];


				for (int i = 0; i < nb_frames; i++) {
					for (int j = 0; j < MFCC_SIZE; j++) //{
						//memcpy(&mfcc[i*MFCC_SIZE], &mdata[i].data[0], sizeof(float)*MFCC_SIZE);
						mfcc[i*MFCC_SIZE + j] = mdata[i].data[j];
						//float num = mdata[i].data[j];
						//char_num = (char*)&num;
						//for (int k = 0; k < sizeof(float); k++) {
						//	char_mfcc[i*MFCC_SIZE + j * sizeof(float) + k] = (char)(num+k);

						//}
						//printf("\n(%i,%i) \t\t%.9f\t\t%.9f", i, j, mfcc[i*MFCC_SIZE + j], mdata[i].data[j]);

					//}
				}
				*nmfcc = MFCC_SIZE;
				*nframes = nb_frames;

				//printf("*nframes: %i \n *nmfcc: %i \n", *nframes, *nmfcc);
				//std::ofstream outfile;
				//outfile.open("tests/test-data/duck_mfcc_just_processed.bin", std::ios::binary | std::ios::out);
				//outfile.write((char*)mfcc, (nb_frames)*(MFCC_SIZE) * sizeof(float));
				//outfile.close();

				delete[] char_mfcc;
				delete[] char_num;

				return nb_frames;
			}

			~MFCC() {
				delete[] mdata;
				delete[] mfcc;
			}

		}; //Class MFCC


		class MFCCpy {
		public:
			float	expavg[MFCC_SIZE];
			float	gmean[MFCC_SIZE];
			float	gsdev[MFCC_SIZE];
			float	gPCA[MFCC_SIZE * 2];
			// global flags
			int		do_amfcc = 0;			// calculate MFCCs using high order autocorrelations
			int		do_zscore = 0;		// normalise data with z-scores
			int		do_addPCA = 0;		// convert z-scores into a few PCA coefficients on each vector
			int		do_zread = 0;			// read zscore mean and stddev from file rather than calculate them
			//Data structures:
			struct mfcc_rec {
				int		ssamp;				// start sample of analysis window
				int		esamp;				// end sample of analysis window
				char	 phone[4];			// assigned phone
				float	vad;				// voice activity
				float	data[MFCC_SIZE];	// MFCC coeffs + deltas
			};
			struct mfcc_rec* mdata = NULL;
			int	mcount = 0;
			int nb_frames = 0;
			//float* mfcc = NULL;
			//Breaking signal:
			float _twindow_size = MFCC_WINSIZE;
			float _tstep_size = MFCC_STEPSIZE;
			int _window_size = 480;
			int _step_size = 160;
			int _nsamp = 0;
			int _srate = 16000;
			int _mfcc_lofreq = MFCC_LOFREQ;
			int _mfcc_hifreq = MFCC_HIFREQ;
			int _nframes = 0;
			int _mfcc_coeff = MFCC_COEFF;

			MFCCpy(int nsamp, int srate=16000, float twindow_size=MFCC_WINSIZE, float tstep_size=MFCC_STEPSIZE, int mfcc_coeff=MFCC_COEFF, int mfcc_lofreq=MFCC_LOFREQ, int mfcc_hifreq=MFCC_HIFREQ) {
				_twindow_size = twindow_size;
				_tstep_size = tstep_size;
				_srate = srate;
				_nsamp = nsamp;

				_window_size = (int)(0.5f+((float)_srate)* ((float)_twindow_size));
				_step_size = (int)(0.5f + ((float)_srate) * ((float)_tstep_size));
				
				
				//printf("\n tstep: %f  step: %i", _tstep_size, _step_size);
				assert(_step_size != 0);
				_nframes = (int)((_nsamp - _window_size + _step_size) / _step_size);

				_mfcc_lofreq = mfcc_lofreq;
				_mfcc_hifreq = mfcc_hifreq;
				_mfcc_coeff = mfcc_coeff;
			}

			// calculate the MFCC coefficients
			int CalcMFCC(float* sig)
			{
			    CMFCC	mfcc(_mfcc_coeff, 1, _window_size, _mfcc_lofreq, _mfcc_hifreq, _srate);

				// get memory table
				assert(_nframes != 0);
				mdata = (struct mfcc_rec*)calloc(_nframes + (int)1, sizeof(struct mfcc_rec));

				// get mfcc coefficients
				mcount = 0;
				for (int i = 0; i < _nsamp - _window_size; i += _step_size) {
					mdata[mcount].ssamp = i;
					mdata[mcount].esamp = i + _window_size - 1;
					if (do_amfcc)
						mfcc.CalcAMFCC(sig + i, mdata[mcount].data);
					else
						mfcc.Calc(sig + i, mdata[mcount].data);
					mcount++;
				}

				// calculate deltas
				if (mcount > 2) {
					for (int i = 1; i < mcount - 1; i++) {
						for (int j = 0; j <= _mfcc_coeff; j++)
							mdata[i].data[j + _mfcc_coeff + 1] = (mdata[i + 1].data[j] - mdata[i - 1].data[j]) / 2;
					}
					for (int j = 0; j <= _mfcc_coeff; j++) {
						mdata[0].data[j + _mfcc_coeff + 1] = mdata[1].data[j + _mfcc_coeff + 1];
						mdata[mcount - 1].data[j + _mfcc_coeff + 1] = mdata[mcount - 2].data[j + _mfcc_coeff + 1];
					}
				}

				return(mcount);
			}

			float	minlevel, curlevel;

			void SpeechActivity()
			{
				int	i;
				float	alpha;

				for (i = 0; i < mcount; i++) {
					float level = mdata[i].data[_mfcc_coeff];
					// update minimum level
					if (level < minlevel) {
						// follow rapidly going down
						alpha = 2 / (1 + 0.05 / _tstep_size);
					}
					else {
						// follow slowly going up
						alpha = 2 / (1 + 10.0 / _tstep_size);
					}
					minlevel = (1 - alpha) * minlevel + alpha * level;
					// update current level
					if (level > curlevel) {
						// follow rapidly going up
						alpha = 2 / (1 + 0.05 / _tstep_size);
					}
					else {
						// follow slowly going down
						alpha = 2 / (1 + 0.1 / _tstep_size);
					}
					curlevel = (1 - alpha) * curlevel + alpha * level;
					// make decision
					mdata[i].vad = (curlevel - minlevel);
					//		printf("lev=%g min=%g cur=%g vad=%g\n",level,minlevel,curlevel,curlevel-minlevel);
				}
			}

			// Normalise MFCC vectors
			void Normalise()
			{
				int	i, j;
				float	alpha = 2.0 / (1 + 10 / _tstep_size);			// exponential average over 10ms

				if (do_zscore) {
					if (do_zread) {
						for (i = 0; i < mcount; i++) {
							for (j = 0; j < MFCC_SIZE; j++)
								mdata[i].data[j] = (float)((mdata[i].data[j] - gmean[j]) / gsdev[j]);
						}
					}
					else {
						for (j = 0; j < MFCC_SIZE; j++) {
							double sum = 0;
							double sumsq = 0;
							int cnt = 0;
							for (i = 0; i < mcount; i++) if (mdata[i].vad > VAD_THRESH) {
								sum += mdata[i].data[j];
								sumsq += mdata[i].data[j] * mdata[i].data[j];
								cnt++;
							}
							if (cnt > 0) {
								double mean = sum / cnt;
								double sdev = sqrt(cnt * sumsq - sum * sum) / cnt; //android/windows issue std::
								for (i = 0; i < mcount; i++)
									mdata[i].data[j] = (float)((mdata[i].data[j] - mean) / sdev);
								gmean[j] = mean;
								gsdev[j] = sdev;
							}
						}
					}
					if (do_addPCA) {
						/*for (i = 0; i < do_addPCA; i++) {
							double pca = 0;
							for (j = 0; j < MFCC_SIZE; j++) {
								pca += (gmean[j] - PCAzscore[0][2 * j])*PCAzscore[i + 1][2 * j];
								pca += (gsdev[j] - PCAzscore[0][2 * j + 1])*PCAzscore[i + 1][2 * j + 1];
							}
							gPCA[i] = (float)pca;
						}*/
					}
				}
				else {
					for (i = 0; i < mcount; i++) {
						if (mdata[i].vad > VAD_THRESH) {
							for (j = 0; j < MFCC_SIZE; j++)
								expavg[j] = (1 - alpha) * expavg[j] + alpha * mdata[i].data[j];
						}
						for (j = 0; j < MFCC_SIZE; j++) {
							mdata[i].data[j] -= expavg[j];
						}
					}
				}
			}

			//------- Signal Processing part -------------------------------------
			// compute MFCC coefficients
			int get_mfcc(float* fsp, float* mfcc, int* nframes, int* nmfcc) {
				CalcMFCC(fsp);

				// Compute speech activity
				minlevel = curlevel = mdata[0].data[_mfcc_coeff];
				SpeechActivity();

				// Normalise
				for (int j = 0; j < MFCC_SIZE; j++) expavg[j] = 0;
				Normalise();

				//mfcc = new float[(nb_frames)*MFCC_SIZE];
				//printf("mfcc nb_frames: %i\n", nb_frames);
				//char* char_mfcc = new char[(nb_frames)*MFCC_SIZE * sizeof(float)];
				//char* char_num = new char[sizeof(float)];


				for (int i = 0; i < _nframes; i++) {
					for (int j = 0; j < MFCC_SIZE; j++) //{
						//memcpy(&mfcc[i*MFCC_SIZE], &mdata[i].data[0], sizeof(float)*MFCC_SIZE);
						mfcc[i * MFCC_SIZE + j] = mdata[i].data[j];
					//float num = mdata[i].data[j];
					//char_num = (char*)&num;
					//for (int k = 0; k < sizeof(float); k++) {
					//	char_mfcc[i*MFCC_SIZE + j * sizeof(float) + k] = (char)(num+k);

					//}
					//printf("\n(%i,%i) \t\t%.9f\t\t%.9f", i, j, mfcc[i*MFCC_SIZE + j], mdata[i].data[j]);

				//}
				}
				*nmfcc = MFCC_SIZE;
				*nframes = _nframes;

				//printf("*nframes: %i \n *nmfcc: %i \n", *nframes, *nmfcc);
				//std::ofstream outfile;
				//outfile.open("tests/test-data/duck_mfcc_just_processed.bin", std::ios::binary | std::ios::out);
				//outfile.write((char*)mfcc, (nb_frames)*(MFCC_SIZE) * sizeof(float));
				//outfile.close();

				//delete[] char_mfcc;
				//delete[] char_num;

				return _nframes;
				}

			~MFCCpy() {
				delete[] mdata;
				//delete[] mfcc;
			}

		}; //Class MFCCpy




	}
}

//#endif


