/// (c) Mark Huckvale, adapted to crossplatform code by David Barbera
#include <fstream>
#include "Resample.h"


#ifdef _WIN32
#include <Windows.h>
typedef const TCHAR* LPCTSTR;
typedef const char* _T;
#else
typedef bool BOOL;
typedef char TCHAR;
#endif

namespace nni {
	namespace riffio {

		struct riff_format_rec {
			unsigned short	format;
			unsigned short	nchan;
			unsigned int	srate;
			unsigned int	brate;
			unsigned short	balign;
			unsigned short	nbits;
			unsigned short	xsize;
		};

		int Wav_metadata(const char* filename, int* nsamp, int* nchan, int* srate, int* nbits, int* format)
		{
			using namespace std;

			struct riff_format_rec rifform;
			char	riffstr[8];
			unsigned int	size, numf;
			//unsigned int len;
			filebuf fp;

			memset(&rifform, 0, sizeof(rifform));

			// try to open file
			if (!fp.open(filename, ios_base::in | ios_base::binary)) {
				printf("couldnt load it\n"); return(-7);
			}
			// check is RIFF
			if ((fp.sgetn(riffstr, 4) != 4) || (strncmp(riffstr, "RIFF", 4) != 0) || (fp.sgetn((char*)&size, 4) != 4) || (fp.sgetn(riffstr, 4) != 4) || (strncmp(riffstr, "WAVE", 4) != 0)) {
				fp.close();
				printf("not good format\n");
				return(-2);
			}

			/* look for format chunk */
			do {
				if (fp.sgetn(riffstr, 4) != 4) { fp.close(); printf("BAD header\n"); return(0); }
				if (strncmp(riffstr, "fmt ", 4) == 0) break;
				if (fp.sgetn((char*)&size, 4) != 4) { fp.close(); printf("Wrong size\n"); return(0); }
				fp.pubseekoff(size, ios_base::cur, ios_base::in);
				/*unsigned int pos = (unsigned int)ios_base::cur;
				if ((streamsize)fp.pubseekoff(size, ios_base::cur, ios_base::in) != (LONG)(pos + size)) {
					fp.close();
					return(0);
				}*/
			} while (1);

			/* read format chunk */
			if ((fp.sgetn((char*)&size, 4) != 4) || (size < 16)) {
				fp.close();
				printf("failed reading format chunk\n");
				return(-3);
			}
			if (fp.sgetn((char*)&rifform, 16) != 16) {
				fp.close();
				printf("rifform field incorrect\n");
				return(-4);
			}
			if (size > 16) {
				//unsigned int pos = (unsigned int)ios_base::cur;
				//unsigned int filePos = (streamsize)
				fp.pubseekoff(size - 16, ios_base::cur, ios_base::in);
				/*if ( filePos != (LONG)(pos + size - 16)) {
					fp.close();
					return(0);
				}*/
			}

			/* look for data chunk */
			do {
				if (fp.sgetn(riffstr, 4) != 4) { fp.close(); return(0); }
				if (strncmp(riffstr, "data", 4) == 0) break;
				if (fp.sgetn((char*)&size, 4) != 4) { fp.close(); return(0); }
				fp.pubseekoff(size, ios_base::cur, ios_base::in);
				/*unsigned int pos = (unsigned int)ios_base::cur;
				if ((streamsize)fp.pubseekoff(size, ios_base::cur, ios_base::in) != (LONG)(pos + size)) {
					fp.close();
					return(0);
				}*/
			} while (1);

			/* get # samples to replay */
			if (fp.sgetn((char*)&numf, 4) != 4) { fp.close(); return(0); }

			// calculate how many samples we need to store
			if ((rifform.nbits != 8) && (rifform.nbits != 16) && (rifform.nbits != 24) && (rifform.nbits != 32)) { fp.close(); return(0); }
			int	bpsamp = rifform.nbits / 8;
			*nsamp = numf / (rifform.nchan * bpsamp);
			*nchan = rifform.nchan;
			*srate = rifform.srate;
			*nbits = rifform.nbits;
			*format = rifform.format;

			fp.close();

			int nsamples= numf / (rifform.nchan * bpsamp);

			return nsamples;
		}

		int extract_audio_data(const char* filename, int* nsamp, int* nchan, int* srate, float* pchan1)
		{
			using namespace std;

			struct riff_format_rec rifform;
			char	riffstr[8];
			unsigned int	size, numf;
			//unsigned int len;
			filebuf fp;

			memset(&rifform, 0, sizeof(rifform));

			// try to open file
			if (!fp.open(filename, ios_base::in | ios_base::binary)) {
				printf("couldnt load it\n"); return(-7);
			}
			// check is RIFF
			if ((fp.sgetn(riffstr, 4) != 4) || (strncmp(riffstr, "RIFF", 4) != 0) || (fp.sgetn((char*)&size, 4) != 4) || (fp.sgetn(riffstr, 4) != 4) || (strncmp(riffstr, "WAVE", 4) != 0)) {
				fp.close();
				printf("not good format\n");
				return(-2);
			}

			/* look for format chunk */
			do {
				if (fp.sgetn(riffstr, 4) != 4) { fp.close(); printf("BAD header\n"); return(0); }
				if (strncmp(riffstr, "fmt ", 4) == 0) break;
				if (fp.sgetn((char*)&size, 4) != 4) { fp.close(); printf("Wrong size\n"); return(0); }
				fp.pubseekoff(size, ios_base::cur, ios_base::in);
				/*unsigned int pos = (unsigned int)ios_base::cur;
				if ((streamsize)fp.pubseekoff(size, ios_base::cur, ios_base::in) != (LONG)(pos + size)) {
					fp.close();
					return(0);
				}*/
			} while (1);

			/* read format chunk */
			if ((fp.sgetn((char*)&size, 4) != 4) || (size < 16)) {
				fp.close();
				printf("failed reading format chunk\n");
				return(-3);
			}
			if (fp.sgetn((char*)&rifform, 16) != 16) {
				fp.close();
				printf("rifform field incorrect\n");
				return(-4);
			}
			if (size > 16) {
				//unsigned int pos = (unsigned int)ios_base::cur;
				//unsigned int filePos = (streamsize)
				fp.pubseekoff(size - 16, ios_base::cur, ios_base::in);
				/*if ( filePos != (LONG)(pos + size - 16)) {
					fp.close();
					return(0);
				}*/
			}

			/* look for data chunk */
			do {
				if (fp.sgetn(riffstr, 4) != 4) { fp.close(); return(0); }
				if (strncmp(riffstr, "data", 4) == 0) break;
				if (fp.sgetn((char*)&size, 4) != 4) { fp.close(); return(0); }
				fp.pubseekoff(size, ios_base::cur, ios_base::in);
				/*unsigned int pos = (unsigned int)ios_base::cur;
				if ((streamsize)fp.pubseekoff(size, ios_base::cur, ios_base::in) != (LONG)(pos + size)) {
					fp.close();
					return(0);
				}*/
			} while (1);

			/* get # samples to replay */
			if (fp.sgetn((char*)&numf, 4) != 4) { fp.close(); return(0); }

			// calculate how many samples we need to store
			if ((rifform.nbits != 8) && (rifform.nbits != 16) && (rifform.nbits != 24) && (rifform.nbits != 32)) { fp.close(); return(0); }
			int	bpsamp = rifform.nbits / 8;
			*nsamp = numf / (rifform.nchan * bpsamp);
			*nchan = rifform.nchan;
			*srate = rifform.srate;

			// allocate memory
			//*pchan1 = (float*)malloc(*nsamp * sizeof(float));
			//if (rifform.nchan > 1) *pchan2 = (float*)malloc(*nsamp * sizeof(float));

			// decode into float buffer
			int i, j, len;
			unsigned char xbuf[1536];		/* must be a multiple of 12 bytes ( 3/4 bytes * 2 channels) */

											// read in the data 1/2/3/4 bytes per samp, 1/2 channels
			for (i = 0; i < *nsamp;) {
				if ((len = fp.sgetn((char*)xbuf, sizeof(xbuf))) <= 0) break;
				for (j = 0; (j < len) && (i < *nsamp);) {
					if (bpsamp == 1)
						(pchan1)[i] = (float)((((int)(xbuf[j]) - 128) << 24) / 2147483648.0);
					else if (bpsamp == 2)
						(pchan1)[i] = (float)(((((int)xbuf[j + 1]) << 24) | (((int)xbuf[j]) << 16)) / 2147483648.0);
					else if (bpsamp == 3)
						(pchan1)[i] = (float)(((((int)xbuf[j + 2]) << 24) | (((int)xbuf[j + 1]) << 16) | (((int)xbuf[j]) << 8)) / 2147483648.0);
					else {
						(pchan1)[i] = *(float*)(xbuf + j);
					}
					j += bpsamp;
					//if (rifform.nchan > 1) {
					//	if (bpsamp == 1)
					//		(*pchan2)[i] = (float)((((int)(xbuf[j]) - 128) << 24) / 2147483648.0);
					//	else if (bpsamp == 2)
					//		(*pchan2)[i] = (float)(((((int)xbuf[j + 1]) << 24) | (((int)xbuf[j]) << 16)) / 2147483648.0);
					//	else if (bpsamp == 3)
					//		(*pchan2)[i] = (float)(((((int)xbuf[j + 2]) << 24) | (((int)xbuf[j + 1]) << 16) | (((int)xbuf[j]) << 8)) / 2147483648.0);
					//	else {
					//		(*pchan2)[i] = *(float*)(xbuf + j);
					//	}
					//	j += bpsamp;
					//}
					i++;
				}
			}

			fp.close();

			return rifform.nchan;
		}

		int LoadWav(const char *filename, int * nsamp, int * nchan, int * srate, float ** pchan1, float ** pchan2)
		{
			using namespace std;

			struct riff_format_rec rifform;
			char	riffstr[8];
			unsigned int	size, numf;
			//unsigned int len;
			filebuf fp;

			memset(&rifform, 0, sizeof(rifform));

			// try to open file
			if (!fp.open(filename, ios_base::in | ios_base::binary)) {
				printf("couldnt load it\n"); return(-7);
			}
			// check is RIFF
			if ((fp.sgetn(riffstr, 4) != 4) || (strncmp(riffstr, "RIFF", 4) != 0) || (fp.sgetn((char*)&size, 4) != 4) || (fp.sgetn(riffstr, 4) != 4) || (strncmp(riffstr, "WAVE", 4) != 0)) {
				fp.close();
				printf("not good format\n");
				return(-2);
			}

			/* look for format chunk */
			do {
				if (fp.sgetn(riffstr, 4) != 4) { fp.close(); printf("BAD header\n"); return(0); }
				if (strncmp(riffstr, "fmt ", 4) == 0) break;
				if (fp.sgetn((char*)&size, 4) != 4) { fp.close(); printf("Wrong size\n"); return(0); }
				fp.pubseekoff(size, ios_base::cur, ios_base::in);
				/*unsigned int pos = (unsigned int)ios_base::cur;
				if ((streamsize)fp.pubseekoff(size, ios_base::cur, ios_base::in) != (LONG)(pos + size)) {
					fp.close();
					return(0);
				}*/
			} while (1);

			/* read format chunk */
			if ((fp.sgetn((char*)&size, 4) != 4) || (size < 16)) {
				fp.close();
				printf("failed reading format chunk\n");
				return(-3);
			}
			if (fp.sgetn((char*)&rifform, 16) != 16) {
				fp.close();
				printf("rifform field incorrect\n");
				return(-4);
			}
			if (size > 16) {
				//unsigned int pos = (unsigned int)ios_base::cur;
				//unsigned int filePos = (streamsize)
				fp.pubseekoff(size - 16, ios_base::cur, ios_base::in);
				/*if ( filePos != (LONG)(pos + size - 16)) {
					fp.close();
					return(0);
				}*/
			}

			/* look for data chunk */
			do {
				if (fp.sgetn(riffstr, 4) != 4) { fp.close(); return(0); }
				if (strncmp(riffstr, "data", 4) == 0) break;
				if (fp.sgetn((char*)&size, 4) != 4) { fp.close(); return(0); }
				fp.pubseekoff(size, ios_base::cur, ios_base::in);
				/*unsigned int pos = (unsigned int)ios_base::cur;
				if ((streamsize)fp.pubseekoff(size, ios_base::cur, ios_base::in) != (LONG)(pos + size)) {
					fp.close();
					return(0);
				}*/
			} while (1);

			/* get # samples to replay */
			if (fp.sgetn((char*)&numf, 4) != 4) { fp.close(); return(0); }

			// calculate how many samples we need to store
			if ((rifform.nbits != 8) && (rifform.nbits != 16) && (rifform.nbits != 24) && (rifform.nbits != 32)) { fp.close(); return(0); }
			int	bpsamp = rifform.nbits / 8;
			*nsamp = numf / (rifform.nchan*bpsamp);
			*nchan = rifform.nchan;
			*srate = rifform.srate;

			// allocate memory
			*pchan1 = (float *)malloc(*nsamp * sizeof(float));
			if (rifform.nchan > 1) *pchan2 = (float *)malloc(*nsamp * sizeof(float));

			// decode into float buffer
			int i, j, len;
			unsigned char xbuf[1536];		/* must be a multiple of 12 bytes ( 3/4 bytes * 2 channels) */

											// read in the data 1/2/3/4 bytes per samp, 1/2 channels
			for (i = 0; i < *nsamp;) {
				if ((len = fp.sgetn((char*)xbuf, sizeof(xbuf))) <= 0) break;
				for (j = 0; (j < len) && (i < *nsamp);) {
					if (bpsamp == 1)
						(*pchan1)[i] = (float)((((int)(xbuf[j]) - 128) << 24) / 2147483648.0);
					else if (bpsamp == 2)
						(*pchan1)[i] = (float)(((((int)xbuf[j + 1]) << 24) | (((int)xbuf[j]) << 16)) / 2147483648.0);
					else if (bpsamp == 3)
						(*pchan1)[i] = (float)(((((int)xbuf[j + 2]) << 24) | (((int)xbuf[j + 1]) << 16) | (((int)xbuf[j]) << 8)) / 2147483648.0);
					else {
						(*pchan1)[i] = *(float *)(xbuf + j);
					}
					j += bpsamp;
					if (rifform.nchan > 1) {
						if (bpsamp == 1)
							(*pchan2)[i] = (float)((((int)(xbuf[j]) - 128) << 24) / 2147483648.0);
						else if (bpsamp == 2)
							(*pchan2)[i] = (float)(((((int)xbuf[j + 1]) << 24) | (((int)xbuf[j]) << 16)) / 2147483648.0);
						else if (bpsamp == 3)
							(*pchan2)[i] = (float)(((((int)xbuf[j + 2]) << 24) | (((int)xbuf[j + 1]) << 16) | (((int)xbuf[j]) << 8)) / 2147483648.0);
						else {
							(*pchan2)[i] = *(float *)(xbuf + j);
						}
						j += bpsamp;
					}
					i++;
				}
			}

			fp.close();

			return rifform.nchan;
		}


		int SaveWav(const char *filename, int nsamp, int nchan, int srate, float * chan1, float * chan2)
		{
			using namespace std;
			struct riff_format_rec rifform;
			unsigned int	size;
			filebuf fp;

			memset(&rifform, 0, sizeof(rifform));

			// try to open file
			if (!fp.open(filename, ios_base::out | ios_base::binary | ios_base::trunc)) return(0);

			fp.sputn("RIFF", 4);
			size = 44 + 2 * nchan*nsamp;
			fp.sputn((char*)&size, 4);
			fp.sputn("WAVE", 4);
			fp.sputn("fmt ", 4);
			size = 16;
			fp.sputn((char*)&size, 4);
			rifform.balign = 2 * nchan;
			rifform.brate = 2 * nchan*srate;
			rifform.format = 1;
			rifform.nbits = 16;
			rifform.nchan = nchan;
			rifform.srate = srate;
			rifform.xsize = 0;
			fp.sputn((char*)&rifform, 16);
			fp.sputn("data", 4);
			size = 2 * nchan*nsamp;
			fp.sputn((char*)&size, 4);

			// correct overload if occurs
			float smax1 = 1, smax2 = 1;
			for (int i = 0; i < nsamp; i++) {
				if (chan1[i] > smax1) smax1 = chan1[i];
				if (chan1[i] < -smax1) smax1 = -chan1[i];
				if (nchan == 2) {
					if (chan2[i] > smax2) smax2 = chan2[i];
					if (chan2[i] < -smax2) smax2 = -chan2[i];
				}
			}

			short	sbuf[1024];
			for (int i = 0; i < nsamp;) {
				int j;
				for (j = 0; (j < 1024)&(i < nsamp);) {
					sbuf[j++] = (short)(32767 * chan1[i] / smax1);
					if (nchan == 2) sbuf[j++] = (short)(32767 * chan2[i] / smax2);
					i++;
				}
				fp.sputn((char*)sbuf, 2 * j);
			}

			fp.close();

			return nchan;
		}

		//int wavfile_to_samples(const char* wavfile_, float *rfsp, int &nsamp_, int &srate_) {
		//	
		//}

	}//riffio namespace
}//nni namespace

//#endif