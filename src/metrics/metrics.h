// align algorithm from Mark Huckvale

#include "alphabet.h"

using namespace std::chrono;


using namespace boost::interprocess;


namespace nni {
	namespace metrics {

		float maxof3(float a, float b, float c)
		{
			if (a > b) {
				if (a > c)
					return(a);
				else
					return(c);
			}
			else {
				if (b > c)
					return(b);
				else
					return(c);
			}
		}

		int maxof3idx(float a, float b, float c)
		{
			if (a > b) {
				if (a > c)
					return(1);
				else
					return(3);
			}
			else {
				if (b > c)
					return(2);
				else
					return(3);
			}
		}
		//Eigen API without index
		template <typename DerivedPh, typename DerivedPo>
		float Align(const Eigen::MatrixBase<DerivedPh>& phones, int nphone, const Eigen::MatrixBase<DerivedPo>& ptable, int pcount)
		{
			float	**htable = (float **)calloc(pcount, sizeof(float *));
			int		**btable = (int **)calloc(pcount, sizeof(int *));

			// allocate memory
			for (int i = 0; i < pcount; i++) {
				htable[i] = (float *)calloc(nphone, sizeof(float));
				btable[i] = (int *)calloc(nphone, sizeof(int));
			}

			// viterbi aligner
#define TINY	-1.0E10;
			float	d, u, l;
			for (int i = 0; i < pcount; i++) {
				for (int j = 0; j < nphone; j++) {
					if ((i == 0) && (j == 0)) {
						d = 0;
						l = TINY;
						u = TINY;
					}
					else if (i == 0) {
						d = TINY;
						l = htable[i][j - 1];
						u = TINY;
					}
					else if (j == 0) {
						d = TINY;
						l = TINY;
						u = htable[i - 1][j];
					}
					else {
						d = htable[i - 1][j - 1];
						l = htable[i][j - 1];
						u = htable[i - 1][j];
					}
					htable[i][j] = log(ptable(i, phones(j))) + maxof3(d, l, u);
					btable[i][j] = maxof3idx(d, l, u);
				}
			}

			// run backwards to get best path
			int i = pcount - 1;
			int j = nphone - 1;
			float finalscore = htable[i][j];

			////If pindex (best path) is required:
			//while (i >= 0) {
			//	pindex[i] = phones(j);
			//	switch (btable[i][j]) {
			//	case 1:	// diagonal
			//		i--;
			//		j--;
			//		break;
			//	case 2: // left
			//		j--;
			//		break;
			//	case 3: // up
			//		i--;
			//		break;
			//	}
			//}

			// free memory
			for (int i = 0; i < pcount; i++) {
				free(htable[i]);
				free(btable[i]);
			}
			free(htable);
			free(btable);

			// return final likelihood
			return finalscore;
		}

		template <typename Derived>
		float Align(int *phones, int nphone, const Eigen::MatrixBase<Derived>& ptable, int *pindex, int pcount)
		{
			float	**htable = (float **)calloc(pcount, sizeof(float *));
			int		**btable = (int **)calloc(pcount, sizeof(int *));

			// allocate memory
			for (int i = 0; i < pcount; i++) {
				htable[i] = (float *)calloc(nphone, sizeof(float));
				btable[i] = (int *)calloc(nphone, sizeof(int));
			}

			// viterbi aligner
#define TINY	-1.0E10;
			float	d, u, l;
			for (int i = 0; i < pcount; i++) {
				for (int j = 0; j < nphone; j++) {
					if ((i == 0) && (j == 0)) {
						d = 0;
						l = TINY;
						u = TINY;
					}
					else if (i == 0) {
						d = TINY;
						l = htable[i][j - 1];
						u = TINY;
					}
					else if (j == 0) {
						d = TINY;
						l = TINY;
						u = htable[i - 1][j];
					}
					else {
						d = htable[i - 1][j - 1];
						l = htable[i][j - 1];
						u = htable[i - 1][j];
					}
					htable[i][j] = log(ptable(i,phones[j])) + maxof3(d, l, u);
					btable[i][j] = maxof3idx(d, l, u);
				}
			}

			// run backwards to get best path
			int i = pcount - 1;
			int j = nphone - 1;
			float finalscore = htable[i][j];
			while (i >= 0) {
				pindex[i] = phones[j];
				switch (btable[i][j]) {
				case 1:	// diagonal
					i--;
					j--;
					break;
				case 2: // left
					j--;
					break;
				case 3: // up
					i--;
					break;
				}
			}

			// free memory
			for (int i = 0; i < pcount; i++) {
				free(htable[i]);
				free(btable[i]);
			}
			free(htable);
			free(btable);

			// return final likelihood
			return finalscore;
		}

		float Align(int *phones, int nphone, float **ptable, int *pindex, int pcount)
		{
			float	**htable = (float **)calloc(pcount, sizeof(float *));
			int		**btable = (int **)calloc(pcount, sizeof(int *));

			// allocate memory
			for (int i = 0; i < pcount; i++) {
				htable[i] = (float *)calloc(nphone, sizeof(float));
				btable[i] = (int *)calloc(nphone, sizeof(int));
			}

			// viterbi aligner
#define TINY	-1.0E10;
			float	d, u, l;
			for (int i = 0; i < pcount; i++) {
				for (int j = 0; j < nphone; j++) {
					if ((i == 0) && (j == 0)) {
						d = 0;
						l = TINY;
						u = TINY;
					}
					else if (i == 0) {
						d = TINY;
						l = htable[i][j - 1];
						u = TINY;
					}
					else if (j == 0) {
						d = TINY;
						l = TINY;
						u = htable[i - 1][j];
					}
					else {
						d = htable[i - 1][j - 1];
						l = htable[i][j - 1];
						u = htable[i - 1][j];
					}
					htable[i][j] = log(ptable[i][ phones[j]]) + maxof3(d, l, u);
					btable[i][j] = maxof3idx(d, l, u);
				}
			}

			// run backwards to get best path
			int i = pcount - 1;
			int j = nphone - 1;
			float finalscore = htable[i][j];
			while (i >= 0) {
				pindex[i] = phones[j];
				switch (btable[i][j]) {
				case 1:	// diagonal
					i--;
					j--;
					break;
				case 2: // left
					j--;
					break;
				case 3: // up
					i--;
					break;
				}
			}

			// free memory
			for (int i = 0; i < pcount; i++) {
				free(htable[i]);
				free(btable[i]);
			}
			free(htable);
			free(btable);

			// return final likelihood
			return finalscore;
		}

		float Align(int *phones, int nphone, float *buffer, int stride, int *pindex, int pcount)
		{
			float	**htable = (float **)calloc(pcount, sizeof(float *));
			int		**btable = (int **)calloc(pcount, sizeof(int *));

			// allocate memory
			for (int i = 0; i < pcount; i++) {
				htable[i] = (float *)calloc(nphone, sizeof(float));
				btable[i] = (int *)calloc(nphone, sizeof(int));
			}

			// viterbi aligner
#define TINY	-1.0E10;
			float	d, u, l;
			for (int i = 0; i < pcount; i++) {
				for (int j = 0; j < nphone; j++) {
					if ((i == 0) && (j == 0)) {
						d = 0;
						l = TINY;
						u = TINY;
					}
					else if (i == 0) {
						d = TINY;
						l = htable[i][j - 1];
						u = TINY;
					}
					else if (j == 0) {
						d = TINY;
						l = TINY;
						u = htable[i - 1][j];
					}
					else {
						d = htable[i - 1][j - 1];
						l = htable[i][j - 1];
						u = htable[i - 1][j];
					}
					htable[i][j] = log( *(buffer+i*stride+phones[j]) )+ maxof3(d, l, u);
					btable[i][j] = maxof3idx(d, l, u);
				}
			}

			// run backwards to get best path
			int i = pcount - 1;
			int j = nphone - 1;
			float finalscore = htable[i][j];
			while (i >= 0) {
				pindex[i] = phones[j];
				switch (btable[i][j]) {
				case 1:	// diagonal
					i--;
					j--;
					break;
				case 2: // left
					j--;
					break;
				case 3: // up
					i--;
					break;
				}
			}

			// free memory
			for (int i = 0; i < pcount; i++) {
				free(htable[i]);
				free(btable[i]);
			}
			free(htable);
			free(btable);

			// return final likelihood
			return finalscore;
		}

		void align_from_file_eigen(const char* ppfile1) {
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

			//method: eigen alignment
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> pp1(posteriors1, T1, out);

			int phones[5]={ 0, 2, 31, 24, 0 };
			int nphones = 5;

			int *tindex = (int *)calloc(T1, sizeof(int));
			high_resolution_clock::time_point start;
			high_resolution_clock::time_point end;

			

			start = high_resolution_clock::now();
			float tscore = Align(phones, nphones, pp1, tindex, T1);
			end = high_resolution_clock::now();

			auto dur_us = duration<double, std::micro>(end - start).count();
			auto dur_ms = duration<double, std::milli>(end - start).count();
			printf("Time: %lfus %lfms\n", dur_us, dur_ms);

			printf("pp.shape: (%i,%i)\n", T1, out);
			printf("tscore: %f\n", tscore);

			//printf("[");
			//for (int i = 0; i < T1; i++) {
			//	printf(" %i",tindex[i]);
			//}
			//printf(" ]");
				
			free(tindex);
			delete[] buffer1;

		}

		void align_from_file_calloc(const char* ppfile1) {
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

			int phones[5] = { 0, 2, 31, 24, 0 };
			int nphones = 5;

			//method: calloc pointers
			float **pp1 = (float**)calloc(T1, sizeof(float*));
			// allocate memory
			for (int i = 0; i < T1; i++) {
				pp1[i] = (float *)calloc(out, sizeof(float));
				
			}

			for (int i = 0; i < T1; i++) {
				for (int j = 0; j < out; j++) {
					pp1[i][j] = *(posteriors1 + i * out + j);
					//printf("%.9e ",pp1[i][j]);
				}
				//printf("\n");
			}

			int *tindex = (int *)calloc(T1, sizeof(int));
			high_resolution_clock::time_point start;
			high_resolution_clock::time_point end;



			start = high_resolution_clock::now();
			float tscore = Align(phones, nphones, pp1, tindex, T1);
			end = high_resolution_clock::now();

			auto dur_us = duration<double, std::micro>(end - start).count();
			auto dur_ms = duration<double, std::milli>(end - start).count();
			printf("Time: %lfus %lfms\n", dur_us, dur_ms);

			printf("pp.shape: (%i,%i)\n", T1, out);
			printf("tscore: %f\n", tscore);

			//printf("[");
			//for (int i = 0; i < T1; i++) {
			//	printf(" %i", tindex[i]);
			//}
			//printf(" ]");

			// free memory for calloc pointers
			for (int i = 0; i < T1; i++) {
				free(pp1[i]);
			}
			free(pp1);

			free(tindex);
			delete[] buffer1;

		}

		void align_from_file_buffer(const char* ppfile1) {
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

			int phones[5] = { 0, 2, 31, 24, 0 };
			int nphones = 5;


			int *tindex = (int *)calloc(T1, sizeof(int));
			high_resolution_clock::time_point start;
			high_resolution_clock::time_point end;



			start = high_resolution_clock::now();
			float tscore = Align(phones, nphones, posteriors1, out, tindex, T1);
			end = high_resolution_clock::now();

			auto dur_us = duration<double, std::micro>(end - start).count();
			auto dur_ms = duration<double, std::milli>(end - start).count();
			printf("Time: %lfus %lfms\n", dur_us, dur_ms);

			printf("pp.shape: (%i,%i)\n", T1, out);
			printf("tscore: %f\n", tscore);

			//printf("[");
			//for (int i = 0; i < T1; i++) {
			//	printf(" %i", tindex[i]);
			//}
			//printf(" ]");

			free(tindex);
			delete[] buffer1;

		}


		void test_align() {
			std::string ppfile("C:/Users/dbarbera/Documents/Repositories/Neural_Inference/tests/test-data/apple_F.pp");

			high_resolution_clock::time_point start;
			high_resolution_clock::time_point end;

			start = high_resolution_clock::now();
			align_from_file_eigen(ppfile.data());
			end = high_resolution_clock::now();

			auto dur_us = duration<double, std::micro>(end - start).count();
			auto dur_ms = duration<double, std::milli>(end - start).count();
			printf("\nTime Eigen Method: %lfus %lfms\n\n", dur_us, dur_ms);

			start = high_resolution_clock::now();
			align_from_file_calloc(ppfile.data());
			end = high_resolution_clock::now();

			dur_us = duration<double, std::micro>(end - start).count();
			dur_ms = duration<double, std::milli>(end - start).count();
			printf("\nTime Calloc Method: %lfus %lfms\n\n", dur_us, dur_ms);

			start = high_resolution_clock::now();
			align_from_file_buffer(ppfile.data());
			end = high_resolution_clock::now();

			dur_us = duration<double, std::micro>(end - start).count();
			dur_ms = duration<double, std::milli>(end - start).count();
			printf("\nTime Buffer Method: %lfus %lfms\n\n", dur_us, dur_ms);


		}

		void reset_transcript(uint8_t *reset, uint8_t *transcript, int len) {
			std::memcpy(transcript, reset, sizeof(uint8_t)*len);
		}

		struct score {
			float tscore;
			uint8_t transcript[256];
		};

		void flush_to_file(std::ofstream& file, score *beam, int beam_size) {
			for (int i = 0; i < beam_size; i++) {
				file.write((char*)&beam[i].tscore, sizeof(float));
				file.write((char*)&(beam[i].transcript[0]), sizeof(uint8_t) * 256);
			}
		}

		void test_data_storage() {
			std::string binfile("C:/Users/dbarbera/Documents/Repositories/Neural_Inference/tests/test-data/data_best10.bin");

			float tscore = -10.3456f;
			uint8_t transcript[256];
			uint8_t reset[256];

			for (int i = 0; i < 256; i++) {
				reset[i] = 0;
			}


			reset_transcript(&reset[0], &transcript[0], 256);
			//printf(transcript)
			transcript[0] = 2;
			transcript[1] = 21;
			transcript[3] = 34;



			score beam[10];

			for (int i = 0; i < 10; i++) {
				beam[i].tscore = (float)i;
				reset_transcript(&reset[0],&(beam[i].transcript[0]),256);
				beam[i].transcript[0] = i;
				beam[i].transcript[1] = i;
				beam[i].transcript[2] = i;
			}

			std::ofstream outfile;

			outfile=std::ofstream(binfile.data(), std::ios::binary | std::ios::out);



			flush_to_file(outfile, &beam[0],10);

			outfile.close();
		}

		void test_pp_align_on_dictionary() {

			using namespace boost::interprocess;

			const char* FileName = "C:/Users/dbarbera/Documents/Repositories/Neural_Inference/tests/test-data/posteriors.bin";
			std::string posteriorsFile("C:/Users/dbarbera/Documents/Repositories/Neural_Inference/tests/test-data/posteriors.bin");
			//std::string posteriorsFile("C:/Users/dbarbera/Documents/Repositories/Neural_Inference/tests/test-data/posteriors2.bin");
			//std::string posteriorsFile("C:/Users/dbarbera/Documents/Repositories/Neural_Inference/tests/test-data/pp.bin");
			//std::string posteriorsFile("C:/Users/dbarbera/Documents/Repositories/Neural_Inference/tests/test-data/ppx2.bin");
			std::string posteriorsOffsetsFile("C:/Users/dbarbera/Documents/Repositories/Neural_Inference/tests/test-data/posteriors_offsets.bin");
			std::string posteriorsLengthsFile("C:/Users/dbarbera/Documents/Repositories/Neural_Inference/tests/test-data/posteriors_lengths.bin");

			//	std::string dictionaryFile("C:/Users/dbarbera/Documents/Repositories/Neural_Inference/tests/test-data/dictionary_dense_int8.bin");
			std::string dictionaryFile("C:/Users/dbarbera/Documents/Repositories/Neural_Inference/tests/test-data/dictionary_dense_int8_sils_clean2.bin");
			//std::string dictionaryFile("C:/Users/dbarbera/Documents/Repositories/Neural_Inference/tests/test-data/dictionary_dense_int.bin");
			
			std::string dictionaryOffsetsFile("C:/Users/dbarbera/Documents/Repositories/Neural_Inference/tests/test-data/dictionary_dense_offsets_sils_clean2.bin");
			std::string dictionaryLengthsFile("C:/Users/dbarbera/Documents/Repositories/Neural_Inference/tests/test-data/dictionary_dense_lengths_sils_clean2.bin");

			std::ifstream fposteriors, fposteriors_offsets, fposteriors_lengths;
			std::ifstream fdictionary, fdictionary_offsets, fdictionary_lengths;

			//Posteriors
			//fposteriors.open(posteriorsFile.data(), std::ios::binary | std::ios::in);
			//fposteriors.seekg(0, fposteriors.end);
			//int fposteriors_length = fposteriors.tellg();
			//fposteriors.seekg(0, fposteriors.beg);
			//	//Adapting to read in chunks:
			//fposteriors_length = 137 * 45 * sizeof(float);
			//printf("char: %u, float: %u\n", sizeof(char), sizeof(float));
			//char *buffer_posteriors = new char[fposteriors_length];
			//fposteriors.read(buffer_posteriors, fposteriors_length * sizeof(char));
			//fposteriors.close();
			//float *posteriors = (float*)buffer_posteriors;
			//Memory Mapped File:

			file_mapping m_filePosteriors(posteriorsFile.data(), read_only);
			mapped_region regionPosteriors(m_filePosteriors, read_only);
			void * addrPosteriors = regionPosteriors.get_address();
			std::size_t sizePosteriors = regionPosteriors.get_size();
			float *posteriors = (float*)addrPosteriors;

			//Posteriors Offsets
			fposteriors_offsets.open(posteriorsOffsetsFile.data(), std::ios::binary | std::ios::in);
			fposteriors_offsets.seekg(0, fposteriors_offsets.end);
			int fposteriors_offsets_length = fposteriors_offsets.tellg();
			fposteriors_offsets.seekg(0, fposteriors_offsets.beg);
			char *buffer_fposteriors_offsets = new char[fposteriors_offsets_length];
			fposteriors_offsets.read(buffer_fposteriors_offsets, fposteriors_offsets_length * sizeof(char));
			fposteriors_offsets.close();
			uint32_t *posteriors_offsets = (uint32_t*)buffer_fposteriors_offsets;

			//Posteriors Lengths
			fposteriors_lengths.open(posteriorsLengthsFile.data(), std::ios::binary | std::ios::in);
			fposteriors_lengths.seekg(0, fposteriors_lengths.end);
			int fposteriors_lengths_length = fposteriors_lengths.tellg();
			fposteriors_lengths.seekg(0, fposteriors_lengths.beg);
			char *buffer_fposteriors_lengths = new char[fposteriors_lengths_length];
			fposteriors_lengths.read(buffer_fposteriors_lengths, fposteriors_lengths_length * sizeof(char));
			fposteriors_lengths.close();
			uint16_t *posteriors_lengths = (uint16_t*)buffer_fposteriors_lengths;

			//Dictionary
			//fdictionary.open(dictionaryFile.data(), std::ios::binary | std::ios::in);
			//fdictionary.seekg(0, fdictionary.end);
			//int fdictionary_length = fdictionary.tellg();
			//fdictionary.seekg(0, fdictionary.beg);
			//char *buffer_fdictionary = new char[fdictionary_length];
			//fdictionary.read(buffer_fdictionary, fdictionary_length * sizeof(char));
			//fdictionary.close();
			//int *dictionary = (int*)buffer_fdictionary;
			//Using Memory Mapped Files:
			file_mapping m_fileDictionary(dictionaryFile.data(), read_only);
			mapped_region regionDictionary(m_fileDictionary, read_only);
			void * addrDictionary = regionDictionary.get_address();
			std::size_t sizeDictionary = regionDictionary.get_size();
			int8_t *dictionary = (int8_t*)addrDictionary;
			//int *dictionary = (int*)addrDictionary;

			//Dictionary Offfsets
			fdictionary_offsets.open(dictionaryOffsetsFile.data(), std::ios::binary | std::ios::in);
			fdictionary_offsets.seekg(0, fdictionary_offsets.end);
			int fdictionary_offsets_length = fdictionary_offsets.tellg();
			fdictionary_offsets.seekg(0, fdictionary_offsets.beg);
			char *buffer_fdictionary_offsets = new char[fdictionary_offsets_length];
			fdictionary_offsets.read(buffer_fdictionary_offsets, fdictionary_offsets_length * sizeof(char));
			fdictionary_offsets.close();
			int *dictionary_offsets = (int*)buffer_fdictionary_offsets;

			//Dictionary Lengths
			fdictionary_lengths.open(dictionaryLengthsFile.data(), std::ios::binary | std::ios::in);
			fdictionary_lengths.seekg(0, fdictionary_lengths.end);
			int fdictionary_lengths_length = fdictionary_lengths.tellg();
			fdictionary_lengths.seekg(0, fdictionary_lengths.beg);
			char *buffer_fdictionary_lengths = new char[fdictionary_lengths_length];
			fdictionary_lengths.read(buffer_fdictionary_lengths, fdictionary_lengths_length * sizeof(char));
			fdictionary_lengths.close();
			int *dictionary_lengths = (int*)buffer_fdictionary_lengths;


			int out = 45;
			int T;
			int n_attempts=fposteriors_lengths_length / sizeof(uint16_t);
			printf("num attempts %i\n", n_attempts);

			//Iterate on posteriors
			
			float *current_pp = posteriors;
			uint32_t *current_pp_offset = posteriors_offsets+1;
			uint16_t *current_pp_lengths = posteriors_lengths;

			current_pp = posteriors + (*current_pp_offset)*out;

			T = (int)*(current_pp_lengths+1);
			printf("first attempt T: %i\n", T);

			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> pp(current_pp, T, out);

			printf("\n pp: (%i,%i)\n", pp.rows(), pp.cols());
			printf("pp(0,0):%f\n", pp(0, 0));

			//Text precision:
			Eigen::IOFormat precision(9);
			std::ofstream file;
			file.open("tests/test-data/check_1st_pp.txt", std::ios::out);
			file << pp.format(precision) << "\n";
			file.close();


			//Iterate on dictionary

			int8_t *current_w = dictionary;
			//int *current_w = dictionary;
			int *current_w_offset = dictionary_offsets;
			int *current_w_lengths = dictionary_lengths;

			int n_entries = fdictionary_lengths_length / sizeof(int);
			printf(" \ndictionary entries: %u\n", n_entries);


			float tscore;
			float best = -100000000.0f;
			int8_t *best_O= dictionary;
			int best_L= *(current_w_lengths);

			struct best_t {
				float tscore= -100000000.0f;;
				int8_t *offset;// = NULL;
				int length=0;
			};

			int N = 1032;
			struct best_t bestN[1032];

		

			for (int e = 0; e < n_entries; e++) {
				current_w = dictionary + *(current_w_offset+e);
				int L = *(current_w_lengths + e);


				//printf("first word L: %u\n", L);

				Eigen::Map<Eigen::Matrix<int8_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> w(current_w, L, 1);
				//Eigen::Map<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> w(current_w, L, 1);

				//printf("\n w: (%i,%i)\n", w.rows(), w.cols());
				//printf("w(0,0):%u    *int: %d\n", w(0, 0), *current_w);

				//Text precision:
				//Eigen::IOFormat precision(1);
				//std::ofstream file;
				//file.open("tests/test-data/check_1st_w.txt", std::ios::out);
				//file << w.cast<short>() << "\n";
				//file.close();

				tscore = Align(w, L, pp, T);
				if (tscore > best) {
					best = tscore;
					best_O = dictionary + *(current_w_offset + e);
					best_L = *(current_w_lengths + e);
				}
					
				//printf("tscore: %f", tscore);

				//BEST N(10)
				int j = 0;
				while (j < N && tscore < bestN[j].tscore) {
					j++;
				}
				if ( j <= (N - 2) ) {
					if (tscore > bestN[j].tscore) {
						for (int k = N-1; k > j; k--) {
							bestN[k].tscore = bestN[k - 1].tscore;
							bestN[k].offset = bestN[k - 1].offset;
							bestN[k].length = bestN[k - 1].length;
						}
						bestN[j].tscore = tscore;
						bestN[j].offset = dictionary + *(current_w_offset + e);
						bestN[j].length = *(current_w_lengths + e);
					}
				}
				else {
						if (tscore > bestN[j].tscore) {
							bestN[j].tscore = tscore;
							bestN[j].offset = dictionary + *(current_w_offset + e);
							bestN[j].length = *(current_w_lengths + e);
						}

				}//BEST N
				//printf("%d, %f, %i\n", e, tscore, j);
				//printf("[");
				//for (int i = 0; i < N; i++)
				//	printf(" %f", bestN[i].tscore);
				//printf("]\n");
			}// dictionary search


			file.open("tests/test-data/check_bestN.txt", std::ios::out);

			//printf("\n");
			//for (int i = 0; i < N; i++) {
			//	printf("%i\t%f\t", i, bestN[i].tscore);
			//	int L = bestN[i].length;
			//	//current_w = bestN[i].offset;
			//	Eigen::Map<Eigen::Matrix<int8_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> w(bestN[i].offset, L, 1);
			//	//printf("%hi", w(2, 0));
			//	for (int j = 0; j < L; j++) {
			//		printf("%hi ", w(j, 0));
			//	}
			//	printf("\n");
			//}

			for (int i = 0; i < N; i++) {
				//printf("%i\t%f\t", i, bestN[i].tscore);
				file << i << "\t" << bestN[i].tscore << "\t";
				int L = bestN[i].length;
				//current_w = bestN[i].offset;
				Eigen::Map<Eigen::Matrix<int8_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> w(bestN[i].offset, L, 1);
				//printf("%hi", w(2, 0));
				for (int j = 0; j < L; j++) {
					//printf("%hi ", w(j, 0));
					file << phonemes[w(j, 0)] << " ";
				}
				//printf("\n");
				file << "\n";
			}

			file.close();

			printf("best score: %f", best);
			Eigen::Map<Eigen::Matrix<int8_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> best_w(best_O, best_L, 1);
			//std::ofstream file;
			file.open("tests/test-data/check_best_w.txt", std::ios::out);
			file << best_w.cast<short>() << "\n";
			file.close();
			

			//delete[] buffer_posteriors;
			delete[] buffer_fposteriors_offsets;
			delete[] buffer_fposteriors_lengths;
			//delete[] buffer_fdictionary;
			delete[] buffer_fdictionary_offsets;
			delete[] buffer_fdictionary_lengths;
		}

		//******************************************************************************************
		//**** Main algorithm *******************

		//Spaces coding:
		enum space_t {semantic=0, phonemic};

		//Structure to output results
		struct sltv1_t {
			float tscore_target;
			float tscore_sem;
			float tscore_phon;
			uint8_t transcript_target[256];
			uint8_t transcript_sem[256];
			uint8_t transcript_phon[256];
		};

		void flush_to_file_sltv1(std::ofstream& file, sltv1_t *features, int num_attempts) {
			for (int i = 0; i < num_attempts; i++) {
				file.write((char*)&features[i].tscore_target, sizeof(float));
				file.write((char*)&(features[i].transcript_target[0]), sizeof(uint8_t) * 256);
				file.write((char*)&features[i].tscore_sem, sizeof(float));
				file.write((char*)&(features[i].transcript_sem[0]), sizeof(uint8_t) * 256);
				file.write((char*)&features[i].tscore_phon, sizeof(float));
				file.write((char*)&(features[i].transcript_phon[0]), sizeof(uint8_t) * 256);
			}
		}
		

		//Vocabulary of words
		std::vector<std::string> target_words({
			"ant",
			"arch",
			"arm",
			"axe",
			"back",
			"badge",
			"bag",
			"ball",
			"band",
			"bat",
			"bath",
			"beach",
			"beak",
			"bean",
			"bear",
			"beard",
			"bed",
			"bee",
			"beer",
			"bell",
			"belt",
			"bench",
			"bib",
			"bike",
			"bin",
			"bird",
			"blind",
			"block",
			"blouse",
			"boar",
			"board",
			"boat",
			"bomb",
			"bone",
			"book",
			"boot",
			"bow",
			"bowl",
			"box",
			"boy",
			"bra",
			"brain",
			"branch",
			"bread",
			"bricks",
			"bride",
			"bridge",
			"broom",
			"brush",
			"bug",
			"bulb",
			"bull",
			"bun",
			"buoy",
			"bus",
			"cage",
			"cake",
			"can",
			"cane",
			"cap",
			"cape",
			"car",
			"card",
			"cart",
			"case",
			"cat",
			"cave",
			"chain",
			"chair",
			"cheese",
			"chef",
			"chick",
			"chips",
			"church",
			"clamp",
			"claw",
			"clock",
			"cloud",
			"clown",
			"club",
			"coach",
			"coal",
			"coat",
			"coins",
			"comb",
			"cone",
			"core",
			"cork",
			"corn",
			"cot",
			"court",
			"cow",
			"crab",
			"crane",
			"crisps",
			"crook",
			"cross",
			"crown",
			"crutch",
			"cup",
			"dart",
			"deer",
			"desk",
			"dice",
			"disc",
			"dog",
			"doll",
			"dome",
			"door",
			"dove",
			"drain",
			"drawers",
			"dress",
			"drill",
			"drum",
			"duck",
			"ear",
			"egg",
			"eye",
			"face",
			"fan",
			"farm",
			"fence",
			"fig",
			"film",
			"fire",
			"fish",
			"fist",
			"flag",
			"flute",
			"fly",
			"foot",
			"fork",
			"fox",
			"frame",
			"fridge",
			"frog",
			"fruit",
			"gas",
			"gate",
			"ghost",
			"gift",
			"girl",
			"glass",
			"globe",
			"glove",
			"glue",
			"gnome",
			"goat",
			"golf",
			"grapes",
			"grass",
			"grave",
			"groom",
			"gun",
			"hair",
			"ham",
			"hand",
			"harp",
			"hat",
			"hawk",
			"hay",
			"head",
			"heart",
			"hedge",
			"heel",
			"hinge",
			"hive",
			"hoe",
			"hoof",
			"hook",
			"horn",
			"horse",
			"hose",
			"house",
			"hutch",
			"ice",
			"ink",
			"iron",
			"jail",
			"jar",
			"jeans",
			"jeep",
			"judge",
			"jug",
			"key",
			"kilt",
			"king",
			"kite",
			"knee",
			"knife",
			"knight",
			"knob",
			"knot",
			"lamb",
			"lamp",
			"leaf",
			"leek",
			"leg",
			"lens",
			"lights",
			"lion",
			"lips",
			"lock",
			"log",
			"lungs",
			"maid",
			"man",
			"map",
			"mask",
			"match",
			"maze",
			"mill",
			"moon",
			"moose",
			"mop",
			"mouse",
			"mouth",
			"mug",
			"nail",
			"neck",
			"nest",
			"net",
			"noose",
			"nose",
			"note",
			"nun",
			"nurse",
			"nut",
			"oar",
			"owl",
			"paint",
			"pan",
			"park",
			"paw",
			"peach",
			"pear",
			"peas",
			"peg",
			"pen",
			"phone",
			"pie",
			"pig",
			"pills",
			"pin",
			"pipe",
			"plane",
			"plant",
			"plate",
			"pliers",
			"plug",
			"pool",
			"pope",
			"pot",
			"pound",
			"pram",
			"priest",
			"purse",
			"queen",
			"raft",
			"rain",
			"rake",
			"ramp",
			"rat",
			"ring",
			"road",
			"rock",
			"rod",
			"roll",
			"roof",
			"root",
			"rope",
			"rose",
			"rug",
			"sack",
			"safe",
			"sail",
			"salt",
			"sand",
			"saw",
			"scale",
			"scarf",
			"school",
			"scoop",
			"screen",
			"screw",
			"scroll",
			"seal",
			"seat",
			"shark",
			"shawl",
			"shed",
			"sheep",
			"shell",
			"shelves",
			"shield",
			"ship",
			"shirt",
			"shoe",
			"shop",
			"shorts",
			"sieve",
			"sign",
			"sink",
			"skate",
			"skirt",
			"skis",
			"skull",
			"sled",
			"sleeve",
			"sleigh",
			"slide",
			"smoke",
			"snail",
			"snake",
			"snow",
			"sock",
			"sole",
			"soup",
			"space",
			"spade",
			"sponge",
			"spoon",
			"spray",
			"spring",
			"square",
			"stairs",
			"stalk",
			"stamp",
			"stand",
			"star",
			"steam",
			"step",
			"stick",
			"stool",
			"stork",
			"stove",
			"straw",
			"stump",
			"suit",
			"sun",
			"swan",
			"sweet",
			"swim",
			"swing",
			"switch",
			"sword",
			"tail",
			"tank",
			"tap",
			"tape",
			"tear",
			"teeth",
			"tent",
			"thread",
			"thumb",
			"tie",
			"tile",
			"till",
			"tin",
			"tire",
			"toad",
			"toast",
			"toe",
			"tongue",
			"tools",
			"top",
			"torch",
			"towel",
			"toys",
			"train",
			"tray",
			"tree",
			"truck",
			"tube",
			"van",
			"vase",
			"vest",
			"vet",
			"wall",
			"wand",
			"watch",
			"wave",
			"well",
			"whale",
			"wheat",
			"wheel",
			"whip",
			"whisk",
			"wig",
			"wind",
			"wine",
			"wing",
			"witch",
			"wolf",
			"wood",
			"wool",
			"worm",
			"wreath",
				"zip" });

		template <typename DerivedPo>
		void process_space(space_t space, const Eigen::MatrixBase<DerivedPo>& pp, int ppl, std::string target_word, sltv1_t *features, int attempt, uint8_t *zeros) {
			//std::string target_word("whale");
			//printf("sems space target word: %s", target_word);
			std::string space_type;
			switch (space) {
				case semantic:
					space_type="semantic";
					break;
				case phonemic:
					space_type="phonemic";
					break;
			}
			std::string dir_dictionaries("tests/test-data/dictionary_files/dictionaries/");
			std::string dir_offsets("tests/test-data/dictionary_files/offsets/");
			std::string dir_lengths("tests/test-data/dictionary_files/lengths/");
			std::string dictionaryFile = dir_dictionaries + space_type + "_space_" + target_word + ".bin";
			   //std::string offsetsFile("tests/test-data/sem_space_offsets_whale.bin");
			std::string offsetsFile(dir_offsets + space_type + "_space_offsets_" + target_word + ".bin");
			std::string lengthsFile = dir_lengths + space_type + "_space_lengths_" + target_word + ".bin";
			//printf("%s\n%s\n%s\n", dictionaryFile.data(), offsetsFile.data(), lengthsFile.data());

			//Map files to memory
			//using namespace boost::interprocess;


			file_mapping m_dictionaryFile(dictionaryFile.data(), read_only);
			//file_mapping m_SEMdictionaryFile("tests/test-data/dictionary_files/dictionaries/dictionary_400words_phonemes.bin", read_only);
			mapped_region m_regionDictionaryFile(m_dictionaryFile, read_only);
			void * addrDictionary = m_regionDictionaryFile.get_address();
			uint8_t *dictionary = (uint8_t*)addrDictionary;

			file_mapping m_offsetsFile(offsetsFile.data(), read_only);
			mapped_region m_regionOffsetsFile(m_offsetsFile, read_only);
			void * addrOffsets = m_regionOffsetsFile.get_address();
			int *offsets = (int*)addrOffsets;

			file_mapping m_lengthsFile(lengthsFile.data(), read_only);
			mapped_region m_regionLengthsFile(m_lengthsFile, read_only);
			void * addrLengths = m_regionLengthsFile.get_address();
			std::size_t dictionary_entries = m_regionLengthsFile.get_size() / sizeof(int);
			int *lengths = (int*)addrLengths;




			//Using streams
			/*std::ifstream infile;
			infile.open(dictionaryFile.data(), std::ios::binary | std::ios::in);
			infile.seekg(0, infile.end);
			int infile_length = infile.tellg();
			infile.seekg(0, infile.beg);
			char *buffer = new char[infile_length];
			infile.read(buffer, infile_length);
			infile.close();
			uint8_t *dictionary = (uint8_t*)buffer;

			std::ifstream infileO;
			infileO.open(offsetsFile.data(), std::ios::binary | std::ios::in);
			infileO.seekg(0, infileO.end);
			int infile_lengthO = infileO.tellg();
			infileO.seekg(0, infileO.beg);
			char *bufferO = new char[infile_lengthO];
			infileO.read(bufferO, infile_lengthO);
			infileO.close();
			int *offsets = (int*)bufferO;

			std::ifstream infileL;
			infileL.open(lengthsFile.data(), std::ios::binary | std::ios::in);
			infileL.seekg(0, infileL.end);
			int infile_lengthL = infileL.tellg();
			infileL.seekg(0, infileL.beg);
			char *bufferL = new char[infile_lengthL];
			infileL.read(bufferL, infile_lengthL);
			infileL.close();
			int *lengths = (int*)bufferL;
			int dictionary_entries = infile_lengthL / sizeof(int);*/





			float tscore = -1000000;
			float best_tscore = -1000000;
			int best_entry = 0;
			//dictionary_entries = 1;
			for (int i = 0; i < dictionary_entries; i++) {
				Eigen::Map<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> sem_phn(dictionary + *(offsets + i), *(lengths + i), 1);
				tscore = Align(sem_phn, *(lengths + i), pp, ppl);
				if (tscore > best_tscore) {
					best_tscore = tscore;
					best_entry = i;
				}
			}
			printf("%s space -->   best tscore: %.9f, best entry: %d\n", space_type.data(), best_tscore, best_entry);

			switch(space) {
				case semantic:
					features[attempt].tscore_sem = best_tscore;
					std::memcpy(&features[attempt].transcript_sem[0], &zeros[0], 256 * sizeof(uint8_t));
					std::memcpy(&features[attempt].transcript_sem[0], (dictionary + *(offsets + best_entry)), (*(lengths + best_entry)) * sizeof(uint8_t));
					break;
				case phonemic:
					features[attempt].tscore_phon = best_tscore;
					std::memcpy(&features[attempt].transcript_phon[0], &zeros[0], 256 * sizeof(uint8_t));
					std::memcpy(&features[attempt].transcript_phon[0], (dictionary + *(offsets + best_entry)), (*(lengths + best_entry)) * sizeof(uint8_t));
					break;
			}
			//delete[] buffer;
			//delete[] bufferO;
			//delete[] bufferL;
		
			//if (space == "semantic") {
			//	features[attempt].tscore_sem = best_tscore;
			//	std::memcpy(&features[attempt].transcript_sem[0], &zeros[0], 256 * sizeof(uint8_t));
			//	std::memcpy(&features[attempt].transcript_sem[0], (dictionary + *(offsets + best_entry)), (*(lengths + best_entry)) * sizeof(uint8_t));
			//} else {
			//	features[attempt].tscore_phon = best_tscore;
			//	std::memcpy(&features[attempt].transcript_phon[0], &zeros[0], 256 * sizeof(uint8_t));
			//	std::memcpy(&features[attempt].transcript_phon[0], (dictionary + *(offsets + best_entry)), (*(lengths + best_entry)) * sizeof(uint8_t));
			//}
			
		}

		void process_attempt(float *ppfloat, int ppl, uint8_t *tphn, int tphn_length, sltv1_t *features, int attempt, int out, std::string target_word) {
			//printf("\n%s\n", target_word.data());

			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> pp(ppfloat, ppl, out);
			//Just for testing everything is allright
			//printf("\n pp: (%i,%i)\n", pp.rows(), pp.cols());
			//printf("pp(0,0):%f\n", pp(0, 0));

			////Text precision:
			//Eigen::IOFormat precision(9);
			//std::ofstream file;
			//file.open("tests/test-data/check_1_pp.txt", std::ios::out);
			//file << pp.format(precision) << "\n";
			//file.close();


			Eigen::Map<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> phonemes(tphn, tphn_length, 1);
			//Just for testing everything is allright
			//std::ofstream filew;
			//filew.open("tests/test-data/check_1_target_phonemes.txt", std::ios::out);
			//filew << phonemes.cast<short>() << "\n";
			//filew.close();

			//initilize
			uint8_t zeros[256];
			for (int i = 0; i < 256; i++)
				zeros[i] = 0;

			float tscore = -1000000;
			float best_tscore = -1000000;
			int best_entry = 0;
	
			//Target Space (one single word: the target word)
			tscore = Align(phonemes, tphn_length, pp, ppl);
			features[attempt].tscore_target = tscore;
			std::memcpy(&features[attempt].transcript_target[0], &zeros[0], 256 * sizeof(uint8_t));
			std::memcpy(&features[attempt].transcript_target[0], &tphn[0], tphn_length * sizeof(uint8_t));

			//Semantic Space
			process_space(semantic, pp, ppl, target_word, &features[0], attempt, &zeros[0]);
			//Phonemic Space
			process_space(phonemic, pp, ppl, target_word, &features[0], attempt, &zeros[0]);


			////Phonemic Space
			//dictionaryFile = dir_dictionaries + "phon_space_" + target_word + ".bin";
			//offsetsFile = dir_dictionaries + "phon_space_offsets_" + target_word + ".bin";
			//lengthsFile = dir_dictionaries + "phon_space_lengths_" + target_word + ".bin";
			////Map files to memory
			//using namespace boost::interprocess;

			//file_mapping m_pdictionaryFile(dictionaryFile.data(), read_only);
			//mapped_region m_pregionDictionaryFile(m_pdictionaryFile, read_only);
			//void * addrpDictionary = m_regionDictionaryFile.get_address();
			//uint8_t *pdictionary = (uint8_t*)addrpDictionary;

			//file_mapping m_poffsetsFile(offsetsFile.data(), read_only);
			//mapped_region m_pregionOffsetsFile(m_poffsetsFile, read_only);
			//void * addrpOffsets = m_pregionOffsetsFile.get_address();
			//int *poffsets = (int*)addrpOffsets;

			//file_mapping m_plengthsFile(lengthsFile.data(), read_only);
			//mapped_region m_pregionLengthsFile(m_lengthsFile, read_only);
			//void * addrpLengths = m_pregionLengthsFile.get_address();
			//std::size_t pdictionary_entries = m_pregionLengthsFile.get_size() / sizeof(int);
			//int *plengths = (int*)addrpLengths;

			//best_tscore = -1000000;
			//best_entry = 0;
			//pdictionary_entries = 1;
			//for (int i = 0; i < pdictionary_entries; i++) {
			//	Eigen::Map<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >> phon_phn(dictionary + *(poffsets + i), *(plengths + i), 1);
			//	tscore = Align(phon_phn, *(lengths + i), pp, ppl);
			//	if (tscore > best_tscore) {
			//		best_tscore = tscore;
			//		best_entry = i;
			//	}
			//}
			//features[attempt].tscore_phon = best_tscore;
			//std::memcpy(&features[attempt].transcript_phon[0], &zeros[0], 256 * sizeof(uint8_t));
			//std::memcpy(&features[attempt].transcript_phon[0], &tphn[0], tphn_length * sizeof(uint8_t));
		}

		void test_SLTfeatures_v1() {

		std::string posteriorsFile("tests/test-data/posteriors_files/posteriors/posteriors_best.bin");
		std::string offsetsFile("tests/test-data/posteriors_files/offsets/posteriors_best_offsets.bin");
		std::string lengthsFile("tests/test-data/posteriors_files/lengths/posteriors_best_lengths.bin");

		std::string word_codesFile("tests/test-data/dictionary_files/dictionaries/dictionary_400words_phonemes.bin");
		std::string word_codesOffsetsFile("tests/test-data/dictionary_files/offsets/dictionary_400words_phonemes_offsets.bin");
		std::string word_codesLengthsFile("tests/test-data/dictionary_files/lengths/dictionary_400words_phonemes_lengths.bin");

		std::string attempts_targetsFile("tests/test-data/dictionary_files/dictionaries/attempts_word_codes.bin");

		//using namespace boost::interprocess;

		//mapping all files into memory:

		file_mapping m_filePosteriors(posteriorsFile.data(), read_only);
		mapped_region m_regionPosteriors(m_filePosteriors, read_only);
		void * addrPosteriors = m_regionPosteriors.get_address();
		std::size_t sizePosteriors = m_regionPosteriors.get_size();
		float *posteriors = (float*)addrPosteriors;

		//int num_attempts = size_posteriors

		file_mapping m_fileOffsets(offsetsFile.data(), read_only);
		mapped_region m_regionOffsets(m_fileOffsets, read_only);
		void * addrOffsets = m_regionOffsets.get_address();
		//std::size_t sizePosteriors = regionPosteriors.get_size();
		int *offsets = (int*)addrOffsets;

		file_mapping m_fileLengths(lengthsFile.data(), read_only);
		mapped_region m_regionLengths(m_fileLengths, read_only);
		void * addrLengths = m_regionLengths.get_address();
		std::size_t num_attempts = m_regionLengths.get_size()/sizeof(int);
		int *lengths = (int*)addrLengths;

		file_mapping m_fileWordCodes(word_codesFile.data(), read_only);
		mapped_region m_regionWordCodes(m_fileWordCodes, read_only);
		void * addrWordCodes = m_regionWordCodes.get_address();
		//std::size_t sizePosteriors = regionPosteriors.get_size();
		uint8_t *wordcodes = (uint8_t*)addrWordCodes;  //These are the numbers that link tonthe phonetic transcript of a particular word, each word one number

		file_mapping m_fileWordCodesOffsets(word_codesOffsetsFile.data(), read_only);
		mapped_region m_regionWordCodesOffsets(m_fileWordCodesOffsets, read_only);
		void * addrWordCodesOffsets = m_regionWordCodesOffsets.get_address();
		//std::size_t sizePosteriors = regionPosteriors.get_size();
		int *wordcodesOffsets = (int*)addrWordCodesOffsets;

		file_mapping m_fileWordCodesLengths(word_codesLengthsFile.data(), read_only);
		mapped_region m_regionWordCodesLengths(m_fileWordCodesLengths, read_only);
		void * addrWordCodesLengths = m_regionWordCodesLengths.get_address();
		//std::size_t sizePosteriors = regionPosteriors.get_size();
		int *wordcodesLengths = (int*)addrWordCodesLengths;

		file_mapping m_fileAttemtpsTargets(attempts_targetsFile.data(), read_only);
		mapped_region m_regionAttemtpsTargets(m_fileAttemtpsTargets, read_only);
		void * addrAttemtpsTargets = m_regionAttemtpsTargets.get_address();
		//std::size_t sizePosteriors = regionPosteriors.get_size();
		int *attemptstargets = (int*)addrAttemtpsTargets;  //These are the phonetic transcripts of the targets words for each attempt
		

		printf("\n attempts: %i \n", num_attempts); //7200
		struct sltv1_t *SLTfeatures; 
		SLTfeatures=(sltv1_t *)calloc(num_attempts, sizeof(sltv1_t));
		


		int a = 1; //attempt idx
		int out = 45; //# ARPAbet phonemes+'sil'

		//iterating posteriors
		//initialisation
		float *current_pp = posteriors;
		//uint32_t *current_pp_offset = offsets;
		int *current_pp_length = lengths;

		//current_pp = posteriors + (*current_pp_offset)*out;
		current_pp =  posteriors + (*(offsets+a))*out;
		current_pp_length = lengths + a;

		//finding the target phonemes
		int current_target_code = *(attemptstargets + a);
		uint8_t *current_target_phonemes = wordcodes + *(wordcodesOffsets + current_target_code);
		int phonemes_length = *(wordcodesLengths + current_target_code);

		
		//testing for one attempt
//		process_attempt(current_pp,*current_pp_length,current_target_phonemes,phonemes_length,&SLTfeatures[0],a,out,target_words[target_code]);

		//printf("\ntarget tscore: %.9f\n", SLTfeatures[a].tscore_target);
		//for (int i = 0; i < 256; i++)
		//	printf(" %u", SLTfeatures[a].transcript_target[i]);

		//printf("\n sem tscore: %.9f\n", SLTfeatures[a].tscore_sem);
		//for (int i = 0; i < 256; i++)
		//	printf(" %u", SLTfeatures[a].transcript_sem[i]);

		//printf("\n phon tscore: %.9f\n", SLTfeatures[a].tscore_phon);
		//for (int i = 0; i < 256; i++)
		//	printf(" %u", SLTfeatures[a].transcript_phon[i]);

		high_resolution_clock::time_point start;
		high_resolution_clock::time_point end;
		//std::ofstream report;
		//report.open("tests/test-data/report-SLT-Timings-2019-07-12.txt", std::ios::out);
		//std::vector<double> timings;
		////main loop. This can be multithreaded
		start = high_resolution_clock::now();
		//for (int a = 0; a < num_attempts; a++) {
		for (int a = 0; a < 1; a++) {
			//start = high_resolution_clock::now();
			//Iterating posteriors
			current_pp = posteriors + (*(offsets + a))*out;
			current_pp_length = lengths + a;
			//Iterating target words
			current_target_code = *(attemptstargets + a);
			current_target_phonemes = wordcodes + *(wordcodesOffsets + current_target_code);
			phonemes_length = *(wordcodesLengths + current_target_code);
			printf("%s\n", target_words[current_target_code].data());
			//Process the naming attempt
			process_attempt(current_pp, *current_pp_length, current_target_phonemes, phonemes_length, &SLTfeatures[0], a, out, target_words[current_target_code]);
			//end = high_resolution_clock::now();
			printf("%d:\t", a);
			//auto dur_us = duration<double, std::micro>(end - start).count();
			//auto dur_ms = duration<double, std::milli>(end - start).count();
			//printf("%d: --- Time: %lfus %lfms\n", a, dur_us, dur_ms);
			//timings.push_back(dur_ms);
			//report << dur_ms << "," << std::endl;
			
			//printf("\ntarget tscore: %.9f\n", SLTfeatures[a].tscore_target);
			//for (int i = 0; i < 256; i++)
			//	printf(" %u", SLTfeatures[a].transcript_target[i]);

			//printf("\n sem tscore: %.9f\n", SLTfeatures[a].tscore_sem);
			//for (int i = 0; i < 256; i++)
			//	printf(" %u", SLTfeatures[a].transcript_sem[i]);

			//printf("\n phon tscore: %.9f\n", SLTfeatures[a].tscore_phon);
			//for (int i = 0; i < 256; i++)
			//	printf(" %u", SLTfeatures[a].transcript_phon[i]);
		}
		end = high_resolution_clock::now();
		auto dur_us = duration<double, std::micro>(end - start).count();
		auto dur_ms = duration<double, std::milli>(end - start).count();
		printf("\n\n --- SLT extraction Time: %lfus %lfms\n", dur_us, dur_ms);

		//report.close();

		//Storing results in a binary file
		start = high_resolution_clock::now();
		std::ofstream outfile;
		std::string binfile("tests/test-data/SLT_v1_results_test.bin");
		outfile = std::ofstream(binfile.data(), std::ios::binary | std::ios::out);

		flush_to_file_sltv1(outfile, &SLTfeatures[0], num_attempts);

		outfile.close();
		
		end = high_resolution_clock::now();
		dur_us = duration<double, std::micro>(end - start).count();
		dur_ms = duration<double, std::milli>(end - start).count();
		printf("\n\n --- SLT storage Time: %lfus %lfms\n", dur_us, dur_ms);



		free(SLTfeatures);

		}//SLT v1

	}//metrics
}//nni