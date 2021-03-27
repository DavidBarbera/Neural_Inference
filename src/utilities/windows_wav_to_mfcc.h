
#include <windows.h>


#define MAX_PATH 260

using namespace std::chrono;




namespace nni {

	namespace utilities {

		std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;

		std::string ExePath() {
			char charbuffer[MAX_PATH];
			wchar_t  wbuffer[MAX_PATH];
			//std::wstring wbuffer;
			GetModuleFileName(NULL, wbuffer, MAX_PATH);
			size_t ret;
			wcstombs_s(&ret, charbuffer,sizeof(charbuffer), wbuffer, sizeof(wbuffer));
			if (ret == 32) charbuffer[MAX_PATH-1] = '\0';
			std::string buffer = converter.to_bytes(wbuffer);
			std::string::size_type pos = std::string(buffer).find_last_of("\\/");
			return std::string(buffer).substr(0, pos);
		}

		int create_mfcc_folder(std::wstring mfcc_folder) {

			if (CreateDirectory(mfcc_folder.c_str(), NULL) ||
				ERROR_ALREADY_EXISTS == GetLastError())
			{
				//// CopyFile(...)
				//CopyFile(Input.c_str(), string(OutputFolder + CopiedFile).c_str(), TRUE);
				printf("mfcc folder created...\n\n");
				return 0;
			}
			else
			{
				printf("Failed to create mfcc directory for the output.\n\n");
				// Failed to create directory.
				return -1;
			}

			
		}

		int create_pp_folder(std::wstring pp_folder) {

			if (CreateDirectory(pp_folder.c_str(), NULL) ||
				ERROR_ALREADY_EXISTS == GetLastError())
			{
				//// CopyFile(...)
				//CopyFile(Input.c_str(), string(OutputFolder + CopiedFile).c_str(), TRUE);
				printf("pp folder created...\n\n");
				return 0;
			}
			else
			{
				printf("Failed to create pp directory for the output.\n\n");
				// Failed to create directory.
				return -1;
			}


		}



		

		void wavs_to_mfccs_files() {

			high_resolution_clock::time_point start = high_resolution_clock::now();

			WIN32_FIND_DATA data;
			// DIRECTORY

			//std::string wav_path = ExePath() + "\\..\\..\\tests\\wavs";
			//std::string list_files = wav_path + "\\*";
			//std::string mfcc_path = wav_path + "\\..\\mfcc";

			std::string wav_path = ExePath() + "\\wavs";
			std::wstring list_files = converter.from_bytes(wav_path) + L"\\*";
			std::wstring mfcc_path = converter.from_bytes(wav_path) + L"\\..\\mfccs";

			HANDLE hFind = FindFirstFile(list_files.data(), &data);


			create_mfcc_folder(mfcc_path);

			if (hFind != INVALID_HANDLE_VALUE) {
				do {
					//std::cout << data.cFileName << std::endl;
					std::wstring cf(data.cFileName);
					std::string current_file = converter.to_bytes(cf);
					//printf("file: %s\twav extension?: %s\n", current_file.data(), current_file.substr(current_file.find_last_of(".")+1)=="wav"?"YES":"NO");
					//printf("%s\n, extension: %s", data.cFileName, PathFindExtensionA(data).data());
					std::string name_only = current_file.substr(0, current_file.find_last_of("."));
					if (current_file.substr(current_file.find_last_of(".") + 1) == "wav") {
						//printf("wav file path:\n%s\nmfcc file path:\n%s\n", (wav_path+"\\"+current_file).data(),(mfcc_path+"\\").data());
						//printf("Just the name: %s\n", name_only.data());

						wav_to_mfcc((wav_path + "\\" + current_file).data(), (converter.to_bytes(mfcc_path) + "\\" + name_only + ".mfcc").data());


					}


				} while (FindNextFile(hFind, &data));
				FindClose(hFind);
			}
			high_resolution_clock::time_point end = high_resolution_clock::now();
			auto dur_ms = duration<double, std::milli>(end - start).count();
			auto dur_us = duration<double, std::micro>(end - start).count();
			//printf("Total time:\n %.2d ms\n %  \u03BCs \n", dur_ms, dur_us);
			std::cout << dur_ms << " ms\n";
			std::cout << dur_us << " micro seconds\n";
		}

		void wavs_to_pp_files() {
			high_resolution_clock::time_point start = high_resolution_clock::now();

			WIN32_FIND_DATA data;
			// DIRECTORY

			//std::string wav_path = ExePath() + "\\..\\..\\tests\\wavs";
			//std::string list_files = wav_path + "\\*";
			//std::string mfcc_path = wav_path + "\\..\\mfcc";

			std::string path = ExePath();
			std::string wav_path = path + "\\wavs";
			std::wstring list_files = converter.from_bytes(wav_path) + L"\\*";
			std::wstring pp_path = converter.from_bytes(path) + L"\\pp";
			std::string model = path + "\\model_weights.bin";

			char *ppv;

			HANDLE hFind = FindFirstFile(list_files.data(), &data);


			create_pp_folder(pp_path);

			if (hFind != INVALID_HANDLE_VALUE) {
				do {
					//std::cout << data.cFileName << std::endl;
					std::wstring cf(data.cFileName);
					std::string current_file = converter.to_bytes(cf);
					//printf("file: %s\twav extension?: %s\n", current_file.data(), current_file.substr(current_file.find_last_of(".")+1)=="wav"?"YES":"NO");
					//printf("%s\n, extension: %s", data.cFileName, PathFindExtensionA(data).data());
					std::string name_only = current_file.substr(0, current_file.find_last_of("."));
					if (current_file.substr(current_file.find_last_of(".") + 1) == "wav") {
						//printf("wav file path:\n%s\nmfcc file path:\n%s\n", (wav_path+"\\"+current_file).data(),(mfcc_path+"\\").data());
						//printf("Just the name: %s\n", name_only.data());

						printf("\n%s", current_file.data());
						auto pp=wav_to_pp((wav_path + "\\" + current_file).data(), model.data() );
						ppv = (char*)pp.data();


						std::ofstream outfile;
						outfile.open((converter.to_bytes(pp_path) + "\\" + name_only + ".pp").data(), std::ios::binary | std::ios::out);
						outfile.write(ppv, pp.size() * sizeof(float));
						outfile.close();

					}


				} while (FindNextFile(hFind, &data));
				FindClose(hFind);
			}
			high_resolution_clock::time_point end = high_resolution_clock::now();
			auto dur_ms = duration<double, std::milli>(end - start).count();
			auto dur_us = duration<double, std::micro>(end - start).count();
			//printf("Total time:\n %.2d ms\n %  \u03BCs \n", dur_ms, dur_us);
			std::cout << "\n";
			std::cout << dur_ms << " ms\n";
			std::cout << dur_us << " micro seconds\n";
		}

		void wavs_to_pp_files_resample_norm() {
			high_resolution_clock::time_point start = high_resolution_clock::now();

			WIN32_FIND_DATA data;
			// DIRECTORY

			//std::string wav_path = ExePath() + "\\..\\..\\tests\\wavs";
			//std::string list_files = wav_path + "\\*";
			//std::string mfcc_path = wav_path + "\\..\\mfcc";

			std::string path = ExePath();
			std::string wav_path = path + "\\wavs";
			std::wstring list_files = converter.from_bytes(wav_path) + L"\\*";
			std::wstring pp_path = converter.from_bytes(path) + L"\\pp";
			std::wstring pp_res_norm = converter.from_bytes(path) + L"\\pp\\res_norm";
			std::wstring pp_res_notnorm = converter.from_bytes(path) + L"\\pp\\res_notnorm";
			std::wstring pp_notres_norm = converter.from_bytes(path) + L"\\pp\\notres_norm";
			std::wstring pp_notres_notnorm = converter.from_bytes(path) + L"\\pp\\notres_notnorm";
			std::string model = path + "\\model_weights.bin";

			char* ppv;

			HANDLE hFind = FindFirstFile(list_files.data(), &data);


			create_pp_folder(pp_path);
			create_pp_folder(pp_res_norm);
			create_pp_folder(pp_res_notnorm);
			create_pp_folder(pp_notres_norm);
			create_pp_folder(pp_notres_notnorm);

			if (hFind != INVALID_HANDLE_VALUE) {
				do {
					//std::cout << data.cFileName << std::endl;
					std::wstring cf(data.cFileName);
					std::string current_file = converter.to_bytes(cf);
					//printf("file: %s\twav extension?: %s\n", current_file.data(), current_file.substr(current_file.find_last_of(".")+1)=="wav"?"YES":"NO");
					//printf("%s\n, extension: %s", data.cFileName, PathFindExtensionA(data).data());
					std::string name_only = current_file.substr(0, current_file.find_last_of("."));
					if (current_file.substr(current_file.find_last_of(".") + 1) == "wav") {
						//printf("wav file path:\n%s\nmfcc file path:\n%s\n", (wav_path+"\\"+current_file).data(),(mfcc_path+"\\").data());
						//printf("Just the name: %s\n", name_only.data());

						//Resampled and norm
						auto pp = wav_to_pp((wav_path + "\\" + current_file).data(), model.data(),true,true);
						ppv = (char*)pp.data();

						std::ofstream outfile;
						outfile.open((converter.to_bytes(pp_res_norm) + "\\" + name_only + "_r_n.pp").data(), std::ios::binary | std::ios::out);
						outfile.write(ppv, pp.size() * sizeof(float));
						outfile.close();

						//Resampled and notnorm
						pp = wav_to_pp((wav_path + "\\" + current_file).data(), model.data(), true, false);
						ppv = (char*)pp.data();

						  //std::ofstream outfile;
						outfile.open((converter.to_bytes(pp_res_notnorm) + "\\" + name_only + "_r_nn.pp").data(), std::ios::binary | std::ios::out);
						outfile.write(ppv, pp.size() * sizeof(float));
						outfile.close();

						//Not resampled and norm
						pp = wav_to_pp((wav_path + "\\" + current_file).data(), model.data(), false, true);
						ppv = (char*)pp.data();

						  //std::ofstream outfile;
						outfile.open((converter.to_bytes(pp_notres_norm) + "\\" + name_only + "_nr_n.pp").data(), std::ios::binary | std::ios::out);
						outfile.write(ppv, pp.size() * sizeof(float));
						outfile.close();

						//Not resampled and not norm
						pp = wav_to_pp((wav_path + "\\" + current_file).data(), model.data(), false, false);
						ppv = (char*)pp.data();

						  //std::ofstream outfile;
						outfile.open((converter.to_bytes(pp_notres_notnorm) + "\\" + name_only + "_nr_nn.pp").data(), std::ios::binary | std::ios::out);
						outfile.write(ppv, pp.size() * sizeof(float));
						outfile.close();
					}


				} while (FindNextFile(hFind, &data));
				FindClose(hFind);
			}
			high_resolution_clock::time_point end = high_resolution_clock::now();
			auto dur_ms = duration<double, std::milli>(end - start).count();
			auto dur_us = duration<double, std::micro>(end - start).count();
			//printf("Total time:\n %.2d ms\n %  \u03BCs \n", dur_ms, dur_us);
			std::cout << dur_ms << " ms\n";
			std::cout << dur_us << " micro seconds\n";
		}


		void dtw_at_scale() {

			float dtw_dist = 2;
			std::string exe_path = ExePath();

			std::ifstream infile;
			infile.open(exe_path+"\\test_in.csv");
			std::ofstream outfile;
			outfile.open(exe_path + "\\test_out.txt");

			std::string a, b;
			std::string line;
			std::string modelfile = exe_path+"\\model_weights.bin";

		    std::string delimiter =",";

			bool resample1, resample2, norm1, norm2;
			
			//
			while (std::getline(infile, line)) {

				a=line.substr(0, line.find(delimiter));
				b=line.erase(0, line.find(delimiter) + delimiter.length())+".wav";

				printf("%s\t%s\n", a.data(), b.data());


				dtw_dist = dtw_wav((exe_path+"\\references\\"+b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(),true, true, false, false);
				outfile << a << "\t" << b << "\t" << dtw_dist << "\n";
				//outfile << a << "\t" << b << "\t" << "0000"<<"\t"<<dtw_dist<<"\n";
				//dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), true, true, true, false);
				//outfile << a << "\t" << b << "\t" << "0001" << "\t" << dtw_dist << "\n";
				//dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), true, true, false, true);
				//outfile << a << "\t" << b << "\t" << "0010" << "\t" << dtw_dist << "\n";
				//dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), true, true, false, false);
				//outfile << a << "\t" << b << "\t" << "0011" << "\t" << dtw_dist << "\n";
				//dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), true, false, true, true);
				//outfile << a << "\t" << b << "\t" << "0100" << "\t" << dtw_dist << "\n";
				//dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), true, false, true, false);
				//outfile << a << "\t" << b << "\t" << "0101" << "\t" << dtw_dist << "\n";
				//dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), true, false, false, true);
				//outfile << a << "\t" << b << "\t" << "0110" << "\t" << dtw_dist << "\n";
				//dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), true, false, false, false);
				//outfile << a << "\t" << b << "\t" << "0111" << "\t" << dtw_dist << "\n";
				//dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), false, true, true, true);
				//outfile << a << "\t" << b << "\t" << "1000" << "\t" << dtw_dist << "\n";
				//dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), false, true, true, false);
				//outfile << a << "\t" << b << "\t" << "1001" << "\t" << dtw_dist << "\n";
				//dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), false, true, false, true);
				//outfile << a << "\t" << b << "\t" << "1010" << "\t" << dtw_dist << "\n";
				//dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), false, true, false, false);
				//outfile << a << "\t" << b << "\t" << "1011" << "\t" << dtw_dist << "\n";
				//dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), false, false, true, true);
				//outfile << a << "\t" << b << "\t" << "1100" << "\t" << dtw_dist << "\n";
				//dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), false, false, true, false);
				//outfile << a << "\t" << b << "\t" << "1101" << "\t" << dtw_dist << "\n";;
				//dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), false, false, false, true);
				//outfile << a << "\t" << b << "\t" << "1110" << "\t" << dtw_dist << "\n";
				//dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), false, false, false, false);
				//outfile << a << "\t" << b << "\t" << "1111" << "\t" << dtw_dist << "\n";

			}

			printf("\n\n\ntest finished !!!!");

			infile.close();
			outfile.close();
			

		}

		void dtw_wav_vs_pp_at_scale() {

			float dtw_dist = 2;
			std::string exe_path = ExePath();

			std::ifstream infile;
			infile.open(exe_path + "\\test_in.csv");
			std::ofstream outfile;
			outfile.open(exe_path + "\\test_out.txt");

			std::string a, b;
			std::string line;
			std::string modelfile = exe_path + "\\model_weights.bin";

			std::string delimiter = ",";

			bool resample1, resample2, norm1, norm2;

			//
			while (std::getline(infile, line)) {

				a = line.substr(0, line.find(delimiter));
				b = line.erase(0, line.find(delimiter) + delimiter.length()) + ".pp";

				printf("%s\t%s\n", a.data(), b.data());


				dtw_dist = c_dtw_wav_vs_pp( (exe_path + "\\candidates\\" + a).data(), (exe_path + "\\references\\" + b).data(), modelfile.data(), true, true);
				outfile << a << "\t" << b << "\t" << dtw_dist << "\n";
				//outfile << a << "\t" << b << "\t" << "00" << "\t" << dtw_dist << "\n";
				//dtw_dist = c_dtw_wav_vs_pp((exe_path + "\\candidates\\" + a).data(), (exe_path + "\\references\\" + b).data(), modelfile.data(), true, false);
				//outfile << a << "\t" << b << "\t" << "01" << "\t" << dtw_dist << "\n";
				//dtw_dist = c_dtw_wav_vs_pp((exe_path + "\\candidates\\" + a).data(), (exe_path + "\\references\\" + b).data(), modelfile.data(), false, true);
				//outfile << a << "\t" << b << "\t" << "10" << "\t" << dtw_dist << "\n";
				//dtw_dist = c_dtw_wav_vs_pp((exe_path + "\\candidates\\" + a).data(), (exe_path + "\\references\\" + b).data(), modelfile.data(), false, false);
				//outfile << a << "\t" << b << "\t" << "11" << "\t" << dtw_dist << "\n";

			}

			printf("\n\n\ntest finished !!!!");

			infile.close();
			outfile.close();


		}

		void dtw_wav_vs_pp_at_scale_reproducibility() {

			float dtw_dist = 2;
			std::string exe_path = ExePath();

			std::ifstream infile;
			infile.open(exe_path + "\\test_in.csv");
			std::ofstream outfile;
			outfile.open(exe_path + "\\test_out.txt");

			std::string a, b;
			std::string line;
			std::string modelfile = exe_path + "\\model_weights.bin";

			std::string delimiter = ",";

			bool resample1, resample2, norm1, norm2;

			//
			while (std::getline(infile, line)) {

				a = line.substr(0, line.find(delimiter));
				b = line.erase(0, line.find(delimiter) + delimiter.length()) + ".pp";

				printf("%s\t%s\n", a.data(), b.data());


				//dtw_dist = c_dtw_wav_vs_pp((exe_path + "\\candidates\\" + a).data(), (exe_path + "\\references\\" + b).data(), modelfile.data(), true, true);
				//outfile << a << "\t" << b << "\t" << dtw_dist << "\n";
				dtw_dist = c_dtw_wav_vs_pp((exe_path + "\\candidates\\" + a).data(), (exe_path + "\\references\\" + b).data(), modelfile.data(), true, true);
				outfile << a << "\t" << b << "\t" << "00" << "\t" << dtw_dist << "\n";
				dtw_dist = c_dtw_wav_vs_pp((exe_path + "\\candidates\\" + a).data(), (exe_path + "\\references\\" + b).data(), modelfile.data(), true, false);
				outfile << a << "\t" << b << "\t" << "01" << "\t" << dtw_dist << "\n";
				dtw_dist = c_dtw_wav_vs_pp((exe_path + "\\candidates\\" + a).data(), (exe_path + "\\references\\" + b).data(), modelfile.data(), false, true);
				outfile << a << "\t" << b << "\t" << "10" << "\t" << dtw_dist << "\n";
				dtw_dist = c_dtw_wav_vs_pp((exe_path + "\\candidates\\" + a).data(), (exe_path + "\\references\\" + b).data(), modelfile.data(), false, false);
				outfile << a << "\t" << b << "\t" << "11" << "\t" << dtw_dist << "\n";

			}

			printf("\n\n\ntest finished !!!!");

			infile.close();
			outfile.close();


		}

		void dtw_at_scale_reproducibility() {

			float dtw_dist = 2;
			std::string exe_path = ExePath();

			std::ifstream infile;
			infile.open(exe_path + "\\test_in.csv");
			std::ofstream outfile;
			outfile.open(exe_path + "\\test_out.txt");

			std::string a, b;
			std::string line;
			std::string modelfile = exe_path + "\\model_weights.bin";

			std::string delimiter = ",";

			bool resample1, resample2, norm1, norm2;

			//
			while (std::getline(infile, line)) {

				a = line.substr(0, line.find(delimiter));
				b = line.erase(0, line.find(delimiter) + delimiter.length()) + ".wav";

				printf("%s\t%s\n", a.data(), b.data());


				//dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), true, true, false, false);
				//outfile << a << "\t" << b << "\t" << dtw_dist << "\n";
				dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), true, true, true, true);
				outfile << a << "\t" << b << "\t" << "0000"<<"\t"<<dtw_dist<<"\n";
				dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), true, true, true, false);
				outfile << a << "\t" << b << "\t" << "0001" << "\t" << dtw_dist << "\n";
				dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), true, true, false, true);
				outfile << a << "\t" << b << "\t" << "0010" << "\t" << dtw_dist << "\n";
				dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), true, true, false, false);
				outfile << a << "\t" << b << "\t" << "0011" << "\t" << dtw_dist << "\n";
				dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), true, false, true, true);
				outfile << a << "\t" << b << "\t" << "0100" << "\t" << dtw_dist << "\n";
				dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), true, false, true, false);
				outfile << a << "\t" << b << "\t" << "0101" << "\t" << dtw_dist << "\n";
				dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), true, false, false, true);
				outfile << a << "\t" << b << "\t" << "0110" << "\t" << dtw_dist << "\n";
				dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), true, false, false, false);
				outfile << a << "\t" << b << "\t" << "0111" << "\t" << dtw_dist << "\n";
				dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), false, true, true, true);
				outfile << a << "\t" << b << "\t" << "1000" << "\t" << dtw_dist << "\n";
				dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), false, true, true, false);
				outfile << a << "\t" << b << "\t" << "1001" << "\t" << dtw_dist << "\n";
				dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), false, true, false, true);
				outfile << a << "\t" << b << "\t" << "1010" << "\t" << dtw_dist << "\n";
				dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), false, true, false, false);
				outfile << a << "\t" << b << "\t" << "1011" << "\t" << dtw_dist << "\n";
				dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), false, false, true, true);
				outfile << a << "\t" << b << "\t" << "1100" << "\t" << dtw_dist << "\n";
				dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), false, false, true, false);
				outfile << a << "\t" << b << "\t" << "1101" << "\t" << dtw_dist << "\n";;
				dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), false, false, false, true);
				outfile << a << "\t" << b << "\t" << "1110" << "\t" << dtw_dist << "\n";
				dtw_dist = dtw_wav((exe_path + "\\references\\" + b).data(), (exe_path + "\\candidates\\" + a).data(), modelfile.data(), false, false, false, false);
				outfile << a << "\t" << b << "\t" << "1111" << "\t" << dtw_dist << "\n";

			}

			printf("\n\n\ntest finished !!!!");

			infile.close();
			outfile.close();


		}

	}//namespace utilities
}//namespace nni
