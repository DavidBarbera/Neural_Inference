#ifndef NNI_INCLUDED
#define NNI_INCLUDED

#include <string>
#include <chrono>
#include <vector>
#include <iostream>
#include <numeric>
#include <functional>
#include <fstream>
#include <ctime>
#include <sstream>
#include <iomanip>
#include <locale>
#include <codecvt>
#include <algorithm>

#include "../third_party/eigen/Eigen/Dense"

namespace nni {

	namespace activations{}
	namespace layers{}
	namespace DTW {}
	namespace mfcc {}
	namespace riffio {}
	namespace resample {}
	namespace gotcha {}
	namespace metrics {}
	
	using namespace activations;
	using namespace layers;
	using namespace DTW;
	using namespace mfcc;
	using namespace riffio;
	using namespace resample;
	using namespace gotcha;
	using namespace metrics;
	

#ifdef _WIN32
	namespace utilities {}
	using namespace utilities;
#endif
}

//Activations
#include "activations/activations.h"

//Layers
#include "layers/bidirectional_gru.h"
#include "layers/batch_normalization.h"
#include "layers/timedistributed_dense.h"

//DTW
#include "dtw/dtw.h"

//MFCC
#include "mfcc/MFCC.h"

//RIFF-I/O
#include "riffio/RIFFIO.h"
#include "riffio/Resample.h"

//Gotcha
#include "gotcha/samples_to_pp.h"
#include "gotcha/wav_to_mfcc.h"

//Metrics
//#include "metrics/metrics.h"

//Utilities
#ifdef _WIN32
#include "utilities/windows_wav_to_mfcc.h"
#endif

#endif