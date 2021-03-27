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
//#include <assert.h> //




#include <algorithm>

//#include <cmath>
//#include <cstring>
//#include <stdlib.h>
//#include <assert.h>
//#include <stdio.h>
//#include <cstdio>
//#include <errno.h>
//#include <string.h>
//#include <math.h>
//#include <ctype.h>



//#include <Eigen/Dense>
//#include <unsupported/Eigen/CXX11/Tensor>
//#include "/../third_party/eigen/Eigen/Dense"

#include "../third_party/eigen/Eigen/Dense"

//#define BOOST_DATE_TIME_NO_LIB
//////In Android this is way trickier than expected, if there is no need for it then just discard metrics.h
//#include <boost/interprocess/file_mapping.hpp>
//#include <boost/interprocess/mapped_region.hpp>
//#include "thirdparty/boost/boost/interprocess/file_mapping.hpp"
//#include "thirdparty/boost/boost/interprocess/mapped_region.hpp"







namespace nni {

	namespace activations{}
	namespace layers{}
	namespace DTW {}
	namespace mfcc {}
	namespace riffio {}
	namespace resample {}
	namespace gotcha {}
	namespace plugin {}
	namespace metrics {}
	
	//namespace tests {}

	//using namespace boost::interprocess;
	using namespace activations;
	using namespace layers;
	using namespace DTW;
	using namespace mfcc;
	using namespace riffio;
	using namespace resample;
	using namespace gotcha;
	using namespace plugin;
	using namespace metrics;
	
	//using namespace tests;
	
	

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

//Plugin stuff
#include "plugin/plugin.h"



//Testing
//#include "tests/test.h"
//#include "tests/test_wav_to_posteriorgrams.h"
//#include "tests/test_samples_to_posteriorgrams.h"

//Utilities
#ifdef _WIN32
#include "utilities/windows_wav_to_mfcc.h"
#endif

#endif