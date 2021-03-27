

namespace nni {

	namespace layers {

		template <typename DerivedOut, typename DerivedIn>
		void timedistributed_dense(Eigen::MatrixBase<DerivedOut>& output, Eigen::MatrixBase<DerivedIn>& input, float* &weights, int in, int out, int T) {
			
			float *buffer_W = weights;
			float *buffer_bias = buffer_W + (in * out);
			weights = buffer_bias + out;

			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> W(buffer_W, in, out);
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> bias(buffer_bias, 1, out);

			//Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> exp(1, out);

			for (int t = 0; t < T; t++) {
				//auto e = (input.row(t)*W + bias).exp();
				////auto invsum = exp.sum().cwiseInverse();
				//output.row(t) = exp.cwiseProduct(invsum);
				//exp(e);
				output.row(t) = softmax(input.row(t)*W + bias);
			}

		}

	}//namespace layers

}//namespace nni