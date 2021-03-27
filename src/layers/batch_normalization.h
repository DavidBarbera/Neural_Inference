

namespace nni {

	namespace layers {

		template <typename DerivedOut, typename DerivedIn>
		void batch_normalization(Eigen::MatrixBase<DerivedOut>& output, Eigen::MatrixBase<DerivedIn>& input, float* &weights, int in, int T) {
			
			float *buffer_gamma = weights;
			float *buffer_beta = buffer_gamma + in;
			float *buffer_m_mean = buffer_beta + in;
			float *buffer_m_var = buffer_m_mean + in;
			weights = buffer_m_var + in;

			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> gamma(buffer_gamma, 1, in);
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> beta(buffer_beta, 1, in);
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> m_mean(buffer_m_mean, 1, in);
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> m_var(buffer_m_var, 1, in);

			float epsilon = 0.001;
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > e_epsilon( 1, in);
			e_epsilon.setConstant(epsilon);

			for (int t = 0; t < T; t++) {
				output.row(t) = (input.row(t) - m_mean).cwiseProduct((m_var + e_epsilon).cwiseSqrt().cwiseInverse()).cwiseProduct(gamma) + beta;
			}
		}

	}//namespace layers

}//namespace nni