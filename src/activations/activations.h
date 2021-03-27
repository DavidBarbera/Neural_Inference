
namespace nni {
	namespace activations {

		inline float hard_sigmoid(const float x) {
			const float z = 0.2f * x + 0.5f;
			return z <= 0.0f ? 0.0f : (z <= 1.0f ? z : 1.0f);
		}

		inline float relu(const float x) { return (x > 0.0f) ? x : 0.0f; }

		inline float one_minus(const float x) { return 1.0f - x; };

		template <typename Derived>
		inline Eigen::Matrix<float, 1, Eigen::Dynamic, Eigen::RowMajor> softmax( const Eigen::MatrixBase<Derived>& m) {

			   using std::exp;

			   Eigen::Matrix<float, 1, Eigen::Dynamic> theta(m.size());
  
			   float sum(0.0);
			   float max_v = m.maxCoeff();

			   for (int i = 0; i < m.size(); ++i) {
				     theta(i) = exp(m(i) - max_v);  
				     sum += theta(i);               
			   }

			   for (int i = 0; i < m.size(); ++i)
				     theta(i) /= sum;

		   return theta;

		}

	}//namespace activations
}//namespace nni