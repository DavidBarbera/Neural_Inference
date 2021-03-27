

namespace nni {

	namespace layers {

		template <typename DerivedOut, typename DerivedIn>
		void bidirectional_gru(Eigen::MatrixBase<DerivedOut>& output, Eigen::MatrixBase<DerivedIn>& input, float* &weights, int units, int in, int T) {
			
			//Forward layer:
			float *buffer_fW_z = (float*)weights;
			float *buffer_fW_r = (float*)buffer_fW_z + (in * units);
			float *buffer_fW_h = (float*)buffer_fW_r + (in * units);
			float *buffer_fU_z = (float*)buffer_fW_h + (in * units);
			float *buffer_fU_r = (float*)buffer_fU_z + (units * units);
			float *buffer_fU_h = (float*)buffer_fU_r + (units * units);
			float *buffer_fb_z = (float*)buffer_fU_h + (units * units);
			float *buffer_fb_r = (float*)buffer_fb_z + (units);
			float *buffer_fb_h = (float*)buffer_fb_r + (units);
			//Backward layer:
			float *buffer_bW_z = (float*)buffer_fb_h + (units);
			float *buffer_bW_r = (float*)buffer_bW_z + (in * units);
			float *buffer_bW_h = (float*)buffer_bW_r + (in * units);
			float *buffer_bU_z = (float*)buffer_bW_h + (in * units);
			float *buffer_bU_r = (float*)buffer_bU_z + (units * units);
			float *buffer_bU_h = (float*)buffer_bU_r + (units * units);
			float *buffer_bb_z = (float*)buffer_bU_h + (units * units);
			float *buffer_bb_r = (float*)buffer_bb_z + (units);
			float *buffer_bb_h = (float*)buffer_bb_r + (units);
			//Return pointer to next layer of weights:
			weights = (float*)buffer_bb_h + (units);

			//As Matrices:
			//Forward layer:
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> e_fW_z(buffer_fW_z, in, units);
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> e_fW_r(buffer_fW_r, in, units);
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> e_fW_h(buffer_fW_h, in, units);
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> e_fU_z(buffer_fU_z, units, units);
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> e_fU_r(buffer_fU_r, units, units);
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> e_fU_h(buffer_fU_h, units, units);
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> e_fb_z(buffer_fb_z, 1, units);
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> e_fb_r(buffer_fb_r, 1, units);
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> e_fb_h(buffer_fb_h, 1, units);
			//Backward layer
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> e_bW_z(buffer_bW_z, in, units);
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> e_bW_r(buffer_bW_r, in, units);
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> e_bW_h(buffer_bW_h, in, units);
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> e_bU_z(buffer_bU_z, units, units);
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> e_bU_r(buffer_bU_r, units, units);
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> e_bU_h(buffer_bU_h, units, units);
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> e_bb_z(buffer_bb_z, 1, units);
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> e_bb_r(buffer_bb_r, 1, units);
			Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> e_bb_h(buffer_bb_h, 1, units);

			//for time-step calculations
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > z(1,units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > r(1,units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > hh(1,units);
			
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > hf(T+1, units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > hb(T + 1, units);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > xt(1,in);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > ht_1(1, units);
			//init
			hf.row(0).setZero();
			hb.row(T).setZero();

			//Forward pass:
			for (int t = 0; t < T; t++) {
				xt = input.row(t);
				ht_1 = hf.row(t);
				z = ((xt*e_fW_z) + (ht_1*e_fU_z) + e_fb_z).unaryExpr([](float x) {return hard_sigmoid(x); }); // shouldn't be sigmoid?
				r = ((xt*e_fW_r) + (ht_1*e_fU_r) + e_fb_r).unaryExpr([](float x) {return hard_sigmoid(x); });
				hh = ((xt*e_fW_h) + ((r.cwiseProduct(ht_1))*e_fU_h) + e_fb_h).unaryExpr([](float x) {return relu(x); });
				// shouldnt be: hh = ((xt*e_fW_h) + ((r.cwiseProduct((ht_1*e_fU_h) + e_fb_h)))).unaryExpr([](float x) {return relu(x); }); ?
				// otherwise we can't reproduce it in PyTorch. Also, PyTorch has the activations fused to sigmoid and tanh
				//https://github.com/pytorch/pytorch/blob/bba30d1bd82492cd8b22c34775e27bf3b3d983c9/aten/src/ATen/native/cuda/RNN.cu
				hf.row(t + 1) = ((z.cwiseProduct(ht_1)) + (z.unaryExpr([](float x) {return one_minus(x); }).cwiseProduct(hh)));
			}

			//Backward pass:
			for (int t = T - 1; t >= 0; t--) {
				xt = input.row(t);
				ht_1 = hb.row(t + 1);
				z = ((xt*e_bW_z) + (ht_1*e_bU_z) + e_bb_z).unaryExpr([](float x) {return hard_sigmoid(x); });
				r = ((xt*e_bW_r) + (ht_1*e_bU_r) + e_bb_r).unaryExpr([](float x) {return hard_sigmoid(x); });
				hh = ((xt*e_bW_h) + ((r.cwiseProduct(ht_1))*e_bU_h) + e_bb_h).unaryExpr([](float x) {return relu(x); });
				hb.row(t) = ((z.cwiseProduct(ht_1)) + (z.unaryExpr([](float x) {return one_minus(x); }).cwiseProduct(hh)));
			}

			output.block(0, 0, T, units) = hf.block(1, 0, T, units);
			output.block(0, units, T, units) = hb.block(0, 0, T, units);

		}

	}//namespace layers

}//namespace nni