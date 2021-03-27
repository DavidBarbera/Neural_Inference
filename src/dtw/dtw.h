
#define inf 1.0E+10

namespace nni {
	namespace DTW {

		template <typename Derived>
		void print_eigen(const Eigen::MatrixBase<Derived>& M) {
			//std::cout << M << std::endl;
			Eigen::IOFormat precision(9);
			std::ofstream file;
			file.open("files/dtw_test.txt", std::ios::out);
			file << M.format(precision) << "\n";
			file.close();
		}

		template <typename DerivedD, typename DerivedR, typename DerivedC>
		void rxc_neglog_distance(Eigen::MatrixBase<DerivedD>& D, const Eigen::MatrixBase<DerivedR>& M, const Eigen::MatrixBase<DerivedC>& N) {
			const int r = D.rows()-1;
			const int c = D.cols()-1;

			for (int i = 0; i < r; i++) {
				for (int j = 0; j < c; j++) {
					//printf("(%i,%i)", i, j);
					D(i+1, j+1) = -std::log(M.row(i).cwiseProduct(N.row(j)).sum()); //test with eigen's cwise log
				}
			}
		}

		inline int min_of2(int a, int b) {
			return (a < b) ? a : b;
		}

		inline int argmax_of3(float a,float b,float c)
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

		template <typename Derived>
		void traceback(const Eigen::MatrixBase<Derived>& D) {
			int i = D.rows() - 2;
			int j = D.cols() - 2;
			std::vector<int> p;
			std::vector<int> q;
			p.push_back(i);
			q.push_back(j);


			int tb = 0;
			while (i > 0 || j > 0) {
				tb = argmax_of3(D(i, j), D(i, j + 1), D(i + 1, j));
				switch(tb){
				case 0:
					i--;
					j--;
					break;
				case 1:
					i--;
					break;
				default:
					j--;
					break;
				}
				p.insert(p.begin(), i);
				q.insert(q.begin(), j);
				
			}

		}

		template <typename DerivedM, typename DerivedN>
		float dtw(const Eigen::MatrixBase<DerivedM>& M, const Eigen::MatrixBase<DerivedN>& N) {
			const int r = M.rows();
			const int c = N.rows();
		
			const int warp = 1;
			std::vector<float> min_list;
			int a_[2],b_[2],c_[2];


			//printf("r=%i,c=%i\n", r, c);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> D0(r + 1, c + 1);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> D1(r, c);
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> C(r, c);
			//init
			D1.setZero();
			D0.setConstant(inf);
			D0(0, 0) = 0.f;
			//init algorithm
			D0.block(1, 1, r, c) = D1;
		    rxc_neglog_distance(D0,M, N); //calculate all neg-log on inner product distances for M and N
			D1 = D0.block(1, 1, r, c);

			//C = D1;
			//algorithm
			for (int i = 0; i < r; i++) {
				for (int j = 0; j < c; j++) {
					min_list.push_back(D0(i, j));
					for (int k = 1; k < (warp + 1); k++) {
						min_list.push_back( D0(min_of2(i + k, r - 1), j) );
						min_list.push_back( D0(i, min_of2(j + k, c - 1)) );
					}
					D1(i, j) = D1(i,j)+*std::min_element(min_list.begin(), min_list.end());
					D0(i+1, j+1) = D1(i, j);
					min_list.clear();
				}
			}
			//print_eigen(D1);
			//printf("DTW-> D[-1,-1]: %.9f", D1(r - 1, c - 1));
			return D1(r-1,c-1)/(r+c);
		}

	}//namespace dtw
}//namespace nni
