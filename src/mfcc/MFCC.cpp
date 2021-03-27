#include "MFCC.h"


namespace nni {
	namespace mfcc {

		CMFCC::CMFCC(int ncoeff, int energy, int winsize, double lofreq, double hifreq, double srate)
		{
			m_ncoeff = ncoeff;
			m_energy = energy;
			m_winsize = winsize;

			// calculate a good FFT size
			m_fftsize = 1024;
			while (m_fftsize < winsize) m_fftsize *= 2;

			// design the filterbank
			designfbank(lofreq, hifreq, m_fftsize, srate);

			// design a hamming window
			m_wsp = new float[winsize];
			for (int i = 0; i < winsize; i++) m_wsp[i] = (0.54 - 0.46*cos(2 * M_PI*i / (m_winsize - 1)));

			// processing buffers
			m_fsp = new float[m_fftsize + 2];
			m_chan = new float[m_nfilter];

			// AMFCC calculation
			// autocorrelation
			//m_ac = new float[m_fftsize];

			// Kaiser window
			////m_kwin = new float[m_fftsize];
			////memset(m_kwin, 0, m_fftsize * sizeof(float));
			////int kskip = (int)(0.5 + 0.002*srate);			// 2ms of lags dropped (see Shannon & Palliwal 2006)
			////KaiserWindow(m_kwin + kskip, winsize - kskip, 11.3);	// 80dB range
		//	TRACE("KaiserWindow=");
		//	for (int i=0;i<winsize;i++) TRACE("%d. %g\n",i,m_kwin[i]);

		}


		CMFCC::~CMFCC(void)
		{
			for (int i = 0; i < m_nfilter; i++) free(m_filtab[i].win);
			free(m_filtab);
			delete[] m_wsp;
			delete[] m_fsp;
			delete[] m_chan;
			//delete[] m_ac;
			//delete[] m_kwin;
		}


		static double mylog10(double val)
		{
			return (val < 1.0E-10) ? -10 : log10(val);
		}
		static double HertzToMel(double hz)
		{
			return (1000.0 / log(2.0)) * log(1.0 + hz / 1000.0);
		}
		static double MelToHertz(double mel)
		{
			return 1000 * (exp((log(2.0) / 1000)*mel) - 1.0);
		}

		/* design filter bank */
		int CMFCC::designfbank(double lf, double hf, int nfft, int srate)
		{
			int	i, j;
			double	fhz, fmel;
			double	spacing = 100;
			double	width2 = 150;

			/* get number of filters required 300 mel width, 100 mel spacing */
			lf = MelToHertz(HertzToMel(lf) + width2);
			hf = MelToHertz(HertzToMel(hf) - width2);
			fhz = lf;
			m_nfilter = 0;
			while (fhz <= hf) {
				fmel = HertzToMel(fhz);
				fhz = MelToHertz(fmel + spacing);
				m_nfilter++;
			}
			m_filtab = (struct filter_rec *)calloc(m_nfilter, sizeof(struct filter_rec));

			/* fit triangular filters between limits */
			lf = MelToHertz(HertzToMel(lf) + width2);
			hf = MelToHertz(HertzToMel(hf) - width2);
			fhz = lf;
			m_nfilter = 0;
			while (fhz <= hf) {
				fmel = HertzToMel(fhz);
				m_filtab[m_nfilter].lidx = (int)(0.5 + nfft * (MelToHertz(fmel - width2) / srate));
				m_filtab[m_nfilter].cidx = (int)(0.5 + nfft * (MelToHertz(fmel) / srate));
				m_filtab[m_nfilter].hidx = (int)(0.5 + nfft * (MelToHertz(fmel + width2) / srate));
				m_filtab[m_nfilter].win = (double *)calloc(nfft / 2, sizeof(double));
				fhz = MelToHertz(fmel + spacing);
				m_nfilter++;
			}

			/* create triangular weighting windows */
			for (i = 0; i < m_nfilter; i++) {
				for (j = m_filtab[i].lidx; j < m_filtab[i].cidx; j++) {
					m_filtab[i].win[j] = 1.0 - (double)(m_filtab[i].cidx - j) / (double)(m_filtab[i].cidx - m_filtab[i].lidx + 1);
				}
				for (j = m_filtab[i].cidx; j <= m_filtab[i].hidx; j++) {
					m_filtab[i].win[j] = 1.0 - (double)(j - m_filtab[i].cidx) / (double)(m_filtab[i].hidx - m_filtab[i].cidx + 1);
				}
			}

			return(m_nfilter);
		}

		/* FOUR1 - Numerical Recipes in C */

		/*
		FFT routine based on one originally written by N. Brenner of Lincoln
		Laboratories.
		On entry: 'data[1..2*nn]' contains 'nn' complex data points, with each
				  complex value occupying two consecutive locations, real
				  value stored in the lower address.  'nn' MUST be an integer
				  power of 2 (this is not checked for).
		On exit:  'data[1..2*nn]' contains the complex transformation output,
				  again stored in pairs of real and imaginary values.
		'isign' determins the direction of the transform:
				   1 for forward DFT,
				  -1 for inverse DFT, scaled by a factor of 'nn';
									  (i.e. not normalised by 1/N).
		*/

#define FORWARD 1
#define REVERSE -1

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

		static void four1(float *data, int nn, int isign)
		{
			int n, mmax, m, j, istep, i;
			double wtemp, wr, wpr, wpi, wi, theta; /* double precision for the */
											  /* trigonometric recurrences. */
			float tempr, tempi;

			n = nn << 1;
			j = 1;

			/* Bit-Reversal routine */
			for (i = 1; i < n; i += 2) {
				if (j > i) {
					SWAP(data[j], data[i]);   /* exchange the two complex numbers */
					SWAP(data[j + 1], data[i + 1]);
				}
				m = n >> 1;
				while (m >= 2 && j > m) {
					j -= m;
					m >>= 1;
				}
				j += m;
			}

			/* Danielson-Lanczos loops */
			mmax = 2;
			while (n > mmax) {
				/* outer loop executed ln(nn) times */
				istep = 2 * mmax;

				/* initialise for the trigonometric recurrence */
				theta = 2 * M_PI / (isign*mmax);
				wtemp = sin(0.5*theta);
				wpr = -2.0*wtemp*wtemp;
				wpi = sin(theta);
				wr = 1.0;
				wi = 0.0;

				/* two nested inner loops */
				for (m = 1; m < mmax; m += 2) {
					for (i = m; i <= n; i += istep) {

						/* Danielson-Lanczos formula: */
						j = i + mmax;
						tempr = (float)(wr*data[j] - wi * data[j + 1]);
						tempi = (float)(wr*data[j + 1] + wi * data[j]);
						data[j] = data[i] - tempr;
						data[j + 1] = data[i + 1] - tempi;
						data[i] += tempr;
						data[i + 1] += tempi;
					}

					/* trigonometric recurrence */
					wr = (wtemp = wr)*wpr - wi * wpi + wr;
					wi = wi * wpr + wtemp * wpi + wi;
				}
				mmax = istep;
			}
		}

#undef SWAP

		/* REALFT - Numerical Recipes in C */
		/*
		This routine takes a single real function and calculates either its
		DFT or its Inverse DFT depending on 'isign' as for 'four1'.
		On Entry: 'data[0..2n-1]' contains the real function of length '2n',
				  where 'n' MUST be an integral power of 2.
		On Exit:  'data[0..2n-1]' contains the positive half of the complex DFT.
				  'data[0]' contains the real valued first component and
				  'data[1]' contains the real valued last component of the
				  complex transform.
		This routine also calculates the inverse transform of a complex data
		array if it is the transform of real data and 'isign' is set to -1.
		(Result in this case must be multiplied by 1/n.)
		*/

		static void REALFFT(float *data, int n, int isign)
		{
			int i, i1, i2, i3, i4, n2p3;
			double c1 = (double)0.5, c2, h1r, h1i, h2r, h2i;
			double wr, wi, wpr, wpi, wtemp, theta;

			/* the code in this routine expects indexes to go 1..2n */
			data--;

			/* Initialise the recurrence. */
			theta = M_PI / (double)n;
			if (isign == 1) {
				c2 = (double)-0.5;
				four1(data, n, 1);                 /* Forward Transform */
			}
			else {
				c2 = (double)0.5;
				theta = -theta;
			}

			wtemp = sin(0.5*theta);
			wpr = -2.0*wtemp*wtemp;
			wpi = sin(theta);
			wr = 1.0 + wpr;
			wi = wpi;
			n2p3 = 2 * n + 3;

			for (i = 2; i <= n / 2; i++) {        /* Case i==1 done below */

			  /* Separate transforms from the data */
				i4 = 1 + (i3 = n2p3 - (i2 = 1 + (i1 = i + i - 1)));
				h1r = c1 * (data[i1] + data[i3]);
				h1i = c1 * (data[i2] - data[i4]);
				h2r = -c2 * (data[i2] + data[i4]);
				h2i = c2 * (data[i1] - data[i3]);

				/* Recombine to form true transform of original data */
				data[i1] = (float)(h1r + wr * h2r - wi * h2i);
				data[i2] = (float)(h1i + wr * h2i + wi * h2r);
				data[i3] = (float)(h1r - wr * h2r + wi * h2i);
				data[i4] = (float)(-h1i + wr * h2i + wi * h2r);

				/* The recurrence Relations */
				wr = (wtemp = wr)*wpr - wi * wpi + wr;
				wi = wi * wpr + wtemp * wpi + wi;
			}

			if (isign == 1) {
				/* Squeeze the first & last data into the original array */
				data[1] = (float)((h1r = data[1]) + data[2]);
				data[2] = (float)(h1r - data[2]);
			}
			else {
				data[1] = (float)(c1*((h1r = data[1]) + data[2]));
				data[2] = (float)(c1*(h1r - data[2]));
				four1(data, n, -1);                /* Inverse Transform */
			}
		}

		//void realft(float data[], unsigned long n, int isign)
		//{
		//	void four1(float data[], unsigned long nn, int isign);
		//	unsigned long i, i1, i2, i3, i4, np3;
		//	float c1 = 0.5, c2, h1r, h1i, h2r, h2i;
		//	double wr, wi, wpr, wpi, wtemp, theta;

		//	theta = 3.141592653589793 / (double)(n >> 1);
		//	if (isign == 1) {
		//		c2 = -0.5;
		//		four1(data, n >> 1, 1);
		//	}
		//	else {
		//		c2 = 0.5;
		//		theta = -theta;
		//	}
		//	wtemp = sin(0.5 * theta);
		//	wpr = -2.0 * wtemp * wtemp;
		//	wpi = sin(theta);
		//	wr = 1.0 + wpr;
		//	wi = wpi;
		//	np3 = n + 3;
		//	for (i = 2; i <= (n >> 2); i++) {
		//		i4 = 1 + (i3 = np3 - (i2 = 1 + (i1 = i + i - 1)));
		//		h1r = c1 * (data[i1] + data[i3]);
		//		h1i = c1 * (data[i2] - data[i4]);
		//		h2r = -c2 * (data[i2] + data[i4]);
		//		h2i = c2 * (data[i1] - data[i3]);
		//		data[i1] = h1r + wr * h2r - wi * h2i;
		//		data[i2] = h1i + wr * h2i + wi * h2r;
		//		data[i3] = h1r - wr * h2r + wi * h2i;
		//		data[i4] = -h1i + wr * h2i + wi * h2r;
		//		wr = (wtemp = wr) * wpr - wi * wpi + wr;
		//		wi = wi * wpr + wtemp * wpi + wi;
		//	}
		//	if (isign == 1) {
		//		data[1] = (h1r = data[1]) + data[2];
		//		data[2] = h1r - data[2];
		//	}
		//	else {
		//		data[1] = c1 * ((h1r = data[1]) + data[2]);
		//		data[2] = c1 * (h1r - data[2]);
		//		four1(data, n >> 1, -1);
		//	}
		//}
		///* (C) Copr. 1986-92 Numerical Recipes Software 7MZ9%"W5:!+). */


		// zeroth order modified Bessel function of the first kind
		const int ITERLIMIT = 20;
		const double CONVERGE = 1.0E-9;

		double CMFCC::besselI0(double p)
		{
			// initialise iterative loop
			double n = 1;
			double t = 1;
			double d = 1;

			// iteration
			int k = 1;
			double v;

			p = p / 2;

			do {
				n = n * p;
				d = d * k;
				v = n / d;
				t = t + v * v;
			} while ((++k < ITERLIMIT) && (v > CONVERGE));

			return t;
		}

		void CMFCC::KaiserWindow(float *win, int len, double alpha)
		{
			int	hlen = len / 2;
			double d = besselI0(alpha);

			if (len & 1) {
				// calculate window 1..hlen+1
				win[hlen] = 1.0;
				for (int n = 1; n <= hlen; n++) {
					double v = n / (double)hlen;
					win[hlen - n] = win[hlen + n] = besselI0(alpha*sqrt(1 - v * v)) / d;
				}
			}
			else {
				// calculate window 0.5..hlen+0.5
				for (int n = 0; n < hlen; n++) {
					double v = 2 * (n + 0.5) / (double)(len - 1);
					win[hlen - 1 - n] = win[hlen + n] = besselI0(alpha*sqrt(1 - v * v)) / d;
				}
			}
		}

		/* compute autocorrelation */
		void CMFCC::autocorrelation(float *x, float *rr, int n)
		{
			int		i, k;
			double	sum;

			for (i = 0; i < n; i++) rr[i] = 0;

			/* compute autocorrelations */
			for (i = 0; i < n; i++) {
				sum = 0;
				for (k = 0; k < (n - i); k++) sum += x[k] * x[k + i];
				rr[i] = (float)sum;
			}
		}

		// calculate the MFCC coefficients
		void CMFCC::Calc(float *sig, float *data)
		{
			int	j;
			/* pre-emphasise and apply Hamming window */
			m_fsp[0] = 0;
			for (j = 1; j < m_winsize; j++)
				m_fsp[j] = (sig[j] - MFCC_PREEMP * sig[j - 1]) * m_wsp[j];
			for (; j < m_fftsize; j++) m_fsp[j] = 0.0;
			/* calculate energy */
			if (m_energy) {
				double sumsq = 0;
				for (j = 0; j < m_winsize; j++) sumsq += m_fsp[j] * m_fsp[j];
				data[m_ncoeff] = mylog10(sumsq / m_winsize);		// deliberately in bels
			}
			/* do forward FFT */
			REALFFT(m_fsp, m_fftsize / 2, FORWARD);
			/* create magnitude */
			for (j = 0; j < m_fftsize / 2; j++) m_fsp[j] = sqrt(m_fsp[2 * j] * m_fsp[2 * j] + m_fsp[2 * j + 1] * m_fsp[2 * j + 1]);
			/* form filter outputs */
			for (j = 0; j < m_nfilter; j++) {
				double sum = 0.0;
				for (int k = m_filtab[j].lidx; k <= m_filtab[j].hidx; k++)
					sum += m_filtab[j].win[k] * m_fsp[k];
				m_chan[j] = mylog10(sum);
			}
			/* cosine transform */
			//DCT - III (Direct Cosine Transform Type III)
			//to decorrelate highly correlated filter bank coefficients
			double fnorm = sqrt(2.0 / m_nfilter);
			for (j = 0; j < m_ncoeff; j++) {
				double omega = M_PI * (j + 1) / m_nfilter;
				double sum = 0;
				for (int k = 0; k < m_nfilter; k++)
					sum += m_chan[k] * cos((k + 0.5) * omega);
				data[j] = sum * fnorm;
			}
			/* liftering */
			for (j = 0; j < m_ncoeff; j++)
				data[j] *= 1.0 + MFCC_LIFTER * sin((j + 1)*M_PI / MFCC_LIFTER) / 2;

		}

		// calculate AMFCC coefficients
		void CMFCC::CalcAMFCC(float *sig, float *data)
		{
			int	j;
			/* pre-emphasise and apply Hamming window */
			m_fsp[0] = 0;
			for (j = 1; j < m_winsize; j++)
				m_fsp[j] = (sig[j] - MFCC_PREEMP * sig[j - 1]) * m_wsp[j];
			for (; j < m_fftsize; j++) m_fsp[j] = 0.0;
			/* calculate energy */
			if (m_energy) {
				double sumsq = 0;
				for (j = 0; j < m_winsize; j++) sumsq += m_fsp[j] * m_fsp[j];
				data[m_ncoeff] = mylog10(sumsq / m_winsize);		// deliberately in bels
			}
			/* calculate autocorrelation */
			autocorrelation(m_fsp, m_ac, m_winsize);
			/* window autocorrelation */
			for (j = 0; j < m_fftsize; j++) m_ac[j] = (float)(m_ac[j] * m_kwin[j]);
			/* do forward FFT */
			REALFFT(m_ac, m_fftsize / 2, FORWARD);
			/* create magnitude */
			for (j = 0; j < m_fftsize / 2; j++) m_fsp[j] = sqrt(m_ac[2 * j] * m_ac[2 * j] + m_ac[2 * j + 1] * m_ac[2 * j + 1]);
			/* form filter outputs */
			for (j = 0; j < m_nfilter; j++) {
				double sum = 0.0;
				for (int k = m_filtab[j].lidx; k <= m_filtab[j].hidx; k++)
					sum += m_filtab[j].win[k] * m_fsp[k];
				m_chan[j] = mylog10(sum);
			}
			/* cosine transform */
			//DCT - III (Direct Cosine Transform Type III)
			double fnorm = sqrt(2.0 / m_nfilter);
			for (j = 0; j < m_ncoeff; j++) {
				double omega = M_PI * (j + 1) / m_nfilter;
				double sum = 0;
				for (int k = 0; k < m_nfilter; k++)
					sum += m_chan[k] * cos((k + 0.5) * omega);
				data[j] = sum * fnorm;
			}
			/* liftering */
			for (j = 0; j < m_ncoeff; j++)
				data[j] *= 1.0 + MFCC_LIFTER * sin((j + 1)*M_PI / MFCC_LIFTER) / 2;

		}
	}

}
