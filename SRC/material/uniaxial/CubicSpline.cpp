// #include "StdAfx.h"
#include "CubicSpline.h"
#include "TriMatrix.h"
 #include <float.h>
#include<math.h>
//#include<conio.h>



        double CubicSpline::EvalT(double x)
        {
			if(xs[0] == 0 && xs[1] == 0 && xs[2] == 0 && xs[3])
			{
				return 10e8;////printf("kkkkkkkkk");
			}
            int i;
            for (i = 0; i < xsL; i++)
            {
                if (xs[i] != 0) break;
            }
			////printf("%d\n\n",i);
            if (i == xsL) return 10e8f;
            double Eps = 0.01f;
            if (x > xs[xsL - 1] - Eps || x < xs[0] + Eps)
            {
                //return 10e8f;
            }
            double t = (Eval(x + Eps) - Eval(x - Eps)) / (Eps * 2);
            return t;
        }

        void CubicSpline::Fit(double* x, int xl, double* y, int yl)
        {
			////printf("jaja?\n");
			xs = new double[xl];
			ys = new double[xl];
			for(int i = 0;i < xl;i++)
			{
				xs[i] = x[i];
				ys[i] = y[i];
			}
			xsL = xl;
			ysL = yl;
            double *dys ,* dxs, *ms;
			dys = new double[xl * 2];
			dxs = new double[xl  * 2];
			ms = new double[xl * 2];
			int dysL = 0,dxsL = 0,msL = 0;
            int length = xl;
            for (int i = 0; i < length - 1; i++)
            {
                double dx = xs[i + 1] - xs[i], dy = ys[i + 1] - ys[i];
                dxs[dxsL] = (dx);dxsL++;
				dys[dysL] = (dy);dysL++;
				ms[msL]= (dy / dx);msL++;
            }

            // Get degree-1 coefficients
            c1s = new double[xl * 2]; c1sc = 0;
            c1s[c1sc] = (ms[0]); c1sc++;
            for (int i = 0; i < dxsL - 1; i++)
            {
                double m = ms[i], mNext = ms[i + 1];
                if (m * mNext <= 0)
                {
                    c1s[c1sc] = (0);c1sc++;
                }
                else
                {
                    double dx = dxs[i], dxNext = dxs[i + 1], common = dx + dxNext;
                    c1s[c1sc] = (3 * common / ((common + dxNext) / m + (common + dx) / mNext)); c1sc++;
                }
            }
            c1s[c1sc] = (ms[msL - 1]); c1sc++;
			c1sL = c1sc;

            // Get degree-2 and degree-3 coefficients
            c2s = new double[xl * 2]; c3s = new double[xl * 2];c2sc = 0;c3sc=0;
            for (int i = 0; i < c1sL - 1; i++)
            {
                double c1 = c1s[i], m = ms[i], invDx = 1 / dxs[i], common = c1 + c1s[i + 1] - m - m;
                c2s[c2sc] = ((m - c1 - common) * invDx); c3s[c3sc] = (common * invDx * invDx);c2sc++;c3sc++;
            }
			c2sL = c2sc;
			c3sL = c3sc;
			sc1sL = 0;
			//sc2sL = c2sc;
			//sc3sL = c3sc;
			////getch();
			Eval(0);
        }


        double CubicSpline::Eval(double x)
        {
			if(xs[0] == 0 && xs[1] == 0 && xs[2] == 0 && xs[3] == 0)return 10e8;
            int i = xsL - 1;
            if (x == xs[i]) { return ys[i]; }

            // Search for the interval x is in, returning the corresponding y if x is one of the original xs
            int low = 0, high = c3sL - 1;
            int mid;
            while (low <= high)
            {
                mid = (int)(0.5 * (low + high));
                double xHere = xs[mid];
                if (xHere < x) { low = mid + 1; }
                else if (xHere > x) { high = mid - 1; }
                else { return ys[mid]; }
            }
			i = (high > 0)?high:0;

            // Interpolate
            double diff = x - xs[i], diffSq = diff * diff;
            return ys[i] + c1s[i] * diff + c2s[i] * diffSq + c3s[i] * diff * diffSq;
        }
	
		CubicSpline::CubicSpline(void)
		{
		}


		CubicSpline::~CubicSpline(void)
		{
		}
