/* tests use of random number generators  */
/*	C.A. Bertulani - 03/20/2000
/**************************************************************/
typedef double Number;

#include<iostream>
#include<math.h>
using namespace std;
int main() {
    /* Create arrays with elements from a pre-defined formula  */
	Number random(long*);
	Number random2(long*);
	Number random3(long*);
	long n,i;
	n = -1;
	for(i=1;i<10;i++){
		cout << random(&n) << endl;
		n = 2*n;
	}
	cout << endl;
    n=-1; 
	for(i=1;i<10;i++){
		cout << random2(&n) << endl;
		n = 2*n;
	}
    n=-1; 
	for(i=1;i<10;i++){
		cout << random3(&n) << endl;
		n = 2*n;
	}
	return 0;
}

/**************************************************************/
/*  routine random   */                                           
/*  Generates random numbers between 0.0 and 1.0. Call with idum
	a negative number to initialize; thereafter, do not alter
	idum between successive deviates in a sequence. RNMX 
	approximates the largest floating value that is less than 1.\
	The period of this routine is ~ 10^8
	C.A. Bertulani - 03/20/2000
*/
/**************************************************************/
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

Number random(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	Number temp;

	if (*idum <= 0 || !iy) {			/* initialize  */
		if (-(*idum) < 1) *idum=1;		/* Be sure to prevent idum=0 */
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {		/* load the shuffle table after 8 warm-ups) */
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;						/* Start here when not initializing  */
	*idum=IA*(*idum-k*IQ)-IR*k;			/* Compute idum=(IA*idum) % IM without */
	if (*idum < 0) *idum += IM;			/* overflows by Schrage's method */
	j=iy/NDIV;							/* Will be in the range 0..NTAB-1. */
	iy=iv[j];							/* Output previously stored value and */
	iv[j] = *idum;						/* refill the suffle table  */ 
	if ((temp=AM*iy) > RNMX) return RNMX; /* Because users don't expect endpoint values */
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

/*     End of random routine   */
/**************************************************************/
/**************************************************************/
/*  routine random2   */
/*  Generates random numbers between 0.0 and 1.0 (exclusive of 
    the end point values. Has a longer period (> 2 x 10^18) than 
	routine random. Call with idum a negative number to 
	initialize; thereafter, do not alter idum between successive 
	deviates in a sequence. RNMX  approximates the largest 
	floating value that is less than 1.
	C.A. Bertulani - 03/20/2000
*/
/**************************************************************/
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

Number random2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	Number temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
/*     End of random2 routine   */
/**************************************************************/
/**************************************************************/
/*  routine random3   */
/*  Generates random numbers between 0.0 and 1.0 (exclusive of 
    the end point values. Call with idum a negative number to 
	initialize, or reinitialize sequence. It is 2 times faster
	than random and 3 times faster than random2.
	Based on an algorithm described by D.E. Knuth, "Seminumerical
    Algorithms" (Wesley, MA). 
	C.A. Bertulani - 03/20/2000
*/
/**************************************************************/
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

Number random3(long *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
/*     End of random3 routine   */
/**************************************************************/