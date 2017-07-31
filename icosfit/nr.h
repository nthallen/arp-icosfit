#ifndef _NR_H_
#define _NR_H_
   
#ifndef _FCOMPLEX_DECLARE_T_
typedef struct FCOMPLEX {ICOS_Float r,i;} fcomplex;
#define _FCOMPLEX_DECLARE_T_
#endif /* _FCOMPLEX_DECLARE_T_ */
   
#ifndef _ARITHCODE_DECLARE_T_
typedef struct {
	unsigned long *ilob,*iupb,*ncumfq,jdif,nc,minint,nch,ncum,nrad;
} arithcode;
#define _ARITHCODE_DECLARE_T_
#endif /* _ARITHCODE_DECLARE_T_ */
   
#ifndef _HUFFCODE_DECLARE_T_
typedef struct {
	unsigned long *icod,*ncod,*left,*right,nch,nodemax;
} huffcode;
#define _HUFFCODE_DECLARE_T_
#endif /* _HUFFCODE_DECLARE_T_ */
   
#include <stdio.h>
   
void addint(double **uf, double **uc, double **res, int nf);
void airy(ICOS_Float x, ICOS_Float *ai, ICOS_Float *bi, ICOS_Float *aip, ICOS_Float *bip);
void amebsa(ICOS_Float **p, ICOS_Float y[], int ndim, ICOS_Float pb[],	ICOS_Float *yb,
	ICOS_Float ftol, ICOS_Float (*funk)(ICOS_Float []), int *iter, ICOS_Float temptr);
void amoeba(ICOS_Float **p, ICOS_Float y[], int ndim, ICOS_Float ftol,
	ICOS_Float (*funk)(ICOS_Float []), int *iter);
ICOS_Float amotry(ICOS_Float **p, ICOS_Float y[], ICOS_Float psum[], int ndim,

ICOS_Float (*funk)(ICOS_Float []), int ihi, ICOS_Float fac);
ICOS_Float amotsa(ICOS_Float **p, ICOS_Float y[], ICOS_Float psum[], int ndim, ICOS_Float pb[],
	ICOS_Float *yb, ICOS_Float (*funk)(ICOS_Float []), int ihi, ICOS_Float *yhi, ICOS_Float
fac);
void anneal(ICOS_Float x[], ICOS_Float y[], int iorder[], int ncity);
double anorm2(double **a, int n);
void arcmak(unsigned long nfreq[], unsigned long nchh, unsigned long
nradd,
	arithcode *acode);
void arcode(unsigned long *ich, unsigned char **codep, unsigned long
*lcode,
	unsigned long *lcd, int isign, arithcode *acode);
void arcsum(unsigned long iin[], unsigned long iout[], unsigned long ja,
	int nwk, unsigned long nrad, unsigned long nc);
void asolve(unsigned long n, double b[], double x[], int itrnsp);
void atimes(unsigned long n, double x[], double r[], int itrnsp);
void avevar(ICOS_Float data[], unsigned long n, ICOS_Float *ave, ICOS_Float *var);
void balanc(ICOS_Float **a, int n);
void banbks(ICOS_Float **a, unsigned long n, int m1, int m2, ICOS_Float **al,
	unsigned long indx[], ICOS_Float b[]);
void bandec(ICOS_Float **a, unsigned long n, int m1, int m2, ICOS_Float **al,
	unsigned long indx[], ICOS_Float *d);
void banmul(ICOS_Float **a, unsigned long n, int m1, int m2, ICOS_Float x[], ICOS_Float
b[]);
void bcucof(ICOS_Float y[], ICOS_Float y1[], ICOS_Float y2[], ICOS_Float y12[], ICOS_Float d1,
	ICOS_Float d2, ICOS_Float **c);
void bcuint(ICOS_Float y[], ICOS_Float y1[], ICOS_Float y2[], ICOS_Float y12[],
	ICOS_Float x1l, ICOS_Float x1u, ICOS_Float x2l, ICOS_Float x2u, ICOS_Float x1,
	ICOS_Float x2, ICOS_Float *ansy, ICOS_Float *ansy1, ICOS_Float *ansy2);
void beschb(double x, double *gam1, double *gam2, double *gampl,
	double *gammi);
ICOS_Float bessi(int n, ICOS_Float x);

ICOS_Float bessi0(ICOS_Float x);
ICOS_Float bessi1(ICOS_Float x);
void bessik(ICOS_Float x, ICOS_Float xnu, ICOS_Float *ri, ICOS_Float *rk, ICOS_Float *rip,
	ICOS_Float *rkp);
ICOS_Float bessj(int n, ICOS_Float x);
ICOS_Float bessj0(ICOS_Float x);
ICOS_Float bessj1(ICOS_Float x);
void bessjy(ICOS_Float x, ICOS_Float xnu, ICOS_Float *rj, ICOS_Float *ry, ICOS_Float *rjp,
	ICOS_Float *ryp);
ICOS_Float bessk(int n, ICOS_Float x);
ICOS_Float bessk0(ICOS_Float x);
ICOS_Float bessk1(ICOS_Float x);
ICOS_Float bessy(int n, ICOS_Float x);
ICOS_Float bessy0(ICOS_Float x);
ICOS_Float bessy1(ICOS_Float x);
ICOS_Float beta(ICOS_Float z, ICOS_Float w);
ICOS_Float betacf(ICOS_Float a, ICOS_Float b, ICOS_Float x);
ICOS_Float betai(ICOS_Float a, ICOS_Float b, ICOS_Float x);
ICOS_Float bico(int n, int k);
void bksub(int ne, int nb, int jf, int k1, int k2, ICOS_Float ***c);
ICOS_Float bnldev(ICOS_Float pp, int n, long *idum);
ICOS_Float brent(ICOS_Float ax, ICOS_Float bx, ICOS_Float cx,
	ICOS_Float (*f)(ICOS_Float), ICOS_Float tol, ICOS_Float *xmin);
void broydn(ICOS_Float x[], int n, int *check,
	void (*vecfunc)(int, ICOS_Float [], ICOS_Float []));
void bsstep(ICOS_Float y[], ICOS_Float dydx[], int nv, ICOS_Float *xx, ICOS_Float htry,
	ICOS_Float eps, ICOS_Float yscal[], ICOS_Float *hdid, ICOS_Float *hnext,

void (*derivs)(ICOS_Float, ICOS_Float [], ICOS_Float []));
void caldat(long julian, int *mm, int *id, int *iyyy);
void chder(ICOS_Float a, ICOS_Float b, ICOS_Float c[], ICOS_Float cder[], int n);
ICOS_Float chebev(ICOS_Float a, ICOS_Float b, ICOS_Float c[], int m, ICOS_Float x);
void chebft(ICOS_Float a, ICOS_Float b, ICOS_Float c[], int n, ICOS_Float (*func)(ICOS_Float));
void chebpc(ICOS_Float c[], ICOS_Float d[], int n);
void chint(ICOS_Float a, ICOS_Float b, ICOS_Float c[], ICOS_Float cint[], int n);
ICOS_Float chixy(ICOS_Float bang);
void choldc(ICOS_Float **a, int n, ICOS_Float p[]);
void cholsl(ICOS_Float **a, int n, ICOS_Float p[], ICOS_Float b[], ICOS_Float x[]);
void chsone(ICOS_Float bins[], ICOS_Float ebins[], int nbins, int knstrn,
	ICOS_Float *df, ICOS_Float *chsq, ICOS_Float *prob);
void chstwo(ICOS_Float bins1[], ICOS_Float bins2[], int nbins, int knstrn,
	ICOS_Float *df, ICOS_Float *chsq, ICOS_Float *prob);
void cisi(ICOS_Float x, ICOS_Float *ci, ICOS_Float *si);
void cntab1(int **nn, int ni, int nj, ICOS_Float *chisq,
	ICOS_Float *df, ICOS_Float *prob, ICOS_Float *cramrv, ICOS_Float *ccc);
void cntab2(int **nn, int ni, int nj, ICOS_Float *h, ICOS_Float *hx, ICOS_Float *hy,
	ICOS_Float *hygx, ICOS_Float *hxgy, ICOS_Float *uygx, ICOS_Float *uxgy, ICOS_Float *uxy);
void convlv(ICOS_Float data[], unsigned long n, ICOS_Float respns[], unsigned long
m,
	int isign, ICOS_Float ans[]);
void copy(double **aout, double **ain, int n);
void correl(ICOS_Float data1[], ICOS_Float data2[], unsigned long n, ICOS_Float ans[]);
void cosft(ICOS_Float y[], int n, int isign);
void cosft1(ICOS_Float y[], int n);
void cosft2(ICOS_Float y[], int n, int isign);
void covsrt(ICOS_Float **covar, int ma, int ia[], int mfit);
void crank(unsigned long n, ICOS_Float w[], ICOS_Float *s);
void cyclic(ICOS_Float a[], ICOS_Float b[], ICOS_Float c[], ICOS_Float alpha, ICOS_Float beta,
	ICOS_Float r[], ICOS_Float x[], unsigned long n);

void daub4(ICOS_Float a[], unsigned long n, int isign);
ICOS_Float dawson(ICOS_Float x);
ICOS_Float dbrent(ICOS_Float ax, ICOS_Float bx, ICOS_Float cx,
	ICOS_Float (*f)(ICOS_Float), ICOS_Float (*df)(ICOS_Float), ICOS_Float tol, ICOS_Float *xmin);
void ddpoly(ICOS_Float c[], int nc, ICOS_Float x, ICOS_Float pd[], int nd);
int decchk(char string[], int n, char *ch);
void derivs(ICOS_Float x, ICOS_Float y[], ICOS_Float dydx[]);
ICOS_Float df1dim(ICOS_Float x);
void dfour1(double data[], unsigned long nn, int isign);
void dfpmin(ICOS_Float p[], int n, ICOS_Float gtol, int *iter, ICOS_Float *fret,
	ICOS_Float (*func)(ICOS_Float []), void (*dfunc)(ICOS_Float [], ICOS_Float []));
ICOS_Float dfridr(ICOS_Float (*func)(ICOS_Float), ICOS_Float x, ICOS_Float h, ICOS_Float *err);
void dftcor(ICOS_Float w, ICOS_Float delta, ICOS_Float a, ICOS_Float b, ICOS_Float endpts[],
	ICOS_Float *corre, ICOS_Float *corim, ICOS_Float *corfac);
void dftint(ICOS_Float (*func)(ICOS_Float), ICOS_Float a, ICOS_Float b, ICOS_Float w,
	ICOS_Float *cosint, ICOS_Float *sinint);
void difeq(int k, int k1, int k2, int jsf, int is1, int isf,
	int indexv[], int ne, ICOS_Float **s, ICOS_Float **y);
void dlinmin(ICOS_Float p[], ICOS_Float xi[], int n, ICOS_Float *fret,
	ICOS_Float (*func)(ICOS_Float []), void (*dfunc)(ICOS_Float [], ICOS_Float[]));
double dpythag(double a, double b);
void drealft(double data[], unsigned long n, int isign);
void dsprsax(double sa[], unsigned long ija[], double x[], double b[],
	unsigned long n);
void dsprstx(double sa[], unsigned long ija[], double x[], double b[],
	unsigned long n);
void dsvbksb(double **u, double w[], double **v, int m, int n, double b[],
	double x[]);
void dsvdcmp(double **a, int m, int n, double w[], double **v);

void eclass(int nf[], int n, int lista[], int listb[], int m);
void eclazz(int nf[], int n, int (*equiv)(int, int));
ICOS_Float ei(ICOS_Float x);
void eigsrt(ICOS_Float d[], ICOS_Float **v, int n);
ICOS_Float elle(ICOS_Float phi, ICOS_Float ak);
ICOS_Float ellf(ICOS_Float phi, ICOS_Float ak);
ICOS_Float ellpi(ICOS_Float phi, ICOS_Float en, ICOS_Float ak);
void elmhes(ICOS_Float **a, int n);
ICOS_Float erfcc(ICOS_Float x);
ICOS_Float erff(ICOS_Float x);
ICOS_Float erffc(ICOS_Float x);
void eulsum(ICOS_Float *sum, ICOS_Float term, int jterm, ICOS_Float wksp[]);
ICOS_Float evlmem(ICOS_Float fdt, ICOS_Float d[], int m, ICOS_Float xms);
ICOS_Float expdev(long *idum);
ICOS_Float expint(int n, ICOS_Float x);
ICOS_Float f1(ICOS_Float x);
ICOS_Float f1dim(ICOS_Float x);
ICOS_Float f2(ICOS_Float y);
ICOS_Float f3(ICOS_Float z);
ICOS_Float factln(int n);
ICOS_Float factrl(int n);
void fasper(ICOS_Float x[], ICOS_Float y[], unsigned long n, ICOS_Float ofac, ICOS_Float
hifac,
	ICOS_Float wk1[], ICOS_Float wk2[], unsigned long nwk, unsigned long *nout,
	unsigned long *jmax, ICOS_Float *prob);
void fdjac(int n, ICOS_Float x[], ICOS_Float fvec[], ICOS_Float **df,
	void (*vecfunc)(int, ICOS_Float [], ICOS_Float []));
void fgauss(ICOS_Float x, ICOS_Float a[], ICOS_Float *y, ICOS_Float dyda[], int na);
void fill0(double **u, int n);
void fit(ICOS_Float x[], ICOS_Float y[], int ndata, ICOS_Float sig[], int mwt,
	ICOS_Float *a, ICOS_Float *b, ICOS_Float *siga, ICOS_Float *sigb, ICOS_Float *chi2, ICOS_Float
*q);
void fitexy(ICOS_Float x[], ICOS_Float y[], int ndat, ICOS_Float sigx[], ICOS_Float sigy[],
	ICOS_Float *a, ICOS_Float *b, ICOS_Float *siga, ICOS_Float *sigb, ICOS_Float *chi2, ICOS_Float
*q);

void fixrts(ICOS_Float d[], int m);
void fleg(ICOS_Float x, ICOS_Float pl[], int nl);
void flmoon(int n, int nph, long *jd, ICOS_Float *frac);
ICOS_Float fmin(ICOS_Float x[]);
void four1(ICOS_Float data[], unsigned long nn, int isign);
void fourew(FILE *file[5], int *na, int *nb, int *nc, int *nd);
void fourfs(FILE *file[5], unsigned long nn[], int ndim, int isign);
void fourn(ICOS_Float data[], unsigned long nn[], int ndim, int isign);
void fpoly(ICOS_Float x, ICOS_Float p[], int np);
void fred2(int n, ICOS_Float a, ICOS_Float b, ICOS_Float t[], ICOS_Float f[], ICOS_Float w[],
	ICOS_Float (*g)(ICOS_Float), ICOS_Float (*ak)(ICOS_Float, ICOS_Float));
ICOS_Float fredin(ICOS_Float x, int n, ICOS_Float a, ICOS_Float b, ICOS_Float t[], ICOS_Float f[], ICOS_Float
w[],
	ICOS_Float (*g)(ICOS_Float), ICOS_Float (*ak)(ICOS_Float, ICOS_Float));
void frenel(ICOS_Float x, ICOS_Float *s, ICOS_Float *c);
void frprmn(ICOS_Float p[], int n, ICOS_Float ftol, int *iter, ICOS_Float *fret,
	ICOS_Float (*func)(ICOS_Float []), void (*dfunc)(ICOS_Float [], ICOS_Float []));
void ftest(ICOS_Float data1[], unsigned long n1, ICOS_Float data2[], unsigned long
n2,
	ICOS_Float *f, ICOS_Float *prob);
ICOS_Float gamdev(int ia, long *idum);
ICOS_Float gammln(ICOS_Float xx);
ICOS_Float gammp(ICOS_Float a, ICOS_Float x);
ICOS_Float gammq(ICOS_Float a, ICOS_Float x);
ICOS_Float gasdev(long *idum);
void gaucof(int n, ICOS_Float a[], ICOS_Float b[], ICOS_Float amu0, ICOS_Float x[], ICOS_Float
w[]);
void gauher(ICOS_Float x[], ICOS_Float w[], int n);

void gaujac(ICOS_Float x[], ICOS_Float w[], int n, ICOS_Float alf, ICOS_Float bet);
void gaulag(ICOS_Float x[], ICOS_Float w[], int n, ICOS_Float alf);
void gauleg(ICOS_Float x1, ICOS_Float x2, ICOS_Float x[], ICOS_Float w[], int n);
void gaussj(ICOS_Float **a, int n, ICOS_Float **b, int m);
void gcf(ICOS_Float *gammcf, ICOS_Float a, ICOS_Float x, ICOS_Float *gln);
ICOS_Float golden(ICOS_Float ax, ICOS_Float bx, ICOS_Float cx, ICOS_Float (*f)(ICOS_Float), ICOS_Float tol,
	ICOS_Float *xmin);
void gser(ICOS_Float *gamser, ICOS_Float a, ICOS_Float x, ICOS_Float *gln);
void hpsel(unsigned long m, unsigned long n, ICOS_Float arr[], ICOS_Float heap[]);
void hpsort(unsigned long n, ICOS_Float ra[]);
void hqr(ICOS_Float **a, int n, ICOS_Float wr[], ICOS_Float wi[]);
void hufapp(unsigned long index[], unsigned long nprob[], unsigned long n,
	unsigned long i);
void hufdec(unsigned long *ich, unsigned char *code, unsigned long lcode,
	unsigned long *nb, huffcode *hcode);
void hufenc(unsigned long ich, unsigned char **codep, unsigned long
*lcode,
	unsigned long *nb, huffcode *hcode);
void hufmak(unsigned long nfreq[], unsigned long nchin, unsigned long
*ilong,
	unsigned long *nlong, huffcode *hcode);
void hunt(ICOS_Float xx[], unsigned long n, ICOS_Float x, unsigned long *jlo);
void hypdrv(ICOS_Float s, ICOS_Float yy[], ICOS_Float dyyds[]);
fcomplex hypgeo(fcomplex a, fcomplex b, fcomplex c, fcomplex z);
void hypser(fcomplex a, fcomplex b, fcomplex c, fcomplex z,
	fcomplex *series, fcomplex *deriv);
unsigned short icrc(unsigned short crc, unsigned char *bufptr,
	unsigned long len, short jinit, int jrev);
unsigned short icrc1(unsigned short crc, unsigned char onech);
unsigned long igray(unsigned long n, int is);
void iindexx(unsigned long n, long arr[], unsigned long indx[]);
void indexx(unsigned long n, ICOS_Float arr[], unsigned long indx[]);
void interp(double **uf, double **uc, int nf);
int irbit1(unsigned long *iseed);
int irbit2(unsigned long *iseed);

void jacobi(ICOS_Float **a, int n, ICOS_Float d[], ICOS_Float **v, int *nrot);
void jacobn(ICOS_Float x, ICOS_Float y[], ICOS_Float dfdx[], ICOS_Float **dfdy, int n);
long julday(int mm, int id, int iyyy);
void kendl1(ICOS_Float data1[], ICOS_Float data2[], unsigned long n, ICOS_Float *tau,
ICOS_Float *z,
	ICOS_Float *prob);
void kendl2(ICOS_Float **tab, int i, int j, ICOS_Float *tau, ICOS_Float *z, ICOS_Float *prob);
void kermom(double w[], double y, int m);
void ks2d1s(ICOS_Float x1[], ICOS_Float y1[], unsigned long n1,
	void (*quadvl)(ICOS_Float, ICOS_Float, ICOS_Float *, ICOS_Float *, ICOS_Float *, ICOS_Float *),
	ICOS_Float *d1, ICOS_Float *prob);
void ks2d2s(ICOS_Float x1[], ICOS_Float y1[], unsigned long n1, ICOS_Float x2[], ICOS_Float
y2[],
	unsigned long n2, ICOS_Float *d, ICOS_Float *prob);
void ksone(ICOS_Float data[], unsigned long n, ICOS_Float (*func)(ICOS_Float), ICOS_Float *d,
	ICOS_Float *prob);
void kstwo(ICOS_Float data1[], unsigned long n1, ICOS_Float data2[], unsigned long
n2,
	ICOS_Float *d, ICOS_Float *prob);
void laguer(fcomplex a[], int m, fcomplex *x, int *its);
void lfit(ICOS_Float x[], ICOS_Float y[], ICOS_Float sig[], int ndat, ICOS_Float a[], int
ia[],
	int ma, ICOS_Float **covar, ICOS_Float *chisq, void (*funcs)(ICOS_Float, ICOS_Float
[], int));
void linbcg(unsigned long n, double b[], double x[], int itol, double tol,
	 int itmax, int *iter, double *err);
void linmin(ICOS_Float p[], ICOS_Float xi[], int n, ICOS_Float *fret,
	ICOS_Float (*func)(ICOS_Float []));
void lnsrch(int n, ICOS_Float xold[], ICOS_Float fold, ICOS_Float g[], ICOS_Float p[], ICOS_Float
x[],
	 ICOS_Float *f, ICOS_Float stpmax, int *check, ICOS_Float (*func)(ICOS_Float []));
void load(ICOS_Float x1, ICOS_Float v[], ICOS_Float y[]);
void load1(ICOS_Float x1, ICOS_Float v1[], ICOS_Float y[]);
void load2(ICOS_Float x2, ICOS_Float v2[], ICOS_Float y[]);
void locate(ICOS_Float xx[], unsigned long n, ICOS_Float x, unsigned long *j);
void lop(double **out, double **u, int n);
void lubksb(ICOS_Float **a, int n, int *indx, ICOS_Float b[]);
void ludcmp(ICOS_Float **a, int n, int *indx, ICOS_Float *d);
void machar(int *ibeta, int *it, int *irnd, int *ngrd,
	int *machep, int *negep, int *iexp, int *minexp, int *maxexp,
	ICOS_Float *eps, ICOS_Float *epsneg, ICOS_Float *xmin, ICOS_Float *xmax);
void matadd(double **a, double **b, double **c, int n);
void matsub(double **a, double **b, double **c, int n);
void medfit(ICOS_Float x[], ICOS_Float y[], int ndata, ICOS_Float *a, ICOS_Float *b, ICOS_Float
*abdev);
void memcof(ICOS_Float data[], int n, int m, ICOS_Float *xms, ICOS_Float d[]);
int metrop(ICOS_Float de, ICOS_Float t);
void mgfas(double **u, int n, int maxcyc);

void mglin(double **u, int n, int ncycle);
ICOS_Float midexp(ICOS_Float (*funk)(ICOS_Float), ICOS_Float aa, ICOS_Float bb, int n);
ICOS_Float midinf(ICOS_Float (*funk)(ICOS_Float), ICOS_Float aa, ICOS_Float bb, int n);
ICOS_Float midpnt(ICOS_Float (*func)(ICOS_Float), ICOS_Float a, ICOS_Float b, int n);
ICOS_Float midsql(ICOS_Float (*funk)(ICOS_Float), ICOS_Float aa, ICOS_Float bb, int n);
ICOS_Float midsqu(ICOS_Float (*funk)(ICOS_Float), ICOS_Float aa, ICOS_Float bb, int n);
void miser(ICOS_Float (*func)(ICOS_Float []), ICOS_Float regn[], int ndim, unsigned long
npts,
	ICOS_Float dith, ICOS_Float *ave, ICOS_Float *var);
void mmid(ICOS_Float y[], ICOS_Float dydx[], int nvar, ICOS_Float xs, ICOS_Float htot,
	int nstep, ICOS_Float yout[], void (*derivs)(ICOS_Float, ICOS_Float[], ICOS_Float[]));
void mnbrak(ICOS_Float *ax, ICOS_Float *bx, ICOS_Float *cx, ICOS_Float *fa, ICOS_Float *fb,
	ICOS_Float *fc, ICOS_Float (*func)(ICOS_Float));
void mnewt(int ntrial, ICOS_Float x[], int n, ICOS_Float tolx, ICOS_Float tolf);
void moment(ICOS_Float data[], int n, ICOS_Float *ave, ICOS_Float *adev, ICOS_Float *sdev,
	ICOS_Float *var, ICOS_Float *skew, ICOS_Float *curt);
void mp2dfr(unsigned char a[], unsigned char s[], int n, int *m);
void mpadd(unsigned char w[], unsigned char u[], unsigned char v[], int
n);
void mpdiv(unsigned char q[], unsigned char r[], unsigned char u[],
	unsigned char v[], int n, int m);
void mpinv(unsigned char u[], unsigned char v[], int n, int m);
void mplsh(unsigned char u[], int n);
void mpmov(unsigned char u[], unsigned char v[], int n);
void mpmul(unsigned char w[], unsigned char u[], unsigned char v[], int n,
	int m);
void mpneg(unsigned char u[], int n);
void mppi(int n);
void mprove(ICOS_Float **a, ICOS_Float **alud, int n, int indx[], ICOS_Float b[],
	ICOS_Float x[]);
void mpsad(unsigned char w[], unsigned char u[], int n, int iv);
void mpsdv(unsigned char w[], unsigned char u[], int n, int iv, int *ir);
void mpsmu(unsigned char w[], unsigned char u[], int n, int iv);
void mpsqrt(unsigned char w[], unsigned char u[], unsigned char v[], int
n,
	int m);
void mpsub(int *is, unsigned char w[], unsigned char u[], unsigned char
v[],
	int n);
void mrqcofi(ICOS_Float x[], ICOS_Float y[], ICOS_Float sig[], int ndata, ICOS_Float a[],
	int ia[], int ma, ICOS_Float **alpha, ICOS_Float beta[], ICOS_Float *chisq,
	void (*funcs)(int, ICOS_Float [], ICOS_Float *, ICOS_Float [], int));
void mrqmini(ICOS_Float x[], ICOS_Float y[], ICOS_Float sig[], int ndata, ICOS_Float a[],
	int ia[], int ma, ICOS_Float **covar, ICOS_Float **alpha, ICOS_Float *chisq,
	void (*funcs)(int, ICOS_Float [], ICOS_Float *, ICOS_Float [], int), ICOS_Float
*alamda);
void mrqcof(ICOS_Float x[], ICOS_Float y[], ICOS_Float sig[], int ndata, ICOS_Float a[],
	int ia[], int ma, ICOS_Float **alpha, ICOS_Float beta[], ICOS_Float *chisq,
	void (*funcs)(ICOS_Float, ICOS_Float [], ICOS_Float *, ICOS_Float [], int));
void mrqmin(ICOS_Float x[], ICOS_Float y[], ICOS_Float sig[], int ndata, ICOS_Float a[],
	int ia[], int ma, ICOS_Float **covar, ICOS_Float **alpha, ICOS_Float *chisq,
	void (*funcs)(ICOS_Float, ICOS_Float [], ICOS_Float *, ICOS_Float [], int), ICOS_Float
*alamda);
void newt(ICOS_Float x[], int n, int *check,
	void (*vecfunc)(int, ICOS_Float [], ICOS_Float []));
void odeint(ICOS_Float ystart[], int nvar, ICOS_Float x1, ICOS_Float x2,
	ICOS_Float eps, ICOS_Float h1, ICOS_Float hmin, int *nok, int *nbad,
	void (*derivs)(ICOS_Float, ICOS_Float [], ICOS_Float []),
	void (*rkqs)(ICOS_Float [], ICOS_Float [], int, ICOS_Float *, ICOS_Float, ICOS_Float,
	ICOS_Float [], ICOS_Float *, ICOS_Float *, void (*)(ICOS_Float, ICOS_Float [], ICOS_Float [])));

void orthog(int n, ICOS_Float anu[], ICOS_Float alpha[], ICOS_Float beta[], ICOS_Float a[],
	ICOS_Float b[]);
void pade(double cof[], int n, ICOS_Float *resid);
void pccheb(ICOS_Float d[], ICOS_Float c[], int n);
void pcshft(ICOS_Float a, ICOS_Float b, ICOS_Float d[], int n);
void pearsn(ICOS_Float x[], ICOS_Float y[], unsigned long n, ICOS_Float *r, ICOS_Float *prob,
	ICOS_Float *z);
void period(ICOS_Float x[], ICOS_Float y[], int n, ICOS_Float ofac, ICOS_Float hifac,
	ICOS_Float px[], ICOS_Float py[], int np, int *nout, int *jmax, ICOS_Float
*prob);
void piksr2(int n, ICOS_Float arr[], ICOS_Float brr[]);
void piksrt(int n, ICOS_Float arr[]);
void pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k,
	ICOS_Float ***c, ICOS_Float **s);
ICOS_Float plgndr(int l, int m, ICOS_Float x);
ICOS_Float poidev(ICOS_Float xm, long *idum);
void polcoe(ICOS_Float x[], ICOS_Float y[], int n, ICOS_Float cof[]);
void polcof(ICOS_Float xa[], ICOS_Float ya[], int n, ICOS_Float cof[]);
void poldiv(ICOS_Float u[], int n, ICOS_Float v[], int nv, ICOS_Float q[], ICOS_Float r[]);
void polin2(ICOS_Float x1a[], ICOS_Float x2a[], ICOS_Float **ya, int m, int n,
	ICOS_Float x1, ICOS_Float x2, ICOS_Float *y, ICOS_Float *dy);
void polint(ICOS_Float xa[], ICOS_Float ya[], int n, ICOS_Float x, ICOS_Float *y, ICOS_Float *dy);
void powell(ICOS_Float p[], ICOS_Float **xi, int n, ICOS_Float ftol, int *iter, ICOS_Float
*fret,
	ICOS_Float (*func)(ICOS_Float []));
void predic(ICOS_Float data[], int ndata, ICOS_Float d[], int m, ICOS_Float future[], int
nfut);
ICOS_Float probks(ICOS_Float alam);
void psdes(unsigned long *lword, unsigned long *irword);
void pwt(ICOS_Float a[], unsigned long n, int isign);

void pwtset(int n);
ICOS_Float pythag(ICOS_Float a, ICOS_Float b);
void pzextr(int iest, ICOS_Float xest, ICOS_Float yest[], ICOS_Float yz[], ICOS_Float dy[],
	int nv);
ICOS_Float qgaus(ICOS_Float (*func)(ICOS_Float), ICOS_Float a, ICOS_Float b);
void qrdcmp(ICOS_Float **a, int n, ICOS_Float *c, ICOS_Float *d, int *sing);
ICOS_Float qromb(ICOS_Float (*func)(ICOS_Float), ICOS_Float a, ICOS_Float b);
ICOS_Float qromo(ICOS_Float (*func)(ICOS_Float), ICOS_Float a, ICOS_Float b,
	ICOS_Float (*choose)(ICOS_Float (*)(ICOS_Float), ICOS_Float, ICOS_Float, int));
void qroot(ICOS_Float p[], int n, ICOS_Float *b, ICOS_Float *c, ICOS_Float eps);
void qrsolv(ICOS_Float **a, int n, ICOS_Float c[], ICOS_Float d[], ICOS_Float b[]);
void qrupdt(ICOS_Float **r, ICOS_Float **qt, int n, ICOS_Float u[], ICOS_Float v[]);
ICOS_Float qsimp(ICOS_Float (*func)(ICOS_Float), ICOS_Float a, ICOS_Float b);
ICOS_Float qtrap(ICOS_Float (*func)(ICOS_Float), ICOS_Float a, ICOS_Float b);
ICOS_Float quad3d(ICOS_Float (*func)(ICOS_Float, ICOS_Float, ICOS_Float), ICOS_Float x1, ICOS_Float x2);
void quadct(ICOS_Float x, ICOS_Float y, ICOS_Float xx[], ICOS_Float yy[], unsigned long nn,
	ICOS_Float *fa, ICOS_Float *fb, ICOS_Float *fc, ICOS_Float *fd);
void quadmx(ICOS_Float **a, int n);
void quadvl(ICOS_Float x, ICOS_Float y, ICOS_Float *fa, ICOS_Float *fb, ICOS_Float *fc, ICOS_Float *fd);
ICOS_Float ran0(long *idum);
ICOS_Float ran1(long *idum);
ICOS_Float ran2(long *idum);
ICOS_Float ran3(long *idum);
ICOS_Float ran4(long *idum);
void rank(unsigned long n, unsigned long indx[], unsigned long irank[]);
void ranpt(ICOS_Float pt[], ICOS_Float regn[], int n);
void ratint(ICOS_Float xa[], ICOS_Float ya[], int n, ICOS_Float x, ICOS_Float *y, ICOS_Float *dy);
void ratlsq(double (*fn)(double), double a, double b, int mm, int kk,
	double cof[], double *dev);
double ratval(double x, double cof[], int mm, int kk);
ICOS_Float rc(ICOS_Float x, ICOS_Float y);
ICOS_Float rd(ICOS_Float x, ICOS_Float y, ICOS_Float z);
void realft(ICOS_Float data[], unsigned long n, int isign);
void rebin(ICOS_Float rc, int nd, ICOS_Float r[], ICOS_Float xin[], ICOS_Float xi[]);
void red(int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, int jmf,
	int ic1, int jc1, int jcf, int kc, ICOS_Float ***c, ICOS_Float **s);
void relax(double **u, double **rhs, int n);
void relax2(double **u, double **rhs, int n);
void resid(double **res, double **u, double **rhs, int n);
ICOS_Float revcst(ICOS_Float x[], ICOS_Float y[], int iorder[], int ncity, int n[]);
void reverse(int iorder[], int ncity, int n[]);
ICOS_Float rf(ICOS_Float x, ICOS_Float y, ICOS_Float z);
ICOS_Float rj(ICOS_Float x, ICOS_Float y, ICOS_Float z, ICOS_Float p);
void rk4(ICOS_Float y[], ICOS_Float dydx[], int n, ICOS_Float x, ICOS_Float h, ICOS_Float yout[],
	void (*derivs)(ICOS_Float, ICOS_Float [], ICOS_Float []));
void rkck(ICOS_Float y[], ICOS_Float dydx[], int n, ICOS_Float x, ICOS_Float h,
	ICOS_Float yout[], ICOS_Float yerr[], void (*derivs)(ICOS_Float, ICOS_Float [], ICOS_Float
[]));
void rkdumb(ICOS_Float vstart[], int nvar, ICOS_Float x1, ICOS_Float x2, int nstep,

void (*derivs)(ICOS_Float, ICOS_Float [], ICOS_Float []));
void rkqs(ICOS_Float y[], ICOS_Float dydx[], int n, ICOS_Float *x,
	ICOS_Float htry, ICOS_Float eps, ICOS_Float yscal[], ICOS_Float *hdid, ICOS_Float *hnext,
	void (*derivs)(ICOS_Float, ICOS_Float [], ICOS_Float []));
void rlft3(ICOS_Float ***data, ICOS_Float **speq, unsigned long nn1,
	unsigned long nn2, unsigned long nn3, int isign);
ICOS_Float rofunc(ICOS_Float b);
void rotate(ICOS_Float **r, ICOS_Float **qt, int n, int i, ICOS_Float a, ICOS_Float b);
void rsolv(ICOS_Float **a, int n, ICOS_Float d[], ICOS_Float b[]);
void rstrct(double **uc, double **uf, int nc);
ICOS_Float rtbis(ICOS_Float (*func)(ICOS_Float), ICOS_Float x1, ICOS_Float x2, ICOS_Float xacc);
ICOS_Float rtflsp(ICOS_Float (*func)(ICOS_Float), ICOS_Float x1, ICOS_Float x2, ICOS_Float xacc);
ICOS_Float rtnewt(void (*funcd)(ICOS_Float, ICOS_Float *, ICOS_Float *), ICOS_Float x1, ICOS_Float x2,
	ICOS_Float xacc);
ICOS_Float rtsafe(void (*funcd)(ICOS_Float, ICOS_Float *, ICOS_Float *), ICOS_Float x1, ICOS_Float x2,
	ICOS_Float xacc);
ICOS_Float rtsec(ICOS_Float (*func)(ICOS_Float), ICOS_Float x1, ICOS_Float x2, ICOS_Float xacc);
void rzextr(int iest, ICOS_Float xest, ICOS_Float yest[], ICOS_Float yz[], ICOS_Float dy[],
int nv);
void savgol(ICOS_Float c[], int np, int nl, int nr, int ld, int m);
void score(ICOS_Float xf, ICOS_Float y[], ICOS_Float f[]);
void scrsho(ICOS_Float (*fx)(ICOS_Float));
ICOS_Float select(unsigned long k, unsigned long n, ICOS_Float arr[]);
ICOS_Float selip(unsigned long k, unsigned long n, ICOS_Float arr[]);
void shell(unsigned long n, ICOS_Float a[]);
void shoot(int n, ICOS_Float v[], ICOS_Float f[]);
void shootf(int n, ICOS_Float v[], ICOS_Float f[]);
void simp1(ICOS_Float **a, int mm, int ll[], int nll, int iabf, int *kp,
	ICOS_Float *bmax);
void simp2(ICOS_Float **a, int m, int n, int *ip, int kp);
void simp3(ICOS_Float **a, int i1, int k1, int ip, int kp);
void simplx(ICOS_Float **a, int m, int n, int m1, int m2, int m3, int *icase,
	int izrov[], int iposv[]);
void simpr(ICOS_Float y[], ICOS_Float dydx[], ICOS_Float dfdx[], ICOS_Float **dfdy,
	int n, ICOS_Float xs, ICOS_Float htot, int nstep, ICOS_Float yout[],
	void (*derivs)(ICOS_Float, ICOS_Float [], ICOS_Float []));
void sinft(ICOS_Float y[], int n);
void slvsm2(double **u, double **rhs);
void slvsml(double **u, double **rhs);
void sncndn(ICOS_Float uu, ICOS_Float emmc, ICOS_Float *sn, ICOS_Float *cn, ICOS_Float *dn);
double snrm(unsigned long n, double sx[], int itol);
void sobseq(int *n, ICOS_Float x[]);

void solvde(int itmax, ICOS_Float conv, ICOS_Float slowc, ICOS_Float scalv[],
	int indexv[], int ne, int nb, int m, ICOS_Float **y, ICOS_Float ***c, ICOS_Float
**s);
void sor(double **a, double **b, double **c, double **d, double **e,
	double **f, double **u, int jmax, double rjac);
void sort(unsigned long n, ICOS_Float arr[]);
void sort2(unsigned long n, ICOS_Float arr[], ICOS_Float brr[]);
void sort3(unsigned long n, ICOS_Float ra[], ICOS_Float rb[], ICOS_Float rc[]);
void spctrm(FILE *fp, ICOS_Float p[], int m, int k, int ovrlap);
void spear(ICOS_Float data1[], ICOS_Float data2[], unsigned long n, ICOS_Float *d, ICOS_Float
*zd,
	ICOS_Float *probd, ICOS_Float *rs, ICOS_Float *probrs);
void sphbes(int n, ICOS_Float x, ICOS_Float *sj, ICOS_Float *sy, ICOS_Float *sjp, ICOS_Float *syp);
void splie2(ICOS_Float x1a[], ICOS_Float x2a[], ICOS_Float **ya, int m, int n, ICOS_Float
**y2a);
void splin2(ICOS_Float x1a[], ICOS_Float x2a[], ICOS_Float **ya, ICOS_Float **y2a, int m, int
n,
	ICOS_Float x1, ICOS_Float x2, ICOS_Float *y);
void spline(ICOS_Float x[], ICOS_Float y[], int n, ICOS_Float yp1, ICOS_Float ypn, ICOS_Float
y2[]);
void splint(ICOS_Float xa[], ICOS_Float ya[], ICOS_Float y2a[], int n, ICOS_Float x, ICOS_Float
*y);
void spread(ICOS_Float y, ICOS_Float yy[], unsigned long n, ICOS_Float x, int m);
void sprsax(ICOS_Float sa[], unsigned long ija[], ICOS_Float x[], ICOS_Float b[],
	unsigned long n);
void sprsin(ICOS_Float **a, int n, ICOS_Float thresh, unsigned long nmax, ICOS_Float
sa[],
	unsigned long ija[]);
void sprspm(ICOS_Float sa[], unsigned long ija[], ICOS_Float sb[], unsigned long
ijb[],
	ICOS_Float sc[], unsigned long ijc[]);
void sprstm(ICOS_Float sa[], unsigned long ija[], ICOS_Float sb[], unsigned long
ijb[],
	ICOS_Float thresh, unsigned long nmax, ICOS_Float sc[], unsigned long
ijc[]);
void sprstp(ICOS_Float sa[], unsigned long ija[], ICOS_Float sb[], unsigned long
ijb[]);
void sprstx(ICOS_Float sa[], unsigned long ija[], ICOS_Float x[], ICOS_Float b[],
	unsigned long n);
void stifbs(ICOS_Float y[], ICOS_Float dydx[], int nv, ICOS_Float *xx,
	ICOS_Float htry, ICOS_Float eps, ICOS_Float yscal[], ICOS_Float *hdid, ICOS_Float *hnext,
	void (*derivs)(ICOS_Float, ICOS_Float [], ICOS_Float []));
void stiff(ICOS_Float y[], ICOS_Float dydx[], int n, ICOS_Float *x,
	ICOS_Float htry, ICOS_Float eps, ICOS_Float yscal[], ICOS_Float *hdid, ICOS_Float *hnext,
	void (*derivs)(ICOS_Float, ICOS_Float [], ICOS_Float []));
void stoerm(ICOS_Float y[], ICOS_Float d2y[], int nv, ICOS_Float xs,
	ICOS_Float htot, int nstep, ICOS_Float yout[],
	void (*derivs)(ICOS_Float, ICOS_Float [], ICOS_Float []));
void svbksb(ICOS_Float **u, ICOS_Float w[], ICOS_Float **v, int m, int n, ICOS_Float b[],
	ICOS_Float x[]);
void svdcmp(ICOS_Float **a, int m, int n, ICOS_Float w[], ICOS_Float **v);
void svdfit(ICOS_Float x[], ICOS_Float y[], ICOS_Float sig[], int ndata, ICOS_Float a[],
	int ma, ICOS_Float **u, ICOS_Float **v, ICOS_Float w[], ICOS_Float *chisq,
	void (*funcs)(ICOS_Float, ICOS_Float [], int));
void svdvar(ICOS_Float **v, int ma, ICOS_Float w[], ICOS_Float **cvm);
void toeplz(ICOS_Float r[], ICOS_Float x[], ICOS_Float y[], int n);
void tptest(ICOS_Float data1[], ICOS_Float data2[], unsigned long n, ICOS_Float *t, ICOS_Float
*prob);
void tqli(ICOS_Float d[], ICOS_Float e[], int n, ICOS_Float **z);

ICOS_Float trapzd(ICOS_Float (*func)(ICOS_Float), ICOS_Float a, ICOS_Float b, int n);
void tred2(ICOS_Float **a, int n, ICOS_Float d[], ICOS_Float e[]);
void tridag(ICOS_Float a[], ICOS_Float b[], ICOS_Float c[], ICOS_Float r[], ICOS_Float u[],
	unsigned long n);
ICOS_Float trncst(ICOS_Float x[], ICOS_Float y[], int iorder[], int ncity, int n[]);
void trnspt(int iorder[], int ncity, int n[]);
void ttest(ICOS_Float data1[], unsigned long n1, ICOS_Float data2[], unsigned long
n2,
	ICOS_Float *t, ICOS_Float *prob);
void tutest(ICOS_Float data1[], unsigned long n1, ICOS_Float data2[], unsigned long
n2,
	ICOS_Float *t, ICOS_Float *prob);
void twofft(ICOS_Float data1[], ICOS_Float data2[], ICOS_Float fft1[], ICOS_Float fft2[],
	unsigned long n);
void vander(double x[], double w[], double q[], int n);
void vegas(ICOS_Float regn[], int ndim, ICOS_Float (*fxn)(ICOS_Float [], ICOS_Float), int
init,
	unsigned long ncall, int itmx, int nprn, ICOS_Float *tgral, ICOS_Float *sd,
	ICOS_Float *chi2a);
void voltra(int n, int m, ICOS_Float t0, ICOS_Float h, ICOS_Float *t, ICOS_Float **f,
	ICOS_Float (*g)(int, ICOS_Float), ICOS_Float (*ak)(int, int, ICOS_Float, ICOS_Float));
void wt1(ICOS_Float a[], unsigned long n, int isign,
	void (*wtstep)(ICOS_Float [], unsigned long, int));
void wtn(ICOS_Float a[], unsigned long nn[], int ndim, int isign,
	void (*wtstep)(ICOS_Float [], unsigned long, int));
void wwghts(ICOS_Float wghts[], int n, ICOS_Float h,
	void (*kermom)(double [], double ,int));
int zbrac(ICOS_Float (*func)(ICOS_Float), ICOS_Float *x1, ICOS_Float *x2);
void zbrak(ICOS_Float (*fx)(ICOS_Float), ICOS_Float x1, ICOS_Float x2, int n, ICOS_Float xb1[],
	ICOS_Float xb2[], int *nb);
ICOS_Float zbrent(ICOS_Float (*func)(ICOS_Float), ICOS_Float x1, ICOS_Float x2, ICOS_Float tol);
void zrhqr(ICOS_Float a[], int m, ICOS_Float rtr[], ICOS_Float rti[]);
ICOS_Float zriddr(ICOS_Float (*func)(ICOS_Float), ICOS_Float x1, ICOS_Float x2, ICOS_Float xacc);
void zroots(fcomplex a[], int m, fcomplex roots[], int polish);
   
#endif /* _NR_H_ */
