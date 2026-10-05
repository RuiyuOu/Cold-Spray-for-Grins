// goh_plugin.cpp -- FreeFEM plugin exposing the C++ GOH kernel inside variational forms.
//   gohSetParams(c, k1, k2, kappa, gamma)          set material parameters (gamma in radians); returns 0.
//   gohPsi(F11,F12,F13,F21,F22,F23,F31,F32,F33, p)         -> energy density
//   gohP  (i, F11,...,F33, p)   i = 0..8,  i = 3*I + J     -> P_IJ
//   gohA  (r, F11,...,F33, p)   r = 0..80, r = 9*(3*i+J) + (3*k+L)  -> dP_iJ/dF_kL
// The kernel is evaluated once per distinct (F,p) and cached, so the 9 or 81
// calls made at one quadrature point share a single kernel evaluation.
#include "ff++.hpp"
// Kernel selection (both implement fibres in the y-z plane, a = (0, cos g, +/- sin g)):
//   default           : goh.hpp      (Mathematica-generated kernel, function goh)
//   -DGOH_USE_REF     : goh_ref.hpp  (hand-written reference kernel, function goh_ref)
#ifdef GOH_USE_REF
#include "goh_ref.hpp"
#define GOH_KERNEL goh_ref
#else
#include "goh.hpp"
#define GOH_KERNEL goh
#endif
#include <cstring>

using namespace Fem2D;

static double gc = 7.64e3, gk1 = 996.6e3, gk2 = 524.6, gkap = 0.226, ggam = 0.8724;

struct Cache {
  bool valid = false;
  double in[10];
  double psi, P[9], A[81];
};
static Cache cache;

static void evalKernel(const double* F, double p) {
  if (cache.valid && std::memcmp(cache.in, F, 9 * sizeof(double)) == 0 && cache.in[9] == p) return;
  GOH_KERNEL(F, p, gc, gk1, gk2, gkap, ggam, &cache.psi, cache.P, cache.A);
  std::memcpy(cache.in, F, 9 * sizeof(double));
  cache.in[9] = p;
  cache.valid = true;
}

static void setParams(double c, double k1, double k2, double kap, double gam) {
  gc = c; gk1 = k1; gk2 = k2; gkap = kap; ggam = gam;
  cache.valid = false;
}

// ---- expression node: mode 0 = psi, 1 = P component, 2 = A component, 3 = set parameters
class GohNode : public E_F0mps {
 public:
  int mode;
  Expression sel;
  Expression a[10];
  GohNode(int m, const basicAC_F0& args) : mode(m), sel(0) {
    int n = 0;
    if (mode == 3) { for (int i = 0; i < 5; i++) a[i] = CastTo<double>(args[n++]); return; }
    if (mode > 0) sel = CastTo<long>(args[n++]);
    for (int i = 0; i < 10; i++) a[i] = CastTo<double>(args[n++]);
  }
  AnyType operator()(Stack s) const {
    if (mode == 3) {
      setParams(GetAny<double>((*a[0])(s)), GetAny<double>((*a[1])(s)), GetAny<double>((*a[2])(s)),
                GetAny<double>((*a[3])(s)), GetAny<double>((*a[4])(s)));
      return SetAny<double>(0.);
    }
    double F[9], p;
    for (int i = 0; i < 9; i++) F[i] = GetAny<double>((*a[i])(s));
    p = GetAny<double>((*a[9])(s));
    evalKernel(F, p);
    double v;
    if (mode == 0) v = cache.psi;
    else {
      long k = GetAny<long>((*sel)(s));
      if (mode == 1) v = (k >= 0 && k < 9)  ? cache.P[k] : 0.;
      else           v = (k >= 0 && k < 81) ? cache.A[k] : 0.;
    }
    return SetAny<double>(v);
  }
  operator aType() const { return atype<double>(); }
};

class GohOp : public OneOperator {
  int mode;
 public:
  E_F0* code(const basicAC_F0& args) const { return new GohNode(mode, args); }
  GohOp(int m) : OneOperator(atype<double>()), mode(m) {
    const int na = (m == 3) ? 5 : ((m > 0) ? 11 : 10);
    n = na;
    t = new aType[na];
    int k = 0;
    if (m == 1 || m == 2) t[k++] = atype<long>();
    for (; k < na; k++) t[k] = atype<double>();
  }
};

static void Load_Init() {
  Global.Add("gohSetParams", "(", new GohOp(3));
  Global.Add("gohPsi", "(", new GohOp(0));
  Global.Add("gohP",   "(", new GohOp(1));
  Global.Add("gohA",   "(", new GohOp(2));
}
LOADFUNC(Load_Init)
