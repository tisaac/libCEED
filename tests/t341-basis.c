/// @file
/// Test interpolation of P0 L2 basis
/// \test Test interpolaton of P0 L2 basis
#include <ceed.h>
#include <math.h>
#include "t340-basis.h"

int main(int argc, char **argv) {
  Ceed ceed;
  const CeedInt num_nodes = 1, Q = 2, dim = 2, num_qpts = Q*Q;
  CeedInt num_comp = 1;
  CeedInt P = num_nodes;
  CeedBasis b;
  CeedScalar q_ref[dim*num_qpts], q_weights[num_qpts];
  CeedScalar interp[P*num_qpts];
  CeedVector X, Y, U;
  const CeedScalar *y, *u;

  CeedInit(argv[1], &ceed);

  buildmats(Q, q_ref, q_weights, interp, CEED_GAUSS);
  CeedBasisCreateL2(ceed, CEED_QUAD, num_comp, num_nodes, num_qpts, interp,
                    q_ref, q_weights, &b);
  // Test GetDiv
  const CeedScalar *interp2;
  // Test GetInterp for Hdiv
  CeedBasisGetInterp(b, &interp2);
  for (CeedInt i=0; i<P*num_qpts; i++) {
    if (fabs(interp[i] - interp2[i]) > 100.*CEED_EPSILON)
      // LCOV_EXCL_START
      printf("%f != %f\n", interp[i], interp2[i]);
    // LCOV_EXCL_STOP
  }

  CeedVectorCreate(ceed, P, &X);
  CeedVectorSetValue(X, 1);
  CeedVectorCreate(ceed, num_qpts, &Y);
  CeedVectorSetValue(Y, 0);
  CeedBasisApply(b, 1, CEED_NOTRANSPOSE, CEED_EVAL_INTERP, X, Y);

  CeedVectorCreate(ceed, P, &U);
  CeedVectorSetValue(U, 0.0);
  CeedBasisApply(b, 1, CEED_TRANSPOSE, CEED_EVAL_INTERP, Y, U);

  CeedVectorGetArrayRead(Y, CEED_MEM_HOST, &y);
  CeedVectorGetArrayRead(U, CEED_MEM_HOST, &u);
  // Check CEED_NOTRANSPOSE
  for (CeedInt i=0; i<num_qpts; i++) {
    if (fabs(1. - y[i]) > 100.*CEED_EPSILON)
      // LCOV_EXCL_START
      printf("%f != %f\n", 1.0, y[i]);
    // LCOV_EXCL_STOP
  }

  // Check CEED_TRANSPOSE
  for (CeedInt i=0; i<P; i++) {
    if (fabs(num_qpts*1. - u[i]) > 100.*CEED_EPSILON)
      // LCOV_EXCL_START
      printf("%f != %f\n", num_qpts*1., u[i]);
    // LCOV_EXCL_STOP
  }
  CeedVectorRestoreArrayRead(Y, &y);
  CeedVectorRestoreArrayRead(U, &u);

  CeedBasisDestroy(&b);
  CeedVectorDestroy(&X);
  CeedVectorDestroy(&Y);
  CeedVectorDestroy(&U);
  CeedDestroy(&ceed);
  return 0;
}
