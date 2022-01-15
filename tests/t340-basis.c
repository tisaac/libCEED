/// @file
/// Test creation and destruction of P0 L2 basis
/// \test Test creation and distruction of P0 L2 basis
#include <ceed.h>
#include "t340-basis.h"

int main(int argc, char **argv) {
  Ceed ceed;
  const CeedInt Q = 2, dim = 2, num_qpts = Q*Q, elem_nodes = 1;
  CeedInt num_comp = 1;
  CeedInt P = elem_nodes;
  CeedBasis b;
  CeedScalar q_ref[dim*num_qpts], q_weights[num_qpts];
  CeedScalar interp[P*num_qpts];

  CeedInit(argv[1], &ceed);

  // Test skipped if using single precision
  if (CEED_SCALAR_TYPE == CEED_SCALAR_FP32) {
    return CeedError(ceed, CEED_ERROR_UNSUPPORTED,
                     "Test not implemented in single precision");
  }

  buildmats(Q, q_ref, q_weights, interp, CEED_GAUSS);
  CeedBasisCreateL2(ceed, CEED_QUAD, num_comp, elem_nodes, num_qpts, interp,
                    q_ref,q_weights, &b);
  CeedBasisView(b, stdout);

  CeedBasisDestroy(&b);
  CeedDestroy(&ceed);
  return 0;
}
