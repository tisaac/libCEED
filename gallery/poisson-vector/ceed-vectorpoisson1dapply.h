// Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory. LLNL-CODE-734707.
// All Rights reserved. See files LICENSE and NOTICE for details.
//
// This file is part of CEED, a collection of benchmarks, miniapps, software
// libraries and APIs for efficient high-order finite element and spectral
// element discretizations for exascale applications. For more information and
// source code availability see http://github.com/ceed.
//
// The CEED research is supported by the Exascale Computing Project 17-SC-20-SC,
// a collaborative effort of two U.S. Department of Energy organizations (Office
// of Science and the National Nuclear Security Administration) responsible for
// the planning and preparation of a capable exascale ecosystem, including
// software, applications, hardware, advanced system engineering and early
// testbed platforms, in support of the nation's exascale computing imperative.

/**
  @brief Ceed QFunction for applying the 1D Poisson operator
           on a vector system with three components
**/

#ifndef vectorpoisson1dapply_h
#define vectorpoisson1dapply_h

CEED_QFUNCTION(Vector3Poisson1DApply)(void *ctx, const CeedInt Q,
                                      const CeedScalar *const *in,
                                      CeedScalar *const *out) {
  // *INDENT-OFF*
  // in[0] is gradient u, shape [1, nc=3, Q]
  // in[1] is quadrature data, size (Q)
  const CeedScalar (*ug)[CEED_Q_VLA] = (const CeedScalar(*)[CEED_Q_VLA])in[0],
                           (*q_data) = in[1];
  // out[0] is output to multiply against gradient v, shape [1, nc=3, Q]
  CeedScalar       (*vg)[CEED_Q_VLA] = (CeedScalar(*)[CEED_Q_VLA])out[0];
  // *INDENT-ON*

  const CeedInt num_comp = 3;

  // Quadrature point loop
  CeedPragmaSIMD
  for (CeedInt i=0; i<Q; i++) {
    for (CeedInt c=0; c<num_comp; c++) {
      vg[c][i] = ug[c][i] * q_data[i];
    }
  } // End of Quadrature Point Loop

  return CEED_ERROR_SUCCESS;
}

#endif // vectorpoisson1dapply_h