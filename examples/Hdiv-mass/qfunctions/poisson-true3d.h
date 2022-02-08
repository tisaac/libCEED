// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
// reserved. See files LICENSE and NOTICE for details.
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

/// @file
/// Compute true solution of the H(div) example using PETSc

#ifndef TRUE3D_H
#define TRUE3D_H

#include <math.h>

// -----------------------------------------------------------------------------
// Compuet true solution
// -----------------------------------------------------------------------------
CEED_QFUNCTION(SetupTrueSoln3D)(void *ctx, const CeedInt Q,
                                const CeedScalar *const *in,
                                CeedScalar *const *out) {
  // *INDENT-OFF*
  // Inputs
  const CeedScalar (*coords) = in[0],
                   (*dxdX)[3][CEED_Q_VLA] = (const CeedScalar(*)[3][CEED_Q_VLA])in[1];
  // Outputs
  CeedScalar (*tr_Hdiv) = out[0];

  CeedScalar x0 = coords[0+0*Q], y0 = coords[0+1*Q], z0 = coords[0+2*Q];
  CeedScalar x1 = coords[1+0*Q], y1 = coords[1+1*Q], z1 = coords[1+2*Q];
  CeedScalar x2 = coords[2+0*Q], y2 = coords[2+1*Q], z2 = coords[2+2*Q];
  CeedScalar x3 = coords[3+0*Q], y3 = coords[3+1*Q], z3 = coords[3+2*Q];
  CeedScalar x4 = coords[4+0*Q], y4 = coords[4+1*Q], z4 = coords[4+2*Q];
  CeedScalar x5 = coords[5+0*Q], y5 = coords[5+1*Q], z5 = coords[5+2*Q];
  CeedScalar x6 = coords[6+0*Q], y6 = coords[6+1*Q], z6 = coords[6+2*Q];
  CeedScalar x7 = coords[7+0*Q], y7 = coords[7+1*Q], z7 = coords[7+2*Q];


  CeedScalar ue0[3] = {x0,y0,z0};//{x0+y0+z0, x0-y0+z0, x0+y0-z0};
  CeedScalar ue1[3] = {x1,y1,z1};//{x1+y1+z1, x1-y1+z1, x1+y1-z1};
  CeedScalar ue2[3] = {x2,y2,z2};//{x2+y2+z2, x2-y2+z2, x2+y2-z2};
  CeedScalar ue3[3] = {x3,y3,z3};//{x3+y3+z3, x3-y3+z3, x3+y3-z3};
  CeedScalar ue4[3] = {x4,y4,z4};//{x4+y4+z4, x4-y4+z4, x4+y4-z4};
  CeedScalar ue5[3] = {x5,y5,z5};//{x5+y5+z5, x5-y5+z5, x5+y5-z5};
  CeedScalar ue6[3] = {x6,y6,z6};//{x6+y6+z6, x6-y6+z6, x6+y6-z6};
  CeedScalar ue7[3] = {x7,y7,z7};//{x7+y7+z7, x7-y7+z7, x7+y7-z7};
  CeedScalar nl[3] = {0., -0.125, 0.};
  CeedScalar nr[3] = {0., 0.125, 0.};
  CeedScalar nbt[3] = {0., 0., -0.125};
  CeedScalar nt[3] = {0., 0., 0.125};
  CeedScalar nf[3] = {0.25, 0., 0.};
  CeedScalar nbk[3] = {-0.25, 0., 0.};
  CeedScalar d0, d1, d2, d3, d4, d5, d6, d7;
  CeedScalar d8, d9, d10, d11, d12, d13, d14, d15;
  CeedScalar d16, d17, d18, d19, d20, d21, d22, d23;

  d0 = ue0[0]*nbt[0]+ue0[1]*nbt[1]+ue0[2]*nbt[2];
  d1 = ue1[0]*nbt[0]+ue1[1]*nbt[1]+ue1[2]*nbt[2];
  d2 = ue1[0]*nbt[0]+ue1[1]*nbt[1]+ue1[2]*nbt[2];
  d3 = ue3[0]*nbt[0]+ue3[1]*nbt[1]+ue3[2]*nbt[2];

  d4 = ue3[0]*nr[0]+ue3[1]*nr[1]+ue3[2]*nr[2];
  d5 = ue2[0]*nr[0]+ue2[1]*nr[1]+ue2[2]*nr[2];
  d6 = ue7[0]*nr[0]+ue7[1]*nr[1]+ue7[2]*nr[2];
  d7 = ue6[0]*nr[0]+ue6[1]*nr[1]+ue6[2]*nr[2];

  d8 = ue4[0]*nt[0]+ue4[1]*nt[1]+ue4[2]*nt[2];
  d9 = ue5[0]*nt[0]+ue5[1]*nt[1]+ue5[2]*nt[2];
  d10 = ue6[0]*nt[0]+ue6[1]*nt[1]+ue6[2]*nt[2];
  d11 = ue7[0]*nt[0]+ue7[1]*nt[1]+ue7[2]*nt[2];

  d12 = ue1[0]*nl[0]+ue1[1]*nl[1]+ue1[2]*nl[2];
  d13 = ue0[0]*nl[0]+ue0[1]*nl[1]+ue0[2]*nl[2];
  d14 = ue5[0]*nl[0]+ue5[1]*nl[1]+ue5[2]*nl[2];
  d15 = ue4[0]*nl[0]+ue4[1]*nl[1]+ue4[2]*nl[2];

  d16 = ue1[0]*nf[0]+ue1[1]*nf[1]+ue1[2]*nf[2];
  d17 = ue3[0]*nf[0]+ue3[1]*nf[1]+ue3[2]*nf[2];
  d18 = ue5[0]*nf[0]+ue5[1]*nf[1]+ue5[2]*nf[2];
  d19 = ue7[0]*nf[0]+ue7[1]*nf[1]+ue7[2]*nf[2];

  d20 = ue0[0]*nbk[0]+ue0[1]*nbk[1]+ue0[2]*nbk[2];
  d21 = ue2[0]*nbk[0]+ue2[1]*nbk[1]+ue2[2]*nbk[2];
  d22 = ue4[0]*nbk[0]+ue4[1]*nbk[1]+ue4[2]*nbk[2];
  d23 = ue6[0]*nbk[0]+ue6[1]*nbk[1]+ue6[2]*nbk[2];

  // bottom
  tr_Hdiv[0] = d0;
  tr_Hdiv[1] = d1;
  tr_Hdiv[2] = d2;
  tr_Hdiv[3] = d3;
  // right
  tr_Hdiv[4] = d4;
  tr_Hdiv[5] = d5;
  tr_Hdiv[6] = d6;
  tr_Hdiv[7] = d7;
  // top
  tr_Hdiv[8] = d8;
  tr_Hdiv[9] = d9;
  tr_Hdiv[10] = d10;
  tr_Hdiv[11] = d11;
  // left
  tr_Hdiv[12] = d12;
  tr_Hdiv[13] = d13;
  tr_Hdiv[14] = d14;
  tr_Hdiv[15] = d15;
  // front
  tr_Hdiv[16] = d16;
  tr_Hdiv[17] = d17;
  tr_Hdiv[18] = d18;
  tr_Hdiv[19] = d19;
  // back
  tr_Hdiv[20] = d20;
  tr_Hdiv[21] = d21;
  tr_Hdiv[22] = d22;
  tr_Hdiv[23] = d23;

  //printf("bottom: %f, %f, %f, %f\n", tr_Hdiv[0], tr_Hdiv[1], tr_Hdiv[2], tr_Hdiv[3]);
  //printf("right: %f, %f, %f, %f\n", tr_Hdiv[4], tr_Hdiv[5], tr_Hdiv[6], tr_Hdiv[7]);
  //printf("top: %f, %f, %f, %f\n", tr_Hdiv[8], tr_Hdiv[9], tr_Hdiv[10], tr_Hdiv[11]);
  //printf("left: %f, %f, %f, %f\n", tr_Hdiv[12], tr_Hdiv[13], tr_Hdiv[14], tr_Hdiv[15]);
  //printf("front: %f, %f, %f, %f\n", tr_Hdiv[16], tr_Hdiv[17], tr_Hdiv[18], tr_Hdiv[19]);
  //printf("back: %f, %f, %f, %f\n", tr_Hdiv[20], tr_Hdiv[21], tr_Hdiv[22], tr_Hdiv[23]);

  return 0;
}
// -----------------------------------------------------------------------------

#endif // End TRUE_H
