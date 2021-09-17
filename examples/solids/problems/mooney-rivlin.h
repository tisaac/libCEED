#ifndef mooney_rivlin_h
#define mooney_rivlin_h

#include <petsc.h>
#include "../include/structs.h"

#ifndef PHYSICS_STRUCT_MR
#define PHYSICS_STRUCT_MR
typedef struct Physics_MR_ *Physics_MR;
struct Physics_MR_ {
  // Material properties for MR
  CeedScalar mu_1;
  CeedScalar mu_2;
  CeedScalar lambda;
};
#endif // PHYSICS_STRUCT_MR

// Create context object
PetscErrorCode PhysicsContext_MR(MPI_Comm comm, Ceed ceed, Units *units,
                                 CeedQFunctionContext *ctx);
PetscErrorCode PhysicsSmootherContext_MR(MPI_Comm comm, Ceed ceed,
    CeedQFunctionContext ctx, CeedQFunctionContext *ctx_smoother);

// Process physics options - Mooney-Rivlin
PetscErrorCode ProcessPhysics_MR(MPI_Comm comm, Physics_MR phys, Units units);

#endif // mooney_rivlin_h
