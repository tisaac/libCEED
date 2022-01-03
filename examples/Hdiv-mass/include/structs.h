#ifndef structs_h
#define structs_h

#include <ceed.h>
#include <petsc.h>

// Application context from user command line options
typedef struct AppCtx_ *AppCtx;
struct AppCtx_ {
  // libCEED arguments
  PetscInt          degree;
  PetscInt          q_extra;
  // Problem type arguments
  PetscFunctionList problems;
  char              problem_name[PETSC_MAX_PATH_LEN];
};

// libCEED data struct
typedef struct CeedData_ *CeedData;
struct CeedData_ {
  CeedBasis            basis_x, basis_u;
  CeedElemRestriction  elem_restr_x, elem_restr_u, elem_restr_geo_data_i,
                       elem_restr_u_i;
  CeedQFunction        qf_residual, qf_error;
  CeedOperator         op_residual, op_error;
  CeedVector           x_ceed, y_ceed;
  CeedQFunctionContext pq2d_context;
};

// 1) poisson-quad2d
#ifndef PHYSICS_POISSONQUAD2D_STRUCT
#define PHYSICS_POISSONQUAD2D_STRUCT
typedef struct PQ2DContext_ *PQ2DContext;
struct PQ2DContext_ {
  CeedScalar kappa;
};
#endif

// 2) poisson-hex3d

// 3) poisson-prism3d

// 4) richard

// Struct that contains all enums and structs used for the physics of all problems
typedef struct Physics_ *Physics;
struct Physics_ {
  PQ2DContext            pq2d_ctx;
};

// PETSc user data
typedef struct User_ *User;
struct User_ {
  MPI_Comm     comm;
  Vec          X_loc, Y_loc;
  CeedVector   x_ceed, y_ceed;
  CeedOperator op_apply, op_error;
  DM           dm;
  Ceed         ceed;
  AppCtx       app_ctx;
  Physics      phys;
};

// Problem specific data
typedef struct {
  CeedQFunctionUser setup_rhs, residual, setup_error, setup_true;
  const char        *setup_rhs_loc, *residual_loc, *setup_error_loc,
        *setup_true_loc;
  CeedQuadMode      quadrature_mode;
  CeedInt           geo_data_size, elem_node;
  PetscErrorCode    (*setup_ctx)(Ceed, CeedData, Physics);

} ProblemData;

#endif // structs_h