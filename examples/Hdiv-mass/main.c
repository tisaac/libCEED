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

//                        libCEED + PETSc Example: Mixed-Poisson in H(div) space
//
// This example demonstrates a simple usage of libCEED with PETSc to solve
//   elasticity problems.
//
// The code uses higher level communication protocols in DMPlex.
//
// Build with: make
// Run with:
//          ./main
//          ./main pc_type svd
const char help[] = "Solve H(div)-mixed problem using PETSc and libCEED\n";

#include "main.h"

int main(int argc, char **argv) {
  // ---------------------------------------------------------------------------
  // Initialize PETSc
  // ---------------------------------------------------------------------------
  PetscInt ierr;
  ierr = PetscInitialize(&argc, &argv, NULL, help);
  if (ierr) return ierr;

  // ---------------------------------------------------------------------------
  // Create structs
  // ---------------------------------------------------------------------------
  AppCtx app_ctx;
  ierr = PetscCalloc1(1, &app_ctx); CHKERRQ(ierr);

  ProblemData *problem_data = NULL;
  ierr = PetscCalloc1(1, &problem_data); CHKERRQ(ierr);

  User user;
  ierr = PetscCalloc1(1, &user); CHKERRQ(ierr);

  CeedData ceed_data;
  ierr = PetscCalloc1(1, &ceed_data); CHKERRQ(ierr);

  Physics phys_ctx;
  ierr = PetscCalloc1(1, &phys_ctx); CHKERRQ(ierr);

  user->app_ctx = app_ctx;
  user->phys    = phys_ctx;

  // ---------------------------------------------------------------------------
  // Process command line options
  // ---------------------------------------------------------------------------
  // -- Register problems to be available on the command line
  ierr = RegisterProblems_Hdiv(app_ctx); CHKERRQ(ierr);

  // -- Process general command line options
  MPI_Comm comm = PETSC_COMM_WORLD;
  ierr = ProcessCommandLineOptions(comm, app_ctx); CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // Choose the problem from the list of registered problems
  // ---------------------------------------------------------------------------
  {
    PetscErrorCode (*p)(ProblemData *, void *);
    ierr = PetscFunctionListFind(app_ctx->problems, app_ctx->problem_name, &p);
    CHKERRQ(ierr);
    if (!p) SETERRQ1(PETSC_COMM_SELF, 1, "Problem '%s' not found",
                       app_ctx->problem_name);
    ierr = (*p)(problem_data, &user); CHKERRQ(ierr);
  }

  // ---------------------------------------------------------------------------
  // Initialize libCEED
  // ---------------------------------------------------------------------------
  // -- Initialize backend
  Ceed ceed;
  CeedInit("/cpu/self/ref/serial", &ceed);
  CeedMemType mem_type_backend;
  CeedGetPreferredMemType(ceed, &mem_type_backend);
  // ---------------------------------------------------------------------------
  // Set-up DM
  // ---------------------------------------------------------------------------
  // PETSc objects
  DM             dm;
  VecType        vec_type;
  ierr = CreateDistributedDM(comm, &dm); CHKERRQ(ierr);
  ierr = DMGetVecType(dm, &vec_type); CHKERRQ(ierr);
  if (!vec_type) { // Not yet set by user -dm_vec_type
    switch (mem_type_backend) {
    case CEED_MEM_HOST: vec_type = VECSTANDARD; break;
    case CEED_MEM_DEVICE: {
      const char *resolved;
      CeedGetResource(ceed, &resolved);
      if (strstr(resolved, "/gpu/cuda")) vec_type = VECCUDA;
      else if (strstr(resolved, "/gpu/hip/occa"))
        vec_type = VECSTANDARD; // https://github.com/CEED/libCEED/issues/678
      else if (strstr(resolved, "/gpu/hip")) vec_type = VECHIP;
      else vec_type = VECSTANDARD;
    }
    }
    ierr = DMSetVecType(dm, vec_type); CHKERRQ(ierr);
  }
  // ---------------------------------------------------------------------------
  // Create global, local solution, local rhs vector
  // ---------------------------------------------------------------------------
  Vec            U_g, U_loc;
  PetscInt       U_l_size, U_g_size, U_loc_size;
  // Create global and local solution vectors
  ierr = DMCreateGlobalVector(dm, &U_g); CHKERRQ(ierr);
  ierr = VecGetSize(U_g, &U_g_size); CHKERRQ(ierr);
  // Local size for matShell
  ierr = VecGetLocalSize(U_g, &U_l_size); CHKERRQ(ierr);
  // Create local unknown vector U_loc
  ierr = DMCreateLocalVector(dm, &U_loc); CHKERRQ(ierr);
  // Local size for libCEED
  ierr = VecGetSize(U_loc, &U_loc_size); CHKERRQ(ierr);

  // Get RHS vector
  Vec rhs_loc;
  PetscScalar *r;
  CeedVector rhs_ceed, target;
  PetscMemType mem_type;
  ierr = VecDuplicate(U_loc, &rhs_loc); CHKERRQ(ierr);
  ierr = VecZeroEntries(rhs_loc); CHKERRQ(ierr);
  ierr = VecGetArrayAndMemType(rhs_loc, &r, &mem_type); CHKERRQ(ierr);
  CeedVectorCreate(ceed, U_l_size, &rhs_ceed);
  CeedVectorSetArray(rhs_ceed, MemTypeP2C(mem_type), CEED_USE_POINTER, r);
  // Get projected true solution
  Vec true_loc;
  PetscScalar *t;
  CeedVector true_ceed;
  PetscMemType t_mem_type;
  ierr = VecDuplicate(U_loc, &true_loc); CHKERRQ(ierr);
  ierr = VecZeroEntries(true_loc); CHKERRQ(ierr);
  ierr = VecGetArrayAndMemType(true_loc, &t, &t_mem_type); CHKERRQ(ierr);
  CeedVectorCreate(ceed, U_l_size, &true_ceed);
  CeedVectorSetArray(true_ceed, MemTypeP2C(t_mem_type), CEED_USE_POINTER, t);

  // ---------------------------------------------------------------------------
  // Setup libCEED
  // ---------------------------------------------------------------------------
  // -- Set up libCEED objects
  ierr = SetupLibceed(dm, ceed, app_ctx, problem_data, U_g_size,
                      U_loc_size, ceed_data, rhs_ceed, &target, true_ceed); CHKERRQ(ierr);

  //CeedVectorView(true_ceed, "%12.8f", stdout);
  // ---------------------------------------------------------------------------
  // Gather RHS
  // ---------------------------------------------------------------------------
  Vec rhs;
  CeedVectorTakeArray(rhs_ceed, MemTypeP2C(mem_type), NULL);
  ierr = VecRestoreArrayAndMemType(rhs_loc, &r); CHKERRQ(ierr);
  ierr = VecDuplicate(U_g, &rhs); CHKERRQ(ierr);
  ierr = VecZeroEntries(rhs); CHKERRQ(ierr);
  ierr = DMLocalToGlobal(dm, rhs_loc, ADD_VALUES, rhs); CHKERRQ(ierr);
  // ---------------------------------------------------------------------------
  // Setup Mat, KSP
  // ---------------------------------------------------------------------------
  user->comm = comm;
  user->dm = dm;
  user->X_loc = U_loc;
  ierr = VecDuplicate(U_loc, &user->Y_loc); CHKERRQ(ierr);
  user->x_ceed = ceed_data->x_ceed;
  user->y_ceed = ceed_data->y_ceed;
  user->op_apply = ceed_data->op_residual;
  user->op_error = ceed_data->op_error;
  user->ceed = ceed;
  // Operator
  Mat mat;
  ierr = MatCreateShell(comm, U_l_size, U_l_size, U_g_size, U_g_size,
                        user, &mat); CHKERRQ(ierr);
  ierr = MatShellSetOperation(mat, MATOP_MULT,
                              (void(*)(void))MatMult_Ceed); CHKERRQ(ierr);
  ierr = MatShellSetVecType(mat, vec_type); CHKERRQ(ierr);

  KSP ksp;
  ierr = KSPCreate(comm, &ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, mat, mat); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
  ierr = KSPSetUp(ksp); CHKERRQ(ierr);
  ierr = KSPSolve(ksp, rhs, U_g); CHKERRQ(ierr);
  //VecView(U_g, PETSC_VIEWER_STDOUT_WORLD);
  // ---------------------------------------------------------------------------
  // Compute pointwise L2 maximum error
  // ---------------------------------------------------------------------------
  CeedScalar l2_error;
  ierr = ComputeError(user, U_g, target, &l2_error); CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // Compute L2 error of projected solution into H(div) space
  // ---------------------------------------------------------------------------
  const CeedScalar *true_array;
  Vec error_vec, true_vec;

  // -- Work vectors
  ierr = VecDuplicate(U_g, &error_vec); CHKERRQ(ierr);
  ierr = VecSet(error_vec, 0.0); CHKERRQ(ierr);
  ierr = VecDuplicate(U_g, &true_vec); CHKERRQ(ierr);
  ierr = VecSet(true_vec, 0.0); CHKERRQ(ierr);

  // -- Assemble global true solution vector
  CeedVectorGetArrayRead(true_ceed, CEED_MEM_HOST, &true_array);
  ierr = VecPlaceArray(user->Y_loc, (PetscScalar *)true_array);
  CHKERRQ(ierr);
  ierr = DMLocalToGlobal(user->dm, user->Y_loc, INSERT_VALUES, true_vec);
  CHKERRQ(ierr);
  ierr = VecResetArray(user->Y_loc); CHKERRQ(ierr);
  CeedVectorRestoreArrayRead(true_ceed, &true_array);

  // -- Compute H(div) projected error
  CeedScalar proj_error;
  ierr = VecWAXPY(error_vec, -1.0, U_g, true_vec); CHKERRQ(ierr);
  ierr = VecNorm(error_vec, NORM_2, &proj_error); CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // Output results
  // ---------------------------------------------------------------------------
  KSPType ksp_type;
  KSPConvergedReason reason;
  PetscReal rnorm;
  PetscInt its;
  ierr = KSPGetType(ksp, &ksp_type); CHKERRQ(ierr);
  ierr = KSPGetConvergedReason(ksp, &reason); CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(ksp, &its); CHKERRQ(ierr);
  ierr = KSPGetResidualNorm(ksp, &rnorm); CHKERRQ(ierr);
  ierr = PetscPrintf(comm,
                     "  KSP:\n"
                     "    KSP Type                            : %s\n"
                     "    KSP Convergence                     : %s\n"
                     "    Total KSP Iterations                : %D\n"
                     "    Final rnorm                         : %e\n"
                     "    L2 Error                            : %e\n",
                     ksp_type, KSPConvergedReasons[reason], its,
                     (double)rnorm, (double)l2_error); CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // Free objects
  // ---------------------------------------------------------------------------

  // Free PETSc objects
  ierr = DMDestroy(&dm); CHKERRQ(ierr);
  ierr = VecDestroy(&U_g); CHKERRQ(ierr);
  ierr = VecDestroy(&U_loc); CHKERRQ(ierr);
  ierr = VecDestroy(&rhs); CHKERRQ(ierr);
  ierr = VecDestroy(&rhs_loc); CHKERRQ(ierr);
  ierr = VecDestroy(&user->Y_loc); CHKERRQ(ierr);
  ierr = MatDestroy(&mat); CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp); CHKERRQ(ierr);

  // -- Function list
  ierr = PetscFunctionListDestroy(&app_ctx->problems); CHKERRQ(ierr);

  // -- Structs
  ierr = PetscFree(app_ctx); CHKERRQ(ierr);
  ierr = PetscFree(problem_data); CHKERRQ(ierr);
  ierr = PetscFree(user); CHKERRQ(ierr);
  ierr = PetscFree(phys_ctx->pq2d_ctx); CHKERRQ(ierr);
  ierr = PetscFree(phys_ctx); CHKERRQ(ierr);

  // Free libCEED objects
  CeedVectorDestroy(&true_ceed);
  CeedVectorDestroy(&rhs_ceed);
  CeedVectorDestroy(&target);
  ierr = CeedDataDestroy(ceed_data); CHKERRQ(ierr);
  CeedDestroy(&ceed);

  return PetscFinalize();
}
