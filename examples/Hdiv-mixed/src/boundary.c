#include "../include/boundary.h"

// ---------------------------------------------------------------------------
// Create boundary label
// ---------------------------------------------------------------------------
PetscErrorCode CreateBCLabel(DM dm, const char name[]) {
  PetscErrorCode ierr;
  DMLabel        label;

  PetscFunctionBeginUser;

  ierr = DMCreateLabel(dm, name); CHKERRQ(ierr);
  ierr = DMGetLabel(dm, name, &label); CHKERRQ(ierr);
  ierr = DMPlexMarkBoundaryFaces(dm, PETSC_DETERMINE, label); CHKERRQ(ierr);
  ierr = DMPlexLabelComplete(dm, label); CHKERRQ(ierr);

  PetscFunctionReturn(0);
};

// ---------------------------------------------------------------------------
// Add Dirichlet boundaries to DM
// ---------------------------------------------------------------------------
PetscErrorCode DMAddBoundariesDirichlet(DM dm) {
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  // BCs given by manufactured solution
  PetscBool  has_label;
  const char *name = "MMS Face Sets";
  PetscInt   face_ids[1] = {1};
  ierr = DMHasLabel(dm, name, &has_label); CHKERRQ(ierr);
  if (!has_label) {
    ierr = CreateBCLabel(dm, name); CHKERRQ(ierr);
  }
  DMLabel label;
  ierr = DMGetLabel(dm, name, &label); CHKERRQ(ierr);
  ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "mms", label, 1, face_ids, 0, 0, NULL,
                       (void(*)(void))BoundaryDirichletMMS, NULL, NULL, NULL); CHKERRQ(ierr);
   

  PetscFunctionReturn(0);
}

#ifndef M_PI
#define M_PI    3.14159265358979323846
#endif
// ---------------------------------------------------------------------------
// Boundary function for manufactured solution
// ---------------------------------------------------------------------------
PetscErrorCode BoundaryDirichletMMS(PetscInt dim, PetscReal t, const PetscReal coords[], 
                                    PetscInt num_comp_u, PetscScalar *u, void *ctx) {
  PetscScalar x = coords[0];
  PetscScalar y = coords[1];
  PetscScalar z = coords[1];

  PetscFunctionBeginUser;

  if (dim == 2) {
    u[0] = -M_PI*cos(M_PI*x) *sin(M_PI*y);
    u[1] = -M_PI*sin(M_PI*x) *cos(M_PI*y);
  } else {
    u[0] = -M_PI*cos(M_PI*x) *sin(M_PI*y) *sin(M_PI*z);
    u[1] = -M_PI*sin(M_PI*x) *cos(M_PI*y) *sin(M_PI*z);
    u[2] = -M_PI*sin(M_PI*x) *sin(M_PI*y) *cos(M_PI*z);
  }
  

  PetscFunctionReturn(0);
}
