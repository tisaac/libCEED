#if PETSC_VERSION_LT(3,14,0)
#  define DMPlexGetClosureIndices(a,b,c,d,e,f,g,h,i) DMPlexGetClosureIndices(a,b,c,d,f,g,i)
#  define DMPlexRestoreClosureIndices(a,b,c,d,e,f,g,h,i) DMPlexRestoreClosureIndices(a,b,c,d,f,g,i)
#endif

#if PETSC_VERSION_LT(3,14,0)
#  define DMAddBoundary(a,b,c,d,e,f,g,h,i,j,k,l) DMAddBoundary(a,b,c,d,e,f,g,h,j,k,l)
#endif

#if PETSC_VERSION_LT(3,14,0)
#  define DMPlexCreateSphereMesh(a,b,c,d,e) DMPlexCreateSphereMesh(a,b,c,e)
#endif

#if PETSC_VERSION_LT(3,14,0)
#  define DMPlexCreateSphereMesh(a,b,c,d,e) DMPlexCreateSphereMesh(a,b,c,e)
#endif