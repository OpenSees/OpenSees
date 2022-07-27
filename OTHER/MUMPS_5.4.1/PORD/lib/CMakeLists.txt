add_library(MUMPS_POORD OBJECT
   bucket.c
   ddbisect.c
   ddcreate.c
   gbipart.c
   gbisect.c
   gelim.c
   graph.c
   interface.c
   minpriority.c
   multisector.c
   nestdiss.c
   sort.c
   symbfac.c
   tree.c
)

target_link_libraries(MUMPS PUBLIC MUMPS_POORD)
   
