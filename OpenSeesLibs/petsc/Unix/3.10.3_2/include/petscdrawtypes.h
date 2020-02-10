#ifndef _PETSCDRAWTYPES_H
#define _PETSCDRAWTYPES_H

/*J
    PetscDrawType - String with the name of a PetscDraw

   Level: beginner

.seealso: PetscDrawSetType(), PetscDraw, PetscViewer, PetscDrawCreate()
J*/
typedef const char* PetscDrawType;
#define PETSC_DRAW_X          "x"
#define PETSC_DRAW_GLUT       "glut"
#define PETSC_DRAW_OPENGLES   "opengles"
#define PETSC_DRAW_NULL       "null"
#define PETSC_DRAW_WIN32      "win32"
#define PETSC_DRAW_TIKZ       "tikz"
#define PETSC_DRAW_IMAGE      "image"

/*S
     PetscDraw - Abstract PETSc object for graphics

   Level: beginner

  Concepts: graphics

.seealso:  PetscDrawCreate(), PetscDrawSetType(), PetscDrawType
S*/
typedef struct _p_PetscDraw* PetscDraw;

/*S
     PetscDrawAxis - Manages X-Y axis

   Level: advanced

  Concepts: graphics, axis

.seealso:  PetscDrawAxisCreate(), PetscDrawAxisSetLimits(), PetscDrawAxisSetColors(), PetscDrawAxisSetLabels()
S*/
typedef struct _p_PetscDrawAxis* PetscDrawAxis;

/*S
     PetscDrawLG - Manages drawing x-y plots

   Level: advanced

  Concepts: graphics, axis

.seealso:  PetscDrawAxisCreate(), PetscDrawLGCreate(), PetscDrawLGAddPoint()
S*/
typedef struct _p_PetscDrawLG*   PetscDrawLG;

/*S
     PetscDrawSP - Manages drawing scatter plots

   Level: advanced

  Concepts: graphics, scatter plots

.seealso:  PetscDrawSPCreate()
S*/
typedef struct _p_PetscDrawSP*   PetscDrawSP;

/*S
     PetscDrawHG - Manages drawing histograms

   Level: advanced

  Concepts: graphics, histograms

.seealso:  PetscDrawHGCreate()
S*/
typedef struct _p_PetscDrawHG*   PetscDrawHG;

/*S
     PetscDrawBar - Manages drawing bar graphs

   Level: advanced

  Concepts: graphics, histograms

.seealso:  PetscDrawBarCreate()
S*/
typedef struct _p_PetscDrawBar*   PetscDrawBar;

#endif
