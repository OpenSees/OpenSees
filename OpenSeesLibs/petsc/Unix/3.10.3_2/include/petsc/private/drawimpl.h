/*
       Abstract data structure and functions for graphics.
*/

#if !defined(_DRAWIMPL_H)
#define _DRAWIMPL_H

#include <petscdraw.h>
#include <petsc/private/petscimpl.h>

PETSC_EXTERN PetscBool PetscDrawRegisterAllCalled;
PETSC_EXTERN PetscErrorCode PetscDrawRegisterAll(void);

struct _PetscDrawOps {
  PetscErrorCode (*setdoublebuffer)(PetscDraw);
  PetscErrorCode (*flush)(PetscDraw);
  PetscErrorCode (*line)(PetscDraw,PetscReal,PetscReal,PetscReal,PetscReal,int);
  PetscErrorCode (*linesetwidth)(PetscDraw,PetscReal);
  PetscErrorCode (*linegetwidth)(PetscDraw,PetscReal*);
  PetscErrorCode (*point)(PetscDraw,PetscReal,PetscReal,int);
  PetscErrorCode (*pointsetsize)(PetscDraw,PetscReal);
  PetscErrorCode (*string)(PetscDraw,PetscReal,PetscReal,int,const char[]);
  PetscErrorCode (*stringvertical)(PetscDraw,PetscReal,PetscReal,int,const char[]);
  PetscErrorCode (*stringsetsize)(PetscDraw,PetscReal,PetscReal);
  PetscErrorCode (*stringgetsize)(PetscDraw,PetscReal*,PetscReal*);
  PetscErrorCode (*setviewport)(PetscDraw,PetscReal,PetscReal,PetscReal,PetscReal);
  PetscErrorCode (*clear)(PetscDraw);
  PetscErrorCode (*rectangle)(PetscDraw,PetscReal,PetscReal,PetscReal,PetscReal,int,int,int,int);
  PetscErrorCode (*triangle)(PetscDraw,PetscReal,PetscReal,PetscReal,PetscReal,PetscReal,PetscReal,int,int,int);
  PetscErrorCode (*ellipse)(PetscDraw,PetscReal,PetscReal,PetscReal,PetscReal,int);
  PetscErrorCode (*getmousebutton)(PetscDraw,PetscDrawButton*,PetscReal *,PetscReal *,PetscReal*,PetscReal*);
  PetscErrorCode (*pause)(PetscDraw);
  PetscErrorCode (*beginpage)(PetscDraw);
  PetscErrorCode (*endpage)(PetscDraw);
  PetscErrorCode (*getpopup)(PetscDraw,PetscDraw*);
  PetscErrorCode (*settitle)(PetscDraw,const char[]);
  PetscErrorCode (*checkresizedwindow)(PetscDraw);
  PetscErrorCode (*resizewindow)(PetscDraw,int,int);
  PetscErrorCode (*destroy)(PetscDraw);
  PetscErrorCode (*view)(PetscDraw,PetscViewer);
  PetscErrorCode (*getsingleton)(PetscDraw,PetscDraw*);
  PetscErrorCode (*restoresingleton)(PetscDraw,PetscDraw*);
  PetscErrorCode (*save)(PetscDraw);
  PetscErrorCode (*getimage)(PetscDraw,unsigned char[][3],unsigned int*,unsigned int*,unsigned char*[]);
  PetscErrorCode (*setcoordinates)(PetscDraw,PetscReal,PetscReal,PetscReal,PetscReal);
  PetscErrorCode (*arrow)(PetscDraw,PetscReal,PetscReal,PetscReal,PetscReal,int);
  PetscErrorCode (*coordinatetopixel)(PetscDraw,PetscReal,PetscReal,int*,int*);
  PetscErrorCode (*pixeltocoordinate)(PetscDraw,int,int,PetscReal*,PetscReal*);
  PetscErrorCode (*pointpixel)(PetscDraw,int,int,int);
  PetscErrorCode (*boxedstring)(PetscDraw,PetscReal,PetscReal,int,int,const char[],PetscReal*,PetscReal*);
};

struct _p_PetscDraw {
  PETSCHEADER(struct _PetscDrawOps);
  PetscReal           pause;       /* sleep time after a synchronized flush */
  PetscReal           port_xl,port_yl,port_xr,port_yr;
  PetscReal           coor_xl,coor_yl,coor_xr,coor_yr;
  PetscReal           currentpoint_x[20],currentpoint_y[20];
  PetscReal           boundbox_xl,boundbox_yl,boundbox_xr,boundbox_yr; /* need to have this for each current point? */
  PetscInt            currentpoint;
  PetscDrawMarkerType markertype;
  char                *title;
  char                *display;
  PetscDraw           popup;
  int                 x,y,h,w;
  char                *savefilename;
  char                *saveimageext;
  char                *savemovieext;
  PetscInt            savefilecount;
  PetscBool           savesinglefile;
  PetscInt            savemoviefps;
  char                *savefinalfilename;
  PetscBool           saveonclear; /* save a new image for every PetscDrawClear() called */
  PetscBool           saveonflush; /* save a new image for every PetscDrawFlush() called */
  void                *data;
};

#endif
