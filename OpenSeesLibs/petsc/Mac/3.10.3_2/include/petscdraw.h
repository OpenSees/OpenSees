/*
  Interface to the PETSc graphics
*/
#if !defined(__PETSCDRAW_H)
#define __PETSCDRAW_H
#include <petscsys.h>
#include <petscdrawtypes.h>

PETSC_EXTERN PetscClassId PETSC_DRAW_CLASSID;

PETSC_EXTERN PetscFunctionList PetscDrawList;
PETSC_EXTERN PetscErrorCode PetscDrawInitializePackage(void);
PETSC_EXTERN PetscErrorCode PetscDrawFinalizePackage(void);
PETSC_EXTERN PetscErrorCode PetscDrawRegister(const char[],PetscErrorCode(*)(PetscDraw));

PETSC_EXTERN PetscErrorCode PetscDrawGetType(PetscDraw,PetscDrawType*);
PETSC_EXTERN PetscErrorCode PetscDrawSetType(PetscDraw,PetscDrawType);
PETSC_EXTERN PetscErrorCode PetscDrawCreate(MPI_Comm,const char[],const char[],int,int,int,int,PetscDraw*);
PETSC_EXTERN PetscErrorCode PetscDrawSetOptionsPrefix(PetscDraw,const char[]);
PETSC_EXTERN PetscErrorCode PetscDrawSetFromOptions(PetscDraw);
PETSC_EXTERN PetscErrorCode PetscDrawSetSave(PetscDraw,const char[]);
PETSC_EXTERN PetscErrorCode PetscDrawSetSaveMovie(PetscDraw,const char[]);
PETSC_EXTERN PetscErrorCode PetscDrawSetSaveFinalImage(PetscDraw,const char[]);
PETSC_EXTERN PetscErrorCode PetscDrawView(PetscDraw,PetscViewer);
PETSC_STATIC_INLINE PetscErrorCode PetscDrawViewFromOptions(PetscDraw A,PetscObject obj,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,obj,name);}

/*
   Number of basic colors in the draw routines, the others are used
   for a uniform colormap.
*/
#define PETSC_DRAW_BASIC_COLORS 33

#define PETSC_DRAW_ROTATE          -1         /* will rotate through the colors, start with 2 */
#define PETSC_DRAW_WHITE            0
#define PETSC_DRAW_BLACK            1
#define PETSC_DRAW_RED              2
#define PETSC_DRAW_GREEN            3
#define PETSC_DRAW_CYAN             4
#define PETSC_DRAW_BLUE             5
#define PETSC_DRAW_MAGENTA          6
#define PETSC_DRAW_AQUAMARINE       7
#define PETSC_DRAW_FORESTGREEN      8
#define PETSC_DRAW_ORANGE           9
#define PETSC_DRAW_VIOLET          10
#define PETSC_DRAW_BROWN           11
#define PETSC_DRAW_PINK            12
#define PETSC_DRAW_CORAL           13
#define PETSC_DRAW_GRAY            14
#define PETSC_DRAW_YELLOW          15

#define PETSC_DRAW_GOLD            16
#define PETSC_DRAW_LIGHTPINK       17
#define PETSC_DRAW_MEDIUMTURQUOISE 18
#define PETSC_DRAW_KHAKI           19
#define PETSC_DRAW_DIMGRAY         20
#define PETSC_DRAW_YELLOWGREEN     21
#define PETSC_DRAW_SKYBLUE         22
#define PETSC_DRAW_DARKGREEN       23
#define PETSC_DRAW_NAVYBLUE        24
#define PETSC_DRAW_SANDYBROWN      25
#define PETSC_DRAW_CADETBLUE       26
#define PETSC_DRAW_POWDERBLUE      27
#define PETSC_DRAW_DEEPPINK        28
#define PETSC_DRAW_THISTLE         29
#define PETSC_DRAW_LIMEGREEN       30
#define PETSC_DRAW_LAVENDERBLUSH   31
#define PETSC_DRAW_PLUM            32
#define PETSC_DRAW_MAXCOLOR       256

/*MC

   PetscDrawRealToColor - Maps a real value within an interval to a color.
   The color is an integer value in the range [PETSC_DRAW_BASIC_COLORS to 255]
   that can be passed to various drawing routines.

   Synopsis:
   #include <petscdraw.h>
   int PetscDrawRealToColor(PetscReal value,PetscReal min,PetscReal max)

   Not Collective

   Input Parameter:
+  value - value to map within the interval [min,max]
.  min - lower end of interval
-  max - upper end of interval

   Notes:
    Values outside the interval [min,max] are clipped.

   Level: intermediate

.seealso: PetscDrawPointPixel(), PetscDrawPoint(), PetscDrawLine(), PetscDrawTriangle(), PetscDrawRectangle()

M*/
PETSC_STATIC_INLINE int PetscDrawRealToColor(PetscReal value,PetscReal min,PetscReal max)
{
  value = PetscClipInterval(value,min,max);
  return PETSC_DRAW_BASIC_COLORS + (int)((255-PETSC_DRAW_BASIC_COLORS)*(value-min)/(max-min));
}

PETSC_EXTERN PetscErrorCode PetscDrawOpenX(MPI_Comm,const char[],const char[],int,int,int,int,PetscDraw*);
PETSC_EXTERN PetscErrorCode PetscDrawOpenGLUT(MPI_Comm,const char[],const char[],int,int,int,int,PetscDraw*);

#define PETSC_DRAW_FULL_SIZE    -3
#define PETSC_DRAW_HALF_SIZE    -4
#define PETSC_DRAW_THIRD_SIZE   -5
#define PETSC_DRAW_QUARTER_SIZE -6

PETSC_EXTERN PetscErrorCode PetscDrawOpenImage(MPI_Comm,const char[],int,int,PetscDraw*);
PETSC_EXTERN PetscErrorCode PetscDrawOpenNull(MPI_Comm,PetscDraw *);
PETSC_EXTERN PetscErrorCode PetscDrawDestroy(PetscDraw*);
PETSC_EXTERN PetscErrorCode PetscDrawIsNull(PetscDraw,PetscBool *);

PETSC_EXTERN PetscErrorCode PetscDrawGetPopup(PetscDraw,PetscDraw*);
PETSC_EXTERN PetscErrorCode PetscDrawScalePopup(PetscDraw,PetscReal,PetscReal);

PETSC_EXTERN PetscErrorCode PetscDrawCheckResizedWindow(PetscDraw);
PETSC_EXTERN PetscErrorCode PetscDrawResizeWindow(PetscDraw,int,int);
PETSC_EXTERN PetscErrorCode PetscDrawGetWindowSize(PetscDraw,int*,int*);
PETSC_EXTERN PetscErrorCode PetscDrawPixelToCoordinate(PetscDraw,int,int,PetscReal*,PetscReal*);
PETSC_EXTERN PetscErrorCode PetscDrawCoordinateToPixel(PetscDraw,PetscReal,PetscReal,int*,int*);

PETSC_EXTERN PetscErrorCode PetscDrawIndicatorFunction(PetscDraw,PetscReal,PetscReal,PetscReal,PetscReal,int,PetscErrorCode (*)(void*,PetscReal,PetscReal,PetscBool*),void*);

PETSC_EXTERN PetscErrorCode PetscDrawLine(PetscDraw,PetscReal,PetscReal,PetscReal,PetscReal,int);
PETSC_EXTERN PetscErrorCode PetscDrawArrow(PetscDraw,PetscReal,PetscReal,PetscReal,PetscReal,int);
PETSC_EXTERN PetscErrorCode PetscDrawLineSetWidth(PetscDraw,PetscReal);
PETSC_EXTERN PetscErrorCode PetscDrawLineGetWidth(PetscDraw,PetscReal*);

/*E
    PetscDrawMarkerType - How a "mark" is indicate in a figure

   Level: intermediate

$  PETSC_MARKER_CROSS - a small pixel based x symbol or the character x if that is not available
$  PETSC_MARKER_PLUS - a small pixel based + symbol or the character + if that is not available
$  PETSC_MARKER_CIRCLE - a small pixel based circle symbol or the character o if that is not available
$  PETSC_MARKER_POINT - the make obtained with PetscDrawPoint()

.seealso: PetscDrawMarker(), PetscDrawSetMarkerType()
E*/
typedef enum {PETSC_DRAW_MARKER_CROSS, PETSC_DRAW_MARKER_POINT,PETSC_DRAW_MARKER_PLUS,PETSC_DRAW_MARKER_CIRCLE} PetscDrawMarkerType;
PETSC_EXTERN const char *const PetscDrawMarkerTypes[];

PETSC_EXTERN PetscErrorCode PetscDrawMarker(PetscDraw,PetscReal,PetscReal,int);
PETSC_EXTERN PetscErrorCode PetscDrawSetMarkerType(PetscDraw,PetscDrawMarkerType);
PETSC_EXTERN PetscErrorCode PetscDrawGetMarkerType(PetscDraw,PetscDrawMarkerType*);

PETSC_EXTERN PetscErrorCode PetscDrawPoint(PetscDraw,PetscReal,PetscReal,int);
PETSC_EXTERN PetscErrorCode PetscDrawPointPixel(PetscDraw,int,int,int);
PETSC_EXTERN PetscErrorCode PetscDrawPointSetSize(PetscDraw,PetscReal);

PETSC_EXTERN PetscErrorCode PetscDrawRectangle(PetscDraw,PetscReal,PetscReal,PetscReal,PetscReal,int,int,int,int);
PETSC_EXTERN PetscErrorCode PetscDrawTriangle(PetscDraw,PetscReal,PetscReal,PetscReal,PetscReal,PetscReal,PetscReal,int,int,int);
PETSC_EXTERN PetscErrorCode PetscDrawEllipse(PetscDraw,PetscReal,PetscReal,PetscReal,PetscReal,int);
PETSC_EXTERN PetscErrorCode PetscDrawTensorContourPatch(PetscDraw,int,int,PetscReal*,PetscReal*,PetscReal,PetscReal,PetscReal*);
PETSC_EXTERN PetscErrorCode PetscDrawTensorContour(PetscDraw,int,int,const PetscReal[],const PetscReal[],PetscReal *);

PETSC_EXTERN PetscErrorCode PetscDrawString(PetscDraw,PetscReal,PetscReal,int,const char[]);
PETSC_EXTERN PetscErrorCode PetscDrawStringCentered(PetscDraw,PetscReal,PetscReal,int,const char[]);
PETSC_EXTERN PetscErrorCode PetscDrawStringBoxed(PetscDraw,PetscReal,PetscReal,int,int,const char[],PetscReal*,PetscReal*);
PETSC_EXTERN PetscErrorCode PetscDrawStringVertical(PetscDraw,PetscReal,PetscReal,int,const char[]);
PETSC_EXTERN PetscErrorCode PetscDrawStringSetSize(PetscDraw,PetscReal,PetscReal);
PETSC_EXTERN PetscErrorCode PetscDrawStringGetSize(PetscDraw,PetscReal*,PetscReal*);

PETSC_EXTERN PetscErrorCode PetscDrawSetViewPort(PetscDraw,PetscReal,PetscReal,PetscReal,PetscReal);
PETSC_EXTERN PetscErrorCode PetscDrawGetViewPort(PetscDraw,PetscReal*,PetscReal*,PetscReal*,PetscReal*);
PETSC_EXTERN PetscErrorCode PetscDrawSplitViewPort(PetscDraw);

PETSC_EXTERN PetscErrorCode PetscDrawSetCoordinates(PetscDraw,PetscReal,PetscReal,PetscReal,PetscReal);
PETSC_EXTERN PetscErrorCode PetscDrawGetCoordinates(PetscDraw,PetscReal*,PetscReal*,PetscReal*,PetscReal*);

PETSC_EXTERN PetscErrorCode PetscDrawSetTitle(PetscDraw,const char[]);
PETSC_EXTERN PetscErrorCode PetscDrawAppendTitle(PetscDraw,const char[]);
PETSC_EXTERN PetscErrorCode PetscDrawGetTitle(PetscDraw,const char *[]);

PETSC_EXTERN PetscErrorCode PetscDrawSetPause(PetscDraw,PetscReal);
PETSC_EXTERN PetscErrorCode PetscDrawGetPause(PetscDraw,PetscReal*);
PETSC_EXTERN PetscErrorCode PetscDrawPause(PetscDraw);
PETSC_EXTERN PetscErrorCode PetscDrawSetDoubleBuffer(PetscDraw);
PETSC_EXTERN PetscErrorCode PetscDrawClear(PetscDraw);
PETSC_EXTERN PetscErrorCode PetscDrawFlush(PetscDraw);
PETSC_EXTERN PetscErrorCode PetscDrawSave(PetscDraw);
PETSC_EXTERN PetscErrorCode PetscDrawSaveMovie(PetscDraw);
PETSC_EXTERN PetscErrorCode PetscDrawBOP(PetscDraw);
PETSC_EXTERN PetscErrorCode PetscDrawEOP(PetscDraw);

PETSC_EXTERN PetscErrorCode PetscDrawSetDisplay(PetscDraw,const char[]);
PETSC_EXTERN PetscErrorCode PetscDrawGetSingleton(PetscDraw,PetscDraw*);
PETSC_EXTERN PetscErrorCode PetscDrawRestoreSingleton(PetscDraw,PetscDraw*);

PETSC_EXTERN PetscErrorCode PetscDrawGetCurrentPoint(PetscDraw,PetscReal*,PetscReal*);
PETSC_EXTERN PetscErrorCode PetscDrawSetCurrentPoint(PetscDraw,PetscReal,PetscReal);
PETSC_EXTERN PetscErrorCode PetscDrawPushCurrentPoint(PetscDraw,PetscReal,PetscReal);
PETSC_EXTERN PetscErrorCode PetscDrawPopCurrentPoint(PetscDraw);
PETSC_EXTERN PetscErrorCode PetscDrawGetBoundingBox(PetscDraw,PetscReal*,PetscReal*,PetscReal*,PetscReal*);

/*E
    PetscDrawButton - Used to determine which button was pressed

   Level: intermediate

.seealso: PetscDrawGetMouseButton()
E*/
typedef enum { PETSC_BUTTON_NONE = 0,
               PETSC_BUTTON_LEFT,PETSC_BUTTON_CENTER,PETSC_BUTTON_RIGHT,
               PETSC_BUTTON_WHEEL_UP,PETSC_BUTTON_WHEEL_DOWN,
               PETSC_BUTTON_LEFT_SHIFT,PETSC_BUTTON_CENTER_SHIFT,PETSC_BUTTON_RIGHT_SHIFT } PetscDrawButton;

PETSC_EXTERN PetscErrorCode PetscDrawGetMouseButton(PetscDraw,PetscDrawButton *,PetscReal*,PetscReal *,PetscReal *,PetscReal *);

PETSC_EXTERN PetscErrorCode PetscDrawZoom(PetscDraw,PetscErrorCode (*)(PetscDraw,void *),void *);

/*S
     PetscDrawViewPorts - Subwindows in a PetscDraw object

   Level: intermediate

  Concepts: graphics

.seealso:  PetscDrawViewPortsCreate(), PetscDrawViewPortsSet()
S*/
typedef struct {
  PetscInt  nports;
  PetscReal *xl;
  PetscReal *xr;
  PetscReal *yl;
  PetscReal *yr;
  PetscDraw draw;
  PetscReal port_xl,port_yl,port_xr,port_yr;   /* original port of parent PetscDraw */
} PetscDrawViewPorts;
PETSC_EXTERN PetscErrorCode PetscDrawViewPortsCreate(PetscDraw,PetscInt,PetscDrawViewPorts**);
PETSC_EXTERN PetscErrorCode PetscDrawViewPortsCreateRect(PetscDraw,PetscInt,PetscInt,PetscDrawViewPorts**);
PETSC_EXTERN PetscErrorCode PetscDrawViewPortsDestroy(PetscDrawViewPorts*);
PETSC_EXTERN PetscErrorCode PetscDrawViewPortsSet(PetscDrawViewPorts*,PetscInt);

PETSC_EXTERN PetscClassId PETSC_DRAWAXIS_CLASSID;

PETSC_EXTERN PetscErrorCode PetscDrawAxisCreate(PetscDraw,PetscDrawAxis*);
PETSC_EXTERN PetscErrorCode PetscDrawAxisDestroy(PetscDrawAxis*);
PETSC_EXTERN PetscErrorCode PetscDrawAxisDraw(PetscDrawAxis);
PETSC_EXTERN PetscErrorCode PetscDrawAxisSetLimits(PetscDrawAxis,PetscReal,PetscReal,PetscReal,PetscReal);
PETSC_EXTERN PetscErrorCode PetscDrawAxisGetLimits(PetscDrawAxis,PetscReal*,PetscReal*,PetscReal*,PetscReal*);
PETSC_EXTERN PetscErrorCode PetscDrawAxisSetHoldLimits(PetscDrawAxis,PetscBool);
PETSC_EXTERN PetscErrorCode PetscDrawAxisSetColors(PetscDrawAxis,int,int,int);
PETSC_EXTERN PetscErrorCode PetscDrawAxisSetLabels(PetscDrawAxis,const char[],const char[],const char[]);

PETSC_EXTERN PetscClassId PETSC_DRAWLG_CLASSID;

PETSC_EXTERN PetscErrorCode PetscDrawLGCreate(PetscDraw,PetscInt,PetscDrawLG*);
PETSC_EXTERN PetscErrorCode PetscDrawLGDestroy(PetscDrawLG*);
PETSC_EXTERN PetscErrorCode PetscDrawLGAddPoint(PetscDrawLG,const PetscReal*,const PetscReal*);
PETSC_EXTERN PetscErrorCode PetscDrawLGAddCommonPoint(PetscDrawLG,const PetscReal,const PetscReal*);
PETSC_EXTERN PetscErrorCode PetscDrawLGAddPoints(PetscDrawLG,PetscInt,PetscReal**,PetscReal**);
PETSC_EXTERN PetscErrorCode PetscDrawLGDraw(PetscDrawLG);
PETSC_EXTERN PetscErrorCode PetscDrawLGSave(PetscDrawLG);
PETSC_EXTERN PetscErrorCode PetscDrawLGView(PetscDrawLG,PetscViewer);
PETSC_EXTERN PetscErrorCode PetscDrawLGReset(PetscDrawLG);
PETSC_EXTERN PetscErrorCode PetscDrawLGSetDimension(PetscDrawLG,PetscInt);
PETSC_EXTERN PetscErrorCode PetscDrawLGGetDimension(PetscDrawLG,PetscInt*);
PETSC_EXTERN PetscErrorCode PetscDrawLGSetLegend(PetscDrawLG,const char *const*);
PETSC_EXTERN PetscErrorCode PetscDrawLGGetAxis(PetscDrawLG,PetscDrawAxis*);
PETSC_EXTERN PetscErrorCode PetscDrawLGGetDraw(PetscDrawLG,PetscDraw*);
PETSC_EXTERN PetscErrorCode PetscDrawLGSetUseMarkers(PetscDrawLG,PetscBool);
PETSC_EXTERN PetscErrorCode PetscDrawLGSetLimits(PetscDrawLG,PetscReal,PetscReal,PetscReal,PetscReal);
PETSC_EXTERN PetscErrorCode PetscDrawLGSetColors(PetscDrawLG,const int[]);
PETSC_EXTERN PetscErrorCode PetscDrawLGSetOptionsPrefix(PetscDrawLG,const char[]);
PETSC_EXTERN PetscErrorCode PetscDrawLGSetFromOptions(PetscDrawLG);

PETSC_EXTERN PetscClassId PETSC_DRAWSP_CLASSID;

PETSC_EXTERN PetscErrorCode PetscDrawSPCreate(PetscDraw,int,PetscDrawSP*);
PETSC_EXTERN PetscErrorCode PetscDrawSPDestroy(PetscDrawSP*);
PETSC_EXTERN PetscErrorCode PetscDrawSPAddPoint(PetscDrawSP,PetscReal*,PetscReal*);
PETSC_EXTERN PetscErrorCode PetscDrawSPAddPoints(PetscDrawSP,int,PetscReal**,PetscReal**);
PETSC_EXTERN PetscErrorCode PetscDrawSPDraw(PetscDrawSP,PetscBool);
PETSC_EXTERN PetscErrorCode PetscDrawSPSave(PetscDrawSP);
PETSC_EXTERN PetscErrorCode PetscDrawSPReset(PetscDrawSP);
PETSC_EXTERN PetscErrorCode PetscDrawSPSetDimension(PetscDrawSP,int);
PETSC_EXTERN PetscErrorCode PetscDrawSPGetAxis(PetscDrawSP,PetscDrawAxis*);
PETSC_EXTERN PetscErrorCode PetscDrawSPGetDraw(PetscDrawSP,PetscDraw*);
PETSC_EXTERN PetscErrorCode PetscDrawSPSetLimits(PetscDrawSP,PetscReal,PetscReal,PetscReal,PetscReal);
PETSC_EXTERN PetscErrorCode PetscDrawLGSPDraw(PetscDrawLG,PetscDrawSP);

PETSC_EXTERN PetscClassId PETSC_DRAWHG_CLASSID;

PETSC_EXTERN PetscErrorCode PetscDrawHGCreate(PetscDraw,int,PetscDrawHG*);
PETSC_EXTERN PetscErrorCode PetscDrawHGDestroy(PetscDrawHG*);
PETSC_EXTERN PetscErrorCode PetscDrawHGAddValue(PetscDrawHG,PetscReal);
PETSC_EXTERN PetscErrorCode PetscDrawHGDraw(PetscDrawHG);
PETSC_EXTERN PetscErrorCode PetscDrawHGSave(PetscDrawHG);
PETSC_EXTERN PetscErrorCode PetscDrawHGView(PetscDrawHG,PetscViewer);
PETSC_EXTERN PetscErrorCode PetscDrawHGReset(PetscDrawHG);
PETSC_EXTERN PetscErrorCode PetscDrawHGGetAxis(PetscDrawHG,PetscDrawAxis*);
PETSC_EXTERN PetscErrorCode PetscDrawHGGetDraw(PetscDrawHG,PetscDraw*);
PETSC_EXTERN PetscErrorCode PetscDrawHGSetLimits(PetscDrawHG,PetscReal,PetscReal,int,int);
PETSC_EXTERN PetscErrorCode PetscDrawHGSetNumberBins(PetscDrawHG,int);
PETSC_EXTERN PetscErrorCode PetscDrawHGSetColor(PetscDrawHG,int);
PETSC_EXTERN PetscErrorCode PetscDrawHGCalcStats(PetscDrawHG, PetscBool );
PETSC_EXTERN PetscErrorCode PetscDrawHGIntegerBins(PetscDrawHG, PetscBool );

PETSC_EXTERN PetscClassId PETSC_DRAWBAR_CLASSID;

PETSC_EXTERN PetscErrorCode PetscDrawBarCreate(PetscDraw,PetscDrawBar*);
PETSC_EXTERN PetscErrorCode PetscDrawBarSetData(PetscDrawBar,PetscInt,const PetscReal[],const char *const*);
PETSC_EXTERN PetscErrorCode PetscDrawBarDestroy(PetscDrawBar*);
PETSC_EXTERN PetscErrorCode PetscDrawBarDraw(PetscDrawBar);
PETSC_EXTERN PetscErrorCode PetscDrawBarSave(PetscDrawBar);
PETSC_EXTERN PetscErrorCode PetscDrawBarSetColor(PetscDrawBar,int);
PETSC_EXTERN PetscErrorCode PetscDrawBarSetLimits(PetscDrawBar,PetscReal,PetscReal);
PETSC_EXTERN PetscErrorCode PetscDrawBarSort(PetscDrawBar,PetscBool,PetscReal);
PETSC_EXTERN PetscErrorCode PetscDrawBarSetFromOptions(PetscDrawBar);
PETSC_EXTERN PetscErrorCode PetscDrawBarGetAxis(PetscDrawBar,PetscDrawAxis*);
PETSC_EXTERN PetscErrorCode PetscDrawBarGetDraw(PetscDrawBar,PetscDraw*);

PETSC_EXTERN PetscErrorCode PetscDrawUtilitySetCmap(const char[],int,unsigned char[],unsigned char[],unsigned char[]);
PETSC_EXTERN PetscErrorCode PetscDrawUtilitySetGamma(PetscReal);

/*
    Handling of X11 I/O errors in parallel
*/
#if defined(PETSC_HAVE_X) && defined(PETSC_HAVE_SETJMP_H)

#include <setjmp.h>

PETSC_EXTERN jmp_buf PetscXIOErrorHandlerJumpBuf;
PETSC_EXTERN void    PetscXIOErrorHandlerJump(void *);
PETSC_EXTERN_TYPEDEF typedef void (*PetscXIOErrorHandler)(void *);
PETSC_EXTERN PetscXIOErrorHandler PetscSetXIOErrorHandler(PetscXIOErrorHandler);

#define PetscDrawCollectiveBegin(draw) 0; do { \
  PetscErrorCode                _Petsc_ierr; \
  jmp_buf                       _Petsc_jmpbuf; \
  volatile PetscXIOErrorHandler _Petsc_xioerrhdl = NULL; \
  PetscBool                     _Petsc_isdrawx, _Petsc_xioerr, _Petsc_xioerr_local = PETSC_FALSE; \
  _Petsc_ierr = PetscObjectTypeCompare((PetscObject)(draw),PETSC_DRAW_X,&_Petsc_isdrawx);CHKERRQ(_Petsc_ierr); \
  if (_Petsc_isdrawx) { \
    _Petsc_ierr = PetscMemcpy(&_Petsc_jmpbuf,&PetscXIOErrorHandlerJumpBuf,sizeof(_Petsc_jmpbuf));CHKERRQ(_Petsc_ierr); \
    _Petsc_xioerrhdl = PetscSetXIOErrorHandler(PetscXIOErrorHandlerJump); \
    if (setjmp(PetscXIOErrorHandlerJumpBuf)) {_Petsc_xioerr_local = PETSC_TRUE; do {_Petsc_ierr = PetscDrawCollectiveEnd(draw);CHKERRQ(_Petsc_ierr);} \
  } do {} while(0)

#define PetscDrawCollectiveEnd(draw) 0; \
  if (_Petsc_isdrawx) { \
    (void)PetscSetXIOErrorHandler(_Petsc_xioerrhdl); \
    _Petsc_ierr = PetscMemcpy(&PetscXIOErrorHandlerJumpBuf,&_Petsc_jmpbuf,sizeof(PetscXIOErrorHandlerJumpBuf));CHKERRQ(_Petsc_ierr); \
    _Petsc_ierr = MPI_Allreduce(&_Petsc_xioerr_local,&_Petsc_xioerr,1,MPIU_BOOL,MPI_LOR,PetscObjectComm((PetscObject)(draw)));CHKERRQ(_Petsc_ierr); \
    if (_Petsc_xioerr) { \
      _Petsc_ierr = PetscDrawSetType((draw),PETSC_DRAW_NULL);CHKERRQ(_Petsc_ierr); \
      PetscFunctionReturn(0); \
    } \
  } } while(0)

#else

#define PetscDrawCollectiveBegin(draw) 0; do {} while(0)
#define PetscDrawCollectiveEnd(draw)   0; do {} while(0)

#endif

#endif
