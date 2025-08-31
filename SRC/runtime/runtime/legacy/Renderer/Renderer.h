/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Description: This file contains the class definition for Renderer.
// Renderer is an abstract base class. An Renderer object is used
// to process 3d graphics.
//
// Written: fmk
// Created: 10/98
// Revision: A
//
#ifndef Renderer_h
#define Renderer_h

class Vector;
class Matrix;
class ColorMap;

class Renderer
{
  public:
    Renderer(ColorMap &theMap);
    Renderer(const char *title, ColorMap &theMap);
    ~Renderer();

    // method to set the color map
    void setColorMap(ColorMap &theMap);

    // method to clear the current image
    int clearImage()  {return -1;}
    int saveImage(const char *imageName);
    int saveImage(const char *title, const char *imageName);

    // methods to be invoked when image processing is to start or is finished
    int startImage()  {return -1;}
    int doneImage()  {return -1;}

    // methods invoked by the objects to display themselves
    int drawPoint(const Vector &, float , int tag = 0, int mode = 0, int width = 1)  {return -1;}
    int drawPoint(const Vector &, const Vector &, int tag = 0, int mode = 0, int width = 1)  {return -1;}

    int drawLine(const Vector &, const Vector &,
                 float V1, float V2, int tag = 0, int mode = 0) {return -1;}
    int drawLine(const Vector &, const Vector &,
                 const Vector &, const Vector &, int tag = 0, int mode = 0) {return -1;}

    int drawCube(const Matrix &, const Vector &, int tag = 0, int mode = 0);

    int drawPolygon(const Matrix &, const Vector &, int tag = 0, int mode = 0)  {return -1;}
    int drawPolygon(const Matrix &, const Matrix &, int tag = 0, int mode = 0)  {return -1;}

    int drawVector(const Vector &position, const Vector &value, double fcator = 1.0, int tag = 0);


    int drawText(const Vector &, char *, int ,
                 char horizontalJustify = 'l',
                 char verticalJustify = 'b') {return -1;}

    //
    // the following are for setting up the viewing system
    //

    // the following are in world coordinates & define view coord system
    int setVRP(float , float , float )  {return -1;} // point on view plane
    int setVPN(float , float , float )  {return -1;} // view plane normal
    int setVUP(float , float , float )  {return -1;} // view-up vector
        
    // the following are in view coordinates        
    int setViewWindow(float, float, float, float)  {return -1;} // view bounds
                               // umin, umax, vmin, vmax

    // location of near and far clipping planes from view plane
    int setPlaneDist(float, float)  {return -1;} 

    int setProjectionMode(const char *)  {return -1;} //parallel or perspective
    int setFillMode(const char *)  {return -1;}    // wire or fill
    int setLineWidth(int)  {return -1;}    // line width

    // eye location if perspective, dirn to +ViewPlane if parallel
    int setPRP(float , float , float )  {return -1;} 

    // the following are in normalized coordinates
    int setPortWindow(float, float, float, float)  {return -1;} // view port
                    // left, right, bottom, top [-1,1,-1,1]
  protected:
    ColorMap *theMap;

  private:
    static int numRenderers;
    static char **theTitles;
    static Renderer **theRenderers;
};

#endif

