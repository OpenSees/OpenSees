/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.3 $
// $Date: 2000-12-19 04:02:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/OpenGlRenderer.cpp,v $
                                                                        
                                                                        
// File: ~/renderer/OpenGLRenderer.C
//
// Written: fmk 
// Created: 10/98
// Revision: A
//
// Description: This file contains the class definition for OpenGLRenderer.
// OpenGLRenderer is an class which diplays using X11 or openGL.
//
// What: "@(#) OpenGLRenderer.h, revA"

#include <OpenGLRenderer.h>
#include <ColorMap.h>
#include <stdio.h>
#include <stdlib.h>
#include <Matrix.h>

#define PARALLEL_MODE 0
#define PERSPECTIVE_MODE 1

//void itoa(int x, char *str);

LONG WINAPI WndProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{ 
    LONG        lRet = 1;
    PAINTSTRUCT ps;

    switch(uMsg) {
    case WM_CREATE:
        break; 

    case WM_DESTROY:
        break; 

    case WM_PAINT: 
        BeginPaint(hWnd, &ps); 
        EndPaint(hWnd, &ps); 
        break; 

    default: 
        lRet = DefWindowProc (hWnd, uMsg, wParam, lParam); 
        break; 
    }

    return lRet;
}


#ifdef _WIN32
int oglSetPixelFormat(HDC hDC, BYTE type, DWORD flags);
HWND oglCreateWindow(char* title, int x, int y, int width, int height,
		     HGLRC *hRC, HDC *hDC);
int oglDestroyWindow(char* title, HWND hWnd, HGLRC hRC, HDC hDC);
int oglCreateBitmap(int width, int height, HGLRC *hRC, HDC *hDC, 
						HBITMAP *theBitmap, BITMAPINFO *info, GLubyte **bits);


int oglDestroyBitmap(HBITMAP *bitmap, HGLRC hRC, HDC hDC);
#else

#endif


#include <db.H>

#include <Vector.h>
#include <WindowDevice.h>

OpenGLRenderer::OpenGLRenderer(char *_title, int _xLoc, int _yLoc, int _width, int _height,
			 ColorMap &_theMap)
  :Renderer(_theMap), winOpen(1), aFile(0), theFile(0), count(-1), bits(0)
{
	    // set the WindowDevices title, height, wdth, xLoc and yLoc
    strcpy(title, _title);

    height = _height;
    width = _width;  
    xLoc = _xLoc;
    yLoc = _yLoc;

#ifdef _WIN32   
    if (winOpen == 0)
      oglDestroyWindow(title,theWND, theHRC, theHDC);      

    theWND = oglCreateWindow(title, xLoc, yLoc, width, height, &theHRC, &theHDC);
    if (theWND == NULL)
      exit(1);
#else

#endif
    winOpen = 0;

    wglMakeCurrent(theHDC, theHRC);
    glClearColor(1.0f,1.0f,1.0f,1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

#ifdef _WIN32
	SwapBuffers(theHDC);
#endif
    glFlush();
}




OpenGLRenderer::OpenGLRenderer(char *_title, int _xLoc, int _yLoc, int _width, 
			       int _height,	
			       ColorMap &_theMap, char *fileName, char *bmpFileName)
  :Renderer(_theMap), winOpen(1), theFile(0), count(-1), bits(0)
{
	    // set the WindowDevices title, height, wdth, xLoc and yLoc
    strcpy(title, _title);

    height = _height;
    width = _width;  
    xLoc = _xLoc;
    yLoc = _yLoc;

#ifdef _WIN32   

#else

#endif
 



    // open the file for  making the movie
    if (fileName != 0) {
	if (strlen(fileName) > MAX_FILENAMELENGTH) 
	    cerr <<"warning - OpenGLRenderer::OpenGLRenderer() - fileName " << fileName << "too long\n";
	else {
	    strcpy(theFileName, fileName);
			/*
		FILE *fp;
		if ((fp = fopen(fileName,"w")) == NULL) {
		*/
	theFile.open(fileName, ios::out);
	if (!theFile) {
	  g3ErrorHandler->warning("WARNING - OpenGLRenderer::OpenGLRenderer() - could not open file %s\n",fileName);	  
	  aFile = 0;
	} else {
      aFile = 1;
	  theFile << title << endl;
	  theFile << xLoc << " " << yLoc << " " << width << " " << height << endl;
	}
     if (winOpen == 0)
      oglDestroyWindow(title,theWND, theHRC, theHDC);      

    theWND = oglCreateWindow(title, xLoc, yLoc, width, height, &theHRC, &theHDC);
    if (theWND == NULL)
      exit(1);
	winOpen = 0;

	    
	}
    }
    
    if (bmpFileName != 0) {
		count = 100000;
	if (strlen(bmpFileName) > MAX_FILENAMELENGTH) 
	    g3ErrorHandler->warning("warning - X11Renderer::X11Renderer() - bmpFileName %s too long, max %d\n",
				    bmpFileName, MAX_FILENAMELENGTH);
	else {
	    strcpy(theBmpFileName, bmpFileName); 
		oglCreateBitmap(width, height, &theHRC, &theHDC, &theBitmap, &info, &bits);
 
	}
    } 
    

    wglMakeCurrent(theHDC, theHRC);
    glClearColor(1.0f,1.0f,1.0f,1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

#ifdef _WIN32
//	SwapBuffers(theHDC);
#endif
	glFlush();
}

OpenGLRenderer::~OpenGLRenderer()
{
#ifdef _WIN32
   if (winOpen == 0) { // we must close the window
		if (count == -1)
       oglDestroyWindow(title,theWND, theHRC, theHDC);
		else 
		oglDestroyBitmap(&theBitmap, theHRC, theHDC);
   }
#else

#endif
   if (aFile == 1)
       theFile.close();
}


int 
OpenGLRenderer::clearImage(void)
{
#ifdef _WIN32
    wglMakeCurrent(theHDC, theHRC);
#else

#endif

    glClearColor(1.0f,1.0f,1.0f,1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
 
    glFlush(); 
    return 0;
}

int 
OpenGLRenderer::startImage(void)
{
#ifdef _WIN32
    wglMakeCurrent(theHDC, theHRC);

#else
  
#endif

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);

    //
    // set up the viewing transformation
    //
    VECTOR u, v, n;

    for (int i=0; i<3; i++) {
	n[i] = vpn[i];
	v[i] = vuv[i];
    }
   
    if (n.Normalize() == 0) {
	cerr << "View::update() - VPN cannot have zero length\n";
	return -1;
    }

    u = v % n;
    if (u.Normalize() == 0) {
	cerr << "View::update() - VUV X VPN cannot have zero length\n";
	return -1;
    }

    v = n % u;
    v.Normalize();

    ViewMat.Set(       u[0],   v[0],   n[0], 0,
         	       u[1],   v[1],   n[1], 0,
		       u[2],   v[2],   n[2], 0,
	      -vrp.Dot(u), -vrp.Dot(v), -vrp.Dot(n), 1.0);

    for (int j=0; j<4; j++)
	for (int k=0; k<4; k++)
	    viewData[j*4+k] = ViewMat.m[j][k];

    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixf(viewData);

   //
   // set up the projection transformation
   //

    float midU, midV, dopU, dopV, dopN, shU, shV, F, B;
    float diffU, diffV;
    midU = (vpwindow[0] + vpwindow[1])/2;
    midV = (vpwindow[2] + vpwindow[3])/2;
    dopN = -cop[2];
    dopU = midU - cop[0];
    dopV = midV - cop[1];

    diffU = vpwindow[1]-vpwindow[0];
    diffV = vpwindow[3]-vpwindow[2];      
  
    shU = dopU/dopN;
    shV = dopV/dopN;
    F = planedist[0];
    B = planedist[2];
    if (F < 0) F = -F; // warning MESSAGES !!!
    if (B < 0) B = -B;
  
    if (projection_mode == PARALLEL_MODE) {
      
	if (F == -B)
	    B = F;

	ProjMat.Set(2.0/diffU,          0.0,          0.0,    0.0,
       		          0.0,    2.0/diffV,         0.0,    0.0,
		  2*shU/diffV,  2*shV/diffV,      1/(F+B),    0.0,
	        -2*midU/diffU, -2*midV/diffV, -F/(F+B),    1.0);
     
    } else { // perspective projection
	float a,b,c,e,f,g,h;
	float diffU, diffV;// prpN;
	
	diffU = vpwindow[1]-vpwindow[0];
	diffV = vpwindow[3]-vpwindow[2];      
	dopN = -dopN; // dopN'
      
	if (dopN > 0) dopN = -dopN; // warning MESSAGE// warning MESSAGE

	a = 2*dopN/(diffU * (dopN - B));
	b = 2*dopN/(diffV * (dopN - B));
	c = -1.0/(dopN - B);

	float zMin = (dopN + F)*c;
	if (zMin >= 0) {  // warning MESSAGE
	    F = -1.0e-8-dopN;
	    zMin = -1.0e8;
	}
     
	e = 1/(1+zMin);
	f = -zMin * e;
	
	g = cop[0] - shU * cop[2];
	h = cop[1] - shV * cop[2];
  
	ProjMat.Set(     a,   0.0,         0.0,    0.0,
		       0.0,     b,         0.0,    0.0,
		     a*shU, b*shV,         e*c,     -c,
		       -a*g,  -b*h, e*c*dopN+f, -c*dopN);  // dopN = -cop[2]!
    }

    for (int jj=0; jj<4; jj++)
	for (int kk=0; kk<4; kk++)
	    projData[jj*4+kk] = ProjMat.m[jj][kk];

    glMatrixMode(GL_PROJECTION);
    glLoadMatrixf(projData);

    if (aFile == 1) {
	theFile << "StartImage\n";
	theFile << "VRP " << vrp[0] << " " << vrp[1] << " " 
	    << vrp[2] << " " << endl;
	theFile << "VPN " << vpn[0] << " " << vpn[1] << " " 
	    << vpn[2] << " " << endl;
	theFile << "VUV " << vuv[0] << " " << vuv[1] << " " 
	    << vuv[2] << " " << endl;
	theFile << "COP " << cop[0] << " " << cop[1] << " " 
	    << cop[2] << " " << endl;
	
	theFile << "PROJECTIONMODE " << projection_mode << endl;
	theFile << "VPWINDOW " << vpwindow[0] << " " << vpwindow[1] << " "
	    << vpwindow[2] << " " << vpwindow[3] << " " << endl;
	theFile << "PLANES " << planedist[0] << " " << planedist[2] << "\n";
	theFile << "PORTWINDOW " << portwindow[0] << " " << portwindow[1] << " "
	    << portwindow[2] << " " << portwindow[3] << " " << endl;

    } 
 
    // done
    return 0;
}


int 
OpenGLRenderer::doneImage(void)
{

 #ifdef _WIN32
	wglMakeCurrent(theHDC, theHRC);
    SwapBuffers(theHDC);
#endif
	glFlush();

    if (aFile == 1) {
	theFile << "DoneImage\n";
    }
    
    if (count != -1) {
	count++;

	this->saveBmpImage();
    }



    
    return 0;
}

int 
OpenGLRenderer::drawLine(const Vector &pos1, const Vector &pos2, 
		       float V1, float V2)
{
	//glLineWidth(10);
    glBegin(GL_LINES);
    float r, g, b;
    r = theMap->getRed(V1);
    g = theMap->getGreen(V1);
    b = theMap->getBlue(V1);

    glColor3f(r,g,b);

    glVertex3f(pos1(0),pos1(1),pos1(2));
    
    if (aFile == 1) {
	theFile << "Line\n" << pos1(0) << " " << pos1(1) << " " << pos1(2) 
	    << " " << r << " " << g << " " << b << " " << endl;
    }

    r = theMap->getRed(V2);
    g = theMap->getGreen(V2);
    b = theMap->getBlue(V2);

    glColor3f(r,g,b);

    glVertex3f(pos2(0),pos2(1),pos2(2));
    glEnd();

    if (aFile == 1) {
	theFile << pos2(0) << " " << pos2(1) << " " << pos2(2) << " " << r 
	    << " " << g << " " << b << " " << endl;
    }

    return 0;  
}


int 
OpenGLRenderer::drawLine(const Vector &end1, const Vector &end2, 
			 const Vector &rgb1, const Vector &rgb2)
{
    glBegin(GL_LINES);
    float r, g, b;
    r = rgb1(0);
    g = rgb1(1);
    b = rgb1(2);
  
    if (aFile == 1) {
	theFile << "Line\n" << end1(0) << " " << end1(1) << " " << end1(2) 
	    << " " << r << " " << g << " " << b << " " << endl;
    }
    glColor3f(r,g,b);
    glVertex3f(end1(0),end1(1),end1(2));

    r = rgb2(0);
    g = rgb2(1);
    b = rgb2(2);
  
    if (aFile == 1) {
	theFile << end2(0) << " " << end2(1) << " " << end2(2) << " " << r 
	    << " " << g << " " << b << " " << endl;
    }
    glColor3f(r,g,b);

    glVertex3f(end2(0),end2(1),end2(2));
    glEnd();

    return 0;
}



int 
OpenGLRenderer::drawPolygon(const Matrix &pos, const Vector &data)

{
#ifdef _G3DEBUG
  if (pos.noCols() != 3) {
    g3ErrorHandler->warning("OpenGLRenderer::drawPolygon - matrix needs 3 cols\n");
    return -1;
  }
  if (pos.noRows() != data.Size()) {
    g3ErrorHandler->warning("OpenGLRenderer::drawPolygon - matrix & vector incompatable\n");
    return -1;
  }
#endif

	double posX, posY, posZ, value;
	float r,g,b;

    glBegin(GL_POLYGON);
    int numRows = pos.noRows();
    for (int i=0; i<numRows; i++) {
      posX = pos(i,0);
      posY = pos(i,1);
      posZ = pos(i,2);
      value = data(i);
      r = theMap->getRed(value);
      g = theMap->getGreen(value);
      b = theMap->getBlue(value);      

    if (aFile == 1) {
	theFile << posX << " " << posY << " " << posZ << " " << r 
	    << " " << g << " " << b << " " << endl;
    }
      glColor3f(r,g,b);
      glVertex3f(posX, posY, posZ);
    }

    glEnd();

    return 0;
}



int 
OpenGLRenderer::drawGText(const Vector &pos, char *text, int length)
{
    MYPOINT *point;

    // add POINTs to the FACE  
    int size = pos.Size();
    float x,y,z;
    if (size == 1) {
	x = pos(0);
	y = 0;
	z = 0;
    } else if (size == 2) {
	x = pos(0);
	y = pos(1);
	z = 0;
    } else {
	x = pos(0);
	y = pos(1);
	z = pos(2);
    }  

    return 0;
}

int 
OpenGLRenderer::drawLText(const Vector &pos, char *text, int length)
{
    int size = pos.Size();
    float x,y;
    if (size == 1) {
	x = pos(0);
	y = 0;
    } else if (size == 2) {
	x = pos(0);
	y = pos(1);
    } else {
	x = pos(0);
	y = pos(1);
    }  

//  theDevice->drawText(x,y, text, length);
  
    return 0;
}


int 
OpenGLRenderer::setVRP(float x, float y, float z)
{
  vrp[0] = x;
  vrp[1] = y;
  vrp[2] = z;

  return 0;
}

int 
OpenGLRenderer::setVPN(float x, float y, float z)
{
  vpn[0] = x;
  vpn[1] = y;
  vpn[2] = z;

  return 0;
}

int 
OpenGLRenderer::setVUP(float x, float y, float z)
{
  vuv[0] = x;
  vuv[1] = y;
  vuv[2] = z;

  return 0;
}

int 
OpenGLRenderer::setViewWindow(float umin, float umax, float vmin, float vmax)
{
  if (umin > umax || vmin > vmax) {
      cerr << "OpenGLRenderer::setViewWindow() - invalid window ";
      cerr << umin << " "<< umax << " "<< vmin << " "<< vmax << endl;
      return -1;
  }

  vpwindow[0] = umin;
  vpwindow[1] = umax;
  vpwindow[2] = vmin;
  vpwindow[3] = vmax;

  return 0;
}

int 
OpenGLRenderer::setPlaneDist(float anear, float afar) 
{

   if ((anear < 0.0) || (afar < 0.0)) {
      cerr << "OpenGLRenderer::setPlaneDist() - invalid planes";
      cerr << anear << " " << afar << endl;
      return -1;
  }
  
  planedist[0] = anear;
  planedist[2] = afar;

  return 0;
}

int 
OpenGLRenderer::setProjectionMode(int newMode)
{
  projection_mode = newMode;
  return 0;
}

int 
OpenGLRenderer::setFillMode(int newMode)
{
//  theScan->setFillMode(newMode);
  return 0;
}

// eye location
int 
OpenGLRenderer::setPRP(float u, float v, float n){
  cop[0] = u;
  cop[1] = v;
  cop[2] = n;

  return 0;
}
    
int 
OpenGLRenderer::setPortWindow(float left, float right, 
			     float bottom, float top)
{
  if (left < -1 || right > 1 || bottom < -1 || top > 1
      || left > right || bottom > top) {
      
      cerr << "OpenGLRenderer::setPortWindow() - bounds invalid ";
      cerr << left << " "<< right << " "<< bottom << " "<< top << endl;
      return -1;
  }

  portwindow[0] = left;
  portwindow[1] = right;
  portwindow[2] = bottom;
  portwindow[3] = top;

  return 0;
}

int 
OpenGLRenderer::saveBmpImage(void)
{
    long i,j, width;
/*
    // get dimensions of the window
    glGetIntegerv(GL_VIEWPORT, viewport);
	*/  
    // create the file name 'bmpFileName$count.BMP'
    char fileName[MAX_FILENAMELENGTH+14];
    char intName[10];
	
    strcpy(fileName, theBmpFileName);
	sprintf(intName,"%d",count);

    strcat(fileName,&intName[1]);        
    strcat(fileName,".BMP");

    // open the file
    FILE *fp;
    if ((fp = fopen(fileName,"wb")) == NULL) {
	g3ErrorHandler->warning("OpenGLRenderer::saveBmpImage() - %s %s\n",
				"could not open file named", fileName);
	count = -1;
	return -1;
    }	
//cerr << "saving: " << fileName << endl;
	
	int bitsize = (info.bmiHeader.biWidth *
        	   info.bmiHeader.biBitCount + 7) / 8 *
		  abs(info.bmiHeader.biHeight);
	int infosize;
	infosize = sizeof(BITMAPINFOHEADER);
    switch (info.bmiHeader.biCompression)
	{
	case BI_BITFIELDS :
            infosize += 12; /* Add 3 RGB doubleword masks */
            if (info.bmiHeader.biClrUsed == 0)
	      break;
	case BI_RGB :
            if (info.bmiHeader.biBitCount > 8 &&
        	info.bmiHeader.biClrUsed == 0)
	      break;
	case BI_RLE8 :
	case BI_RLE4 :
            if (info.bmiHeader.biClrUsed == 0)
              infosize += (1 << info.bmiHeader.biBitCount) * 4;
	    else
              infosize += info.bmiHeader.biClrUsed * 4;
	    break;
	}

    int size = sizeof(BITMAPFILEHEADER) + infosize + bitsize;

    // check the bit map header info has been created, if not create one
/*    if (info == 0) {
	if ((info = (BITMAPINFO *)malloc(sizeof(BITMAPINFOHEADER))) < 0) {
	    g3ErrorHandler->warning("OpenGLRenderer::saveBmpImage() - %s\n",
				    "out of memory creating BITMAPINFO struct");
	    count = -1;
	    return -2;
	}
	info->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	info->bmiHeader.biPlanes = 1;
	info->bmiHeader.biBitCount = 24;
	info->bmiHeader.biCompression = BI_RGB;
	info->bmiHeader.biXPelsPerMeter = 2952;
	info->bmiHeader.biYPelsPerMeter = 2952;
	info->bmiHeader.biClrUsed = 0;
	info->bmiHeader.biClrImportant = 0;
    }
    
    // determine the number of bits needed to save the image
    width = viewport[2]*3;
    width = (width+3) & ~3;
    bitsize = width * viewport[3];

    // check the bits pointeris of correct size, if not
    // delete the old and create a new one
    if (bitsize != currentBitSize) {
	if (currentBitSize != 0)
	    free (bits);
	if ((bits = (GLubyte *)calloc(bitsize, 1)) == 0) {
	    g3ErrorHandler->warning("OpenGLRenderer::saveBmpImage() - %s\n",
				    "out of memory creating BITMAPINFO struct");
	    count = -1;
	    free (info);
	    info = 0;
	    return -3;	    
	}
	currentBitSize = bitsize;
    }

    // set the info for the bit map header
    info->bmiHeader.biWidth = viewport[2];
    info->bmiHeader.biHeight = viewport[3];    
    info->bmiHeader.biSizeImage = currentBitSize;        
*/
    // read the pixels from the frame buffer
    glFinish();
    glPixelStorei(GL_PACK_ALIGNMENT, 4);
    glPixelStorei(GL_PACK_ROW_LENGTH, 0);
    glPixelStorei(GL_PACK_SKIP_ROWS, 0);
    glPixelStorei(GL_PACK_SKIP_PIXELS, 0);

  if (bits == 0) cerr << "BITS ZERO\n";

    glReadPixels(0, 0, info.bmiHeader.biWidth, info.bmiHeader.biHeight,
		GL_BGR_EXT, GL_UNSIGNED_BYTE, bits);
	

   currentBitSize = info.bmiHeader.biWidth * info.bmiHeader.biHeight;
    // create a header for the BMP file
    BITMAPFILEHEADER header;
    header.bfType      = 'MB'; /* Non-portable... sigh */
    header.bfSize      = size;
    header.bfReserved1 = 0;
    header.bfReserved2 = 0;
    header.bfOffBits   = sizeof(BITMAPFILEHEADER) + infosize;

    if (fwrite(&header, 1, sizeof(BITMAPFILEHEADER), fp) < sizeof(BITMAPFILEHEADER))
	{
    // write the header to the file
//    if (fwrite(&header, 1, sizeof(BITMAPFILEHEADER), fp) < sizeof(BITMAPFILEHEADER)) {
	    g3ErrorHandler->warning("OpenGLRenderer::saveBmpImage() - %s\n",
				    "failed to write BITMAPHEADER");
	    fclose(fp);
	    return -4;
	}
    if (fwrite(&info, 1, infosize, fp) < infosize)
        {
    // write the bit map information to the file
//    if (fwrite(&info, 1, sizeof(BITMAPINFOHEADER), fp) < sizeof(BITMAPINFOHEADER)) {
	    g3ErrorHandler->warning("OpenGLRenderer::saveBmpImage() - %s\n",
				    "failed to write BITMAPINFOHEADER");
	    fclose(fp);
	    return -5;	    
	}    
    if (fwrite(bits, 1, bitsize, fp) < bitsize)
        {
    // now we write the bits
    //if (fwrite(bits, 1, currentBitSize, fp) < currentBitSize) {
	g3ErrorHandler->warning("OpenGLRenderer::saveBmpImage() - %s\n",
				    "failed to write BITMAPINFOHEADER");
	fclose(fp);
	return -6;	
    }        

    // if get here we are done .. close file and return ok
    fclose(fp);
    return 0;
}
/*   
char itoc(int x)
{
 if (x == 1) return '1';
 if (x == 2) return '2';
 if (x == 3) return '3';
 if (x == 4) return '4';
 if (x == 5) return '5';
 if (x == 6) return '6';
 if (x == 7) return '7';
 if (x == 8) return '8';
 if (x == 9) return '9';
 return '0';
}

void
itoa(int x, char *str)
{
  int y=x;
  while (y >= 10) 
    y = y/10;
  cerr << y << " ";
  str[0] = itoc(y);
  str[1] = '\0';
  if (x >= 10) {
    int z = x/10;
    z = x - 10*z;
    itoa(z,&str[1]);
  }
}
*/
