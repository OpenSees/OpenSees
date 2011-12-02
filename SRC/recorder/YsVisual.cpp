// YsVisual.cpp: implementation
//
//////////////////////////////////////////////////////////////////////
#include "YsVisual.h"
#include <InelasticYS2DGNL.h>

#ifdef _GLX
#include <OpenGLRenderer.h>
#endif
#ifdef _X11
#include <X11Renderer.h>
#endif

#include <Renderer.h>
#include <PlainMap.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

YsVisual::YsVisual(Element *theEle, const char *title, double scale,
                   int xLoc, int yLoc, int width, int height)
  :Recorder(RECORDER_TAGS_YsVisual),
   pView(0), theMap(0), displayMode(0)
{
  const char *argv[1] = {"ysVisual"};
  Information eleInfo(1);
  Response *theResponse = 0;
  
  ysEle = theEle;
  
  theResponse = theEle->setResponse(argv, 1, eleInfo);
  
  if(theResponse)
    displayMode = theResponse->getResponse();
  else
    opserr << "YsVisual only works on elements derived from InelasticYS2DGNL\n";
  
#ifdef _GLX
  theMap = new PlainMap();
  pView =  new OpenGLRenderer(title, xLoc, yLoc, width, height, *theMap);
#endif
  
#ifdef _X11
  if(!pView) //_GLX and _X11 both defined
    {
      theMap = new PlainMap();
      pView =  new X11Renderer(title, xLoc, yLoc, width, height, *theMap);
    }
#endif	
  
  if(pView)
    {
      pView->setVRP(0.0, 0.0, 0.0);
      pView->setVPN(0.0, 0.0, 1.0);
      pView->setVUP(0.0, 1.0, 0.0);
      pView->setFillMode("wire");             // wire mode
      pView->setPlaneDist(1.0, -1.0);
      pView->setPRP(0.0, 0.0, 10.0);
      pView->setPortWindow(-1, 1, -1, 1);  // use the whole window
      
      pView->setViewWindow(-scale, scale, -scale, scale);
      
		pView->clearImage();
		pView->doneImage();
	}
	else
		opserr << "WARNING: YsVisual::createView - Renderer not available\n";

}

YsVisual::~YsVisual()
{
	if(theMap)
		delete theMap;
	if(pView)
		delete pView;
}

int  YsVisual::record(int cTag, double d)
{
int res = 0;

	if(pView)
	{
    	pView->clearImage();
    	pView->startImage();
    	ysEle->displaySelf(*pView, displayMode, 1);
    	pView->doneImage();
	}
	else 
		res = -1;

	return res;
}

int YsVisual::playback(int cTag)
{
  return 0;
}

int 
YsVisual::restart(void)
{
  return 0;
}
