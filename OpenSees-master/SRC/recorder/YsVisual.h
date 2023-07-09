// YsVisual.h: interface
//
//////////////////////////////////////////////////////////////////////

#if !defined(YSVISUAL)
#define YSVISUAL

#include <Recorder.h>
#include <Element.h>
#include <Response.h>

class ColorMap;
class Renderer;

class YsVisual : public Recorder
{
public:
	virtual ~YsVisual();

	YsVisual(Element* theEle, const char *title, double scale,
	         int xLoc, int yLoc, int width, int height);
	

    int record(int cTag, double d);
    int playback(int cTag);
    int restart(void);    

  private:

  int      displayMode;
  Element  *ysEle;
  Renderer *pView;
  ColorMap *theMap;
};

#endif // !defined
