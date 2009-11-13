#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>
#include <X11/Intrinsic.h>
#include <X11/StringDefs.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <X11/Xresource.h>
#include "icon"


/* number of colors you want plus paul's original 16!!! */
#define NUM_COLOR_CELLS  (64+16)

struct mse_event_rec
  {
    int xclicked;
    int yclicked;
    int buttonNumber;
  };

Display *dpy = NULL;
static Window theWindow, iconWindow;
static Cursor cursorPen;
static XFontStruct *textFontInfo=NULL;
static GC colorPenGC, fillGC, fountainPenGC, calligPenGC;
static XGCValues gcv;

/* hmt's default */
/* int windowh = 950,windoww = 1300; */
/* hmt workshop debug */
/*int windowh = 615,windoww = 650; */
/* int vfx = 0, vfy = 1097;  */


/* int windowh = 896,windoww = 952;  */  /*  PAUL's DEFAULT !!!*/
/* int windowh = 706,windoww = 690;  */  /* small */
/* int windowh = 826,windoww = 890;  */
   int windowh = 950,windoww =1022;

int vfx = 0, vfy = 0;




int char_spacing = 6;
int line_spacing = 15;
int text_screenh = 150;
int left_margin = 000;
char Line[10][90];
int active_line = 0;
int top_line = 0;
int num_lines = 10;

static XPoint colorPenPosition, fountainPenPosition;
static XColor colorCells[NUM_COLOR_CELLS];
static Colormap cmap;
static Pixmap IconPix, dimBorder;
static int ForeColor, BackColor;

static struct mse_event_rec mouseRecord;
static int lineEntered=0;
static int scroll=0;
static char input_line[81];
static int input_line_idx=0;

extern void xmainloop();
extern int Quit();
extern void error();


/* checkerboard used in dimBorder */
#define check_width 16
#define check_height 16
static char check_bits[] = {
    0x55, 0x55, 0xaa, 0xaa, 0x55, 0x55, 0xaa, 0xaa, 0x55, 0x55, 0xaa, 0xaa,
    0x55, 0x55, 0xaa, 0xaa, 0x55, 0x55, 0xaa, 0xaa, 0x55, 0x55, 0xaa, 0xaa,
    0x55, 0x55, 0xaa, 0xaa, 0x55, 0x55, 0xaa, 0xaa };

#define MOUSE 0
#define LINE  1


void Random_Seed(unsigned seed);
int  Random_Integer(int lb, int ub);


int 
x_num_colors(void)
{
  return(NUM_COLOR_CELLS);
}


   
void resetXPRINT()  
  {
    fountainPenPosition.x = left_margin;
    fountainPenPosition.y = windowh - text_screenh;
    XDrawImageString(dpy,theWindow,fountainPenGC,
		     fountainPenPosition.x,fountainPenPosition.y,"_",1);    
    XFlush(dpy);
  }

void SetUpEnv () 
  {
    char *argvec[1];
    int argcount=1,i,j;  
    XSizeHints szhint;
    long atol();
    XSetWindowAttributes attributes;
    XEvent event;
    Arg arg[5];
    char *font;


/*  argvec[0]="Nekton";    pff 7/29/92  for photo ops.... */
/*  argvec[0]="      ";                                   */
    argvec[0]="Nekton";
    for (i=0;i<80;i++)
      input_line[i]=' ';
    for (i=0;i<10;i++)
      for (j=0;j<90;j++)
	Line[i][j]= ' ';
    
    if ((dpy = XOpenDisplay(NULL)) == NULL) 
	 error("Can't open display \n Did you setenv DISPLAY to your localHost \n and xhost + on the local display? \n");

    szhint.width = szhint.max_width = szhint.min_width = windoww;
    szhint.height = szhint.max_height = szhint.min_height = windowh;
    szhint.x = 10;
    szhint.y = 10;
    szhint.flags = PSize|PMinSize|PMaxSize|USPosition; 
    
    ForeColor = WhitePixel (dpy, DefaultScreen(dpy));
    BackColor = BlackPixel (dpy, DefaultScreen(dpy));
    dimBorder = XCreatePixmapFromBitmapData(dpy,DefaultRootWindow(dpy),
				      check_bits,check_width,check_height,
				      ForeColor,BackColor,
				      DefaultDepth(dpy,DefaultScreen(dpy)));
    theWindow = XCreateSimpleWindow (dpy,RootWindow(dpy,DefaultScreen(dpy)),
				     vfx,vfy,szhint.width,szhint.height,2,
				     ForeColor,BackColor);  

    if (!theWindow) error("Can't Open Window\n");
    XSelectInput(dpy,theWindow,KeyPressMask|ButtonPressMask);


    IconPix = XCreateBitmapFromData (dpy, DefaultRootWindow(dpy),
				     icon_bits,icon_width,icon_height);

    /*                         pff 7/29/92  for photo ops.... */
    /* XSetStandardProperties (dpy,theWindow,"      ",NULL,IconPix, */
    /* XSetStandardProperties (dpy,theWindow,"Nekton",NULL,IconPix, */
    XSetStandardProperties (dpy,theWindow,"Nekton",NULL,IconPix,
			    argvec,argcount,&szhint);
    XMapWindow(dpy,theWindow);
    XRaiseWindow(dpy,theWindow);
    XSync(dpy,False);

    cursorPen = XCreateFontCursor (dpy,XC_plus);
    XDefineCursor(dpy, theWindow, cursorPen);
    XFreeCursor(dpy,cursorPen);
    
    gcv.line_style = LineSolid;
    gcv.foreground = WhitePixel (dpy, DefaultScreen(dpy));
    gcv.background = BlackPixel (dpy, DefaultScreen(dpy));


    colorPenGC = XCreateGC(dpy,theWindow,
			   GCLineStyle|GCForeground|GCBackground,&gcv);

    fountainPenGC=XCreateGC(dpy,theWindow,GCForeground|GCBackground,&gcv);

    gcv.fill_style = FillSolid;
    fillGC = XCreateGC(dpy,theWindow,GCFillStyle|GCForeground,&gcv);

    /* standard fonts */
    /* Extra Helvetica font sizes big(--24) to small(--8) */
    /*
      "-Adobe-Helvetica-Medium-R-Normal--24-240-75-75-P-130-ISO8859-1";
      "-Adobe-Helvetica-Medium-R-Normal--18-180-75-75-P-98-ISO8859-1";
      "-Adobe-Helvetica-Medium-R-Normal--14-140-75-75-P-78-ISO8859-1";
      "-Adobe-Helvetica-Medium-R-Normal--12-120-75-75-P-67-ISO8859-1";
      "-Adobe-Helvetica-Medium-R-Normal--10-100-75-75-P-57-ISO8859-1";
      "-Adobe-Helvetica-Medium-R-Normal--8-80-75-75-P-47-ISO8859-1";

      "-*-*-*-*-*--*-*-*-*-*-*-*-*";
    */

    /* use 12pt for big window */
    if (windoww >= 700) 
      {
	font = "-Adobe-Helvetica-Medium-R-Normal--12-120-75-75-P-67-ISO8859-1";
        textFontInfo = XLoadQueryFont(dpy,font);
	if (!textFontInfo)
	  {
	    font = "-*-Helvetica-Medium-R-Normal--12-120-75-75-P-67-ISO8859-1";
	    textFontInfo = XLoadQueryFont(dpy,font);
	  }
      }
    /* use 8pt for small window */
    else 
      {
	font = "-Adobe-Helvetica-Medium-R-Normal--8-80-75-75-P-47-ISO8859-1";
	textFontInfo = 	XLoadQueryFont(dpy,font);
	if (!textFontInfo)
	  {
	    font = "-*-Helvetica-Medium-R-Normal--8-80-75-75-P-47-ISO8859-1";
	    textFontInfo = XLoadQueryFont(dpy,font);
	  }
      }

    /* no Helvetica-Medium-R-Normal 12||8 ... look for something */
    if (!textFontInfo) 
      {
	font = "-*-Helvetica-Medium-R-Normal--*-*-*-*-*-*-*-*";
	textFontInfo = 	XLoadQueryFont(dpy,font);
      }

    /* ok anything in Medium-R-Normal */
    if (!textFontInfo) 
      {
	font = "-*-*-Medium-R-Normal--*-*-*-*-*-*-*-*";
	textFontInfo = 	XLoadQueryFont(dpy,font);
      }

    /* really desperate now */
    if (!textFontInfo) 
      {
	font = "-*-*-*-*-*--*-*-*-*-*-*-*-*";
	textFontInfo = 	XLoadQueryFont(dpy,font);
      }

    if (!textFontInfo) 
      {error("XLoadQueryFont error in SetUpEnv\n");}

    printf("%s\n",font);

    gcv.font = textFontInfo->fid;

    calligPenGC= XCreateGC(dpy,theWindow,
			   GCFont|GCForeground|GCBackground,&gcv);
    resetXPRINT();

#ifdef TITAN
    while(1){
      XDrawImageString(dpy,theWindow,fountainPenGC,
		       fountainPenPosition.x,fountainPenPosition.y,
		       "Please Hit Any Key to Begin",27);

      if(XCheckMaskEvent(dpy,KeyPressMask|ButtonPressMask,&event)==True)
	break;
    }
      XDrawImageString(dpy,theWindow,fountainPenGC,
		       fountainPenPosition.x,fountainPenPosition.y,
		       "                           ",27);
#endif    

    XFlush(dpy);
  }

#if (CRAY|TITAN)
void GETWIN(width,height)
#else
void getwin_(width,height)
#endif
     float *width,*height;
{
  *width  = (float) windoww;
  *height = (float) windowh;
}

#if (CRAY|TITAN)
void XPLANES(nplanes)
#else
void xplanes_(nplanes)
#endif
     int *nplanes;
  {
    *nplanes=XDisplayPlanes(dpy,XDefaultScreen(dpy));
  }

#if (CRAY|TITAN)     
void XWAITM(button,xmouse,ymouse)
#else
void xwaitm_(button,xmouse,ymouse)
#endif
     int *button;
     float *xmouse,*ymouse;
  {
    while(1){
      xmainloop(MOUSE);
      if (mouseRecord.yclicked<=(windowh-text_screenh))
	/* y position is not inside text region,
	 * so mouse click was valid. */
	break;
      }
    *ymouse = (float)fabs((double)(mouseRecord.yclicked-
				   (windowh-text_screenh)));
    *xmouse = (float)mouseRecord.xclicked;
    *button = mouseRecord.buttonNumber;
    XFlush(dpy);
  }

#if (CRAY|TITAN)
void BEEP()
#else
void beep_()
#endif
  {
    XBell(dpy,0);
    XFlush(dpy);
  }

#ifdef CRAY
void GETS(line)
     char *line;
#else
#ifdef TITAN
void GETS(line)
     char **line;
#else
void gets_(line)
     char *line;
#endif
#endif
  {
    int i;
    xmainloop(LINE);
    for (i=0;i<strlen(input_line)&&i<80;i++)
#ifdef TITAN      
      line[0][i]=input_line[i];
#else
      line[i]=input_line[i];
#endif
    for (i=0;i<80;i++)
      input_line[i]=' ';
    lineEntered=0;
    input_line_idx=0;
  }
     
#if (CRAY|TITAN)
void XDRAWTO(x,y)
#else
void xdrawto_(x,y)
#endif
     float *x,*y;
  {
    XDrawLine(dpy,theWindow,colorPenGC,
	      colorPenPosition.x,colorPenPosition.y, 
	      (short) *x, (short) *y);
    colorPenPosition.x = (short) *x;
    colorPenPosition.y = (short) *y;
    XFlush(dpy);
  }
     
#if (CRAY|TITAN)
void XMOVE(x,y)
#else
void xmove_(x,y)
#endif
     float *x,*y;
  {
    colorPenPosition.x = (short) *x;
    colorPenPosition.y = (short) *y;
    XFlush(dpy);
  }
     
#if (CRAY|TITAN)
void XCMAP(red,green,blue,num)
#else
void xcmap_(red,green,blue,num)
#endif
     float red[], green[], blue[];
     int *num;
  {
    int i;
    cmap = XDefaultColormap(dpy,XDefaultScreen(dpy));     
    for (i=0;i<=*num;i++)
	{
	  colorCells[i].red = (short) red[i];
	  colorCells[i].green = (short) green[i];
	  colorCells[i].blue = (short) blue[i];
	  colorCells[i].flags = DoRed|DoGreen|DoBlue;
	  XAllocColor(dpy,cmap,&colorCells[i]);
	}
    XSync(dpy,False);
  }


#if (CRAY|TITAN)
void hmt_XCMAP(red,green,blue,num)
#else
void hmt_xcmap_(red,green,blue,num)
#endif
     float red[], green[], blue[];
     int *num;
  {
    int i;
    int r,g,b;
    int ct;

    cmap = XDefaultColormap(dpy,XDefaultScreen(dpy));     

    /* paul's 16 colors */
    for (i=0;i<=*num;i++)
	{
	  if (*num==1)
	    {
	      colorCells[i].red = (short) 65535;
	      colorCells[i].green = (short) 65535;
	      colorCells[i].blue = (short) 65535;
	    }
	  else
	    {
	      colorCells[i].red = (short) red[i];
	      colorCells[i].green = (short) green[i];
	      colorCells[i].blue = (short) blue[i];
	    }
	  colorCells[i].flags = DoRed|DoGreen|DoBlue;
	  printf("(%4u,%4u,%4u)\n",colorCells[i].red,
		 colorCells[i].green,colorCells[i].blue);
	  XAllocColor(dpy,cmap,&colorCells[i]);
	}

    /* extras */
    Random_Seed(1729);
    for (i=*num+1;i<NUM_COLOR_CELLS;i++, ct++)
	{
	  r = Random_Integer(16383,65535);
	  g = Random_Integer(00000,65537);
	  b = Random_Integer(00000,65535);
	  colorCells[i].red = (short) r;
	  colorCells[i].green = (short) g;
	  colorCells[i].blue = (short) b;
	  colorCells[i].flags = DoRed|DoGreen|DoBlue;
	  /*
	  printf("(%4u,%4u,%4u)\n",colorCells[i].red,
		 colorCells[i].green,colorCells[i].blue);
	  */
	  XAllocColor(dpy,cmap,&colorCells[i]);
	}
    XSync(dpy,False);
    /*
    printf("(sizeof(us)=%d\n",sizeof(unsigned short));
    */
  }




void
Random_Seed(unsigned seed)
{
  void srand(unsigned seed);

  srand(seed);
}


int 
Random_Integer(int lb, int ub)
{
  int rand(void);
  double d;
  int k;


  d = rand()/((double) RAND_MAX + 1);
  k = (int) (d*(ub-lb+1));

  return(lb+k);
}
	     


#if (CRAY|TITAN)
void HMT_XSETFILL(icolor)
#else
void hmt_xsetfill_(icolor)
#endif    
     int *icolor;
  {
    /*    printf("icolor=%d, ng=%d\n",*icolor,*ng%NUM_COLOR_CELLS); */
    /*    gcv.foreground = colorCells[*icolor].pixel; */
    /*    printf("xsetfill() :: icolor=%d\n",*icolor); */

    gcv.foreground = colorCells[*icolor%NUM_COLOR_CELLS].pixel;
    XChangeGC(dpy,fillGC,GCForeground,&gcv);
    XFlush(dpy);
  }



#if (CRAY|TITAN)
void XSETFILL(icolor)
#else
void xsetfill_(icolor)
#endif    
     int *icolor;
  {
    gcv.foreground = colorCells[*icolor].pixel;
    XChangeGC(dpy,fillGC,GCForeground,&gcv);
    XFlush(dpy);
  }
     
#if (CRAY|TITAN)
void XSETPEN(icolor)
#else 
void xsetpen_(icolor)
#endif
     int *icolor;
  {
    /* printf("xsetpen() :: icolor=%d\n",*icolor); */
    gcv.foreground = colorCells[*icolor%NUM_COLOR_CELLS].pixel;
    XChangeGC(dpy,colorPenGC,GCForeground,&gcv);
    XFlush(dpy);
  }
     
#if (CRAY|TITAN)
void XPOLY(xfill,yfill,nvertices,boundary)
#else
void xpoly_(xfill,yfill,nvertices,boundary)
#endif
     float *xfill, *yfill;
     int *nvertices, *boundary;
  {
    int i;
    XPoint *vertices;
    vertices = (XPoint *) calloc(*nvertices,sizeof(XPoint));
    for (i=0;i<*nvertices;i++) 
	{
	  vertices[i].x = (short) xfill[i];
	  vertices[i].y = (short) yfill[i];
	}
    XFillPolygon(dpy,theWindow,
		 fillGC,vertices,*nvertices,Nonconvex,CoordModeOrigin);

    if (*boundary)  
	{
	  XDrawLines(dpy,theWindow,
		     colorPenGC,vertices,*nvertices,CoordModeOrigin);
	  XDrawLine(dpy,theWindow,colorPenGC,
		    vertices[*nvertices-1].x,vertices[*nvertices-1].y,
		    vertices[0].x,vertices[0].y);
	}
    XFlush(dpy);
  }

#if (CRAY|TITAN)
void XCLEAR()
#else
void xclear_()
#endif
  { 
    int i;
    XPoint corner[4];
    corner[0].x=0;
    corner[0].y=0;
    corner[1].x=windoww;
    corner[1].y=0;
    corner[2].x=corner[1].x;
    corner[2].y=windowh-text_screenh-line_spacing-line_spacing;
    corner[3].x=0;
    corner[3].y=corner[2].y;
    i=0;
#if (CRAY|TITAN)
    XSETFILL(&i);
#else
    xsetfill_(&i);
#endif
    XFillPolygon(dpy,theWindow,fillGC,corner,4,Convex,CoordModeOrigin);
    XFlush(dpy);
  }

#if (CRAY|TITAN)
void XCLEARA()
#else
void xcleara_()
#endif
  {
    int i;
    XPoint corner[4];
    corner[0].x=0;
    corner[0].y=0;
    corner[1].x=windoww-(int)((double)windoww/1.3);
    corner[1].y=0;
    corner[2].x=corner[1].x;
    corner[2].y=windowh-text_screenh-line_spacing-line_spacing;
    corner[3].x=0;
    corner[3].y=corner[2].y;
    i=0;
#if (CRAY|TITAN)
    XSETFILL(&i);
#else
    xsetfill_(&i);
#endif
    XFillPolygon(dpy,theWindow,fillGC,corner,4,Convex,CoordModeOrigin);
    XFlush(dpy);
  }     

#ifdef TITAN
void XPUTXT(x,y,fontsize,string,strlength)
     float *x,*y,*fontsize;
     char **string;
     int *strlength;
#else
#ifdef CRAY
void XPUTXT(x,y,fontsize,string,strlength)
     float *x,*y,*fontsize;
     char *string;
     int *strlength;
#else
void xputxt_(x,y,fontsize,string,strlength)
     float *x,*y,*fontsize;
     char *string;
     int *strlength;
#endif
#endif
  {
#ifdef TITAN
    XDrawString(dpy,theWindow,
		     calligPenGC,(int)*x,(int)*y,*string,*strlength);
#else
    XDrawString(dpy,theWindow,
		     calligPenGC,(int)*x,(int)*y,string,*strlength);
#endif
    XFlush(dpy);
  }
     
void advanceLine()
  {
    int i;
    if (scroll)	{
      top_line = (top_line+1)%num_lines;
      active_line = (active_line+1)%num_lines;
      for (i=0;i<num_lines;i++)  
	XDrawImageString(dpy,theWindow,fountainPenGC,
			 left_margin,windowh-text_screenh+i*line_spacing,
			 Line[(top_line+i)%num_lines],80);
      fountainPenPosition.x=left_margin;
      XDrawImageString(dpy,theWindow,fountainPenGC,
		       fountainPenPosition.x,fountainPenPosition.y,
		       "_                                                                                                                                                                                                                                                              ",150);
      XFlush(dpy);
      return;
    }
    active_line++;
    if (active_line == num_lines) {
      active_line--;
      scroll=1;
      advanceLine();
      return;
    }      
    fountainPenPosition.x=left_margin;
    if (fountainPenPosition.y < windowh)
      fountainPenPosition.y += line_spacing;
    if (fountainPenPosition.y >= windowh)
      fountainPenPosition.y =  windowh - text_screenh;
    XDrawImageString(dpy,theWindow,fountainPenGC,
		     fountainPenPosition.x,fountainPenPosition.y,
		     "_                                                                                                                                                                                                                                                              ",150);
    XFlush(dpy); 
  }
     
#ifdef TITAN
void PUTS(message,nchars)
     char **message;
     int *nchars;
#else
#ifdef CRAY
void PUTS(message,nchars)
     char *message;
     int *nchars;
#else
void puts_(message,nchars)
     char *message;
     int *nchars;
#endif
#endif
  {
    int i;
#ifdef TITAN
    XDrawImageString(dpy,theWindow,fountainPenGC,
		     fountainPenPosition.x,fountainPenPosition.y,
		     *message,*nchars);
    for (i=0;i < *nchars;i++)
      Line[active_line][i]= message[0][i];
#else
    XDrawImageString(dpy,theWindow,fountainPenGC,
		     fountainPenPosition.x,fountainPenPosition.y,
		     message,*nchars);
    for (i=0;i < *nchars;i++)
      Line[active_line][i]= message[i];
#endif
    for (i = *nchars;i<80;i++)
      Line[active_line][i]= ' ';
    advanceLine();
    XFlush(dpy); 
  }

void showChar(c)
     char c;
  {
    XDrawImageString(dpy,theWindow,fountainPenGC,
		     fountainPenPosition.x,fountainPenPosition.y,&c,1);
    fountainPenPosition.x += char_spacing;
    XDrawImageString(dpy,theWindow,fountainPenGC,
		     fountainPenPosition.x,fountainPenPosition.y,"_",1);
    XFlush(dpy);
  }
     

void backSpace()
  {
    
    if (input_line_idx == 0) return;
    XDrawImageString(dpy,theWindow,fountainPenGC,
		     fountainPenPosition.x,fountainPenPosition.y,"  ",2);
    fountainPenPosition.x -= char_spacing;
    if (fountainPenPosition.x < 10) fountainPenPosition.x = 10;
    XDrawImageString(dpy,theWindow,fountainPenGC,
		     fountainPenPosition.x,fountainPenPosition.y,"_",1);    
    XFlush(dpy);
    input_line[--input_line_idx]=' ';
  }
     



void TypeChar(c)
     char c;
  {
    int i;

    if (c=='\003') Quit();
    switch (c) 
	{
	case '\r':
	case '\n':
          XDrawImageString(dpy,theWindow,fountainPenGC,
			   fountainPenPosition.x,fountainPenPosition.y,"  ",2);
	  for (i=0;i<80;i++)
	    Line[active_line][i]=input_line[i];
	  advanceLine();
	  lineEntered=1;
	  return;
	
	case '\177':
	case '\010':
	  backSpace();
	  return;
	}
    
    if (c >= 32 && c <=127)  {  
      if (input_line_idx == 80) 
	XBell(dpy,0);
      else {
	showChar(c);
	input_line[input_line_idx++]=c;
      }
     }
    else XBell(dpy,0);
  }



void xmainloop (code)
     int code;
  {
    Window wind;
    static XEvent event;

    XSync(dpy,True);
    switch (code)
	{
	case MOUSE:
	    {
	      XButtonPressedEvent *but_event; 
	      for (;;) {
		XMaskEvent(dpy,ButtonPressMask|KeyPressMask,&event);
		if (event.type==ButtonPress)
		  break;
	      }
	      but_event = (XButtonPressedEvent *) &event;
	      wind = but_event->window;
	      
	      if (wind==theWindow)
		{
		      if ((but_event->button & 0xff) == Button1)
			  mouseRecord.buttonNumber=1;
		      else if ((but_event->button & 0xff) == Button2)
			  mouseRecord.buttonNumber=2;
		      else if ((but_event->button & 0xff) == Button3)
			  mouseRecord.buttonNumber=3;
		      mouseRecord.xclicked=but_event->x;
		      mouseRecord.yclicked=but_event->y;
		 }
	    }
	    break;
	    
	case LINE:
	      {
		XKeyPressedEvent *key_event;
		char *st, keybuf[10];
		int i;
	      keyIn:
		for (;;) {
		  XMaskEvent(dpy,ButtonPressMask|KeyPressMask,&event);
		  if (event.type==KeyPress)
		    break;
		}
		key_event = (XKeyPressedEvent *) &event;
		wind=key_event->window;
		if (wind==theWindow)
		    {
		      i=XLookupString(key_event,keybuf,sizeof(keybuf),NULL,NULL);
		      for (st=keybuf;i>0;i--)
			TypeChar(*st++);
		      if (!lineEntered) goto keyIn;
		    } 
	      } 
	    break;
	  default:
	    printf("unexpected code in xmainloop %d\n",code);
	    exit(-1);
	    break;
	  }
  }
    
int Quit()
  {
    XDestroyWindow(dpy,theWindow);
    XSync(dpy,0);
    exit(0);
  }
     
#if (CRAY|TITAN)
void CLOSEW()
#else
void closew_()
#endif
  {
    Quit();
  }

void error (mesg)
     char *mesg;
  {
    fprintf(stderr, mesg);
    exit(0);
  }


/*************************************
  main program for PRENEK and POSTNEK
 *************************************/
#if (CRAY|TITAN)
    void MAINC()
#else
    void mainc_()
#endif
  {
    /* Pop up XWindows */
    SetUpEnv ();

    /* call the Fortran code */
#if (CRAY | TITAN)
    FPREP();
#else
    fprep_ ();      /* run Fortran Main code */
#endif

    Quit ();  
  }

#ifdef CRAY
void SYSTEM( string )
char *string;
{
	system( string );
}

#include <sys/stat.h>
#include <errno.h>

void RENAME( from, to )
char *from, *to;
{
	rename( from, to );
}
#endif
