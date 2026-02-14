
#define Bufsz 256

FILE  *Fps;                          /* postscript file pointer */
extern float Fctr, Orgx, Orgy, Curx, Cury;  /* factor and coordinates */
float Xwid = 29.7;                   /* paper size */
float Ywid = 21.0;
float Xfact, Yfact;                  /* factor to draw */
FILE  *Fp;                           /* file pointer */
char  Buf[Bufsz];                    /* buffer area */
char  Sbuf[Bufsz];                   /* string buffer area */
char  *Fnamb = "redraw.tmp";         /* work file name */
char  Fnam[256];                     /* redraw data file name */
char  Onam[256];                     /* Postscript file name */
int   Unitf = 10;                    /* mm --> cm */
int   Lwid, Ldsh;
int   Modinp;                        /* mode */
int   Ipmode;                        /* =1,2,4 figs/paper */
int   Font1, Font2;
float Gray;
float Rgb1, Rgb2, Rgb3, Rgb4;
int   Xxe, Yye, Setf, Rot;
int   Wfs, Wfm;
char  Opnmd[2];                      /* "w"/"a" */
char  Prcmd;                         /* =1st,2nd fig */
float Flw;                           /* factor for line width */
char  printknd = 0;                  /* =1:color printer, =0: gray printer */

