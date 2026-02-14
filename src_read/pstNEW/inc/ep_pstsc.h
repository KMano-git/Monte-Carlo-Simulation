
#define Bufsz 256

extern FILE  *Fps;                          /* postscript file pointer */
extern float Fctr, Orgx, Orgy, Curx, Cury;  /* factor and coordinates */
extern float Xwid;                          /* paper size */
extern float Ywid;
extern float Xfact, Yfact;                  /* factor to draw */
extern FILE  *Fp;                           /* file pointer */
extern char  Buf[];                         /* buffer area */
extern char  Sbuf[];                        /* string buffer area */
extern char  *Fnamb;                        /* work file name */
extern char  Fnam[];                        /* redraw data file name */
extern char  Onam[];                        /* Postscript file name */
extern int   Unitf;                         /* mm --> cm */
extern int   Lwid, Ldsh;
extern int   Modinp, Ipmode;
extern int   Font1, Font2;
extern float Gray;
extern float Rgb1, Rgb2, Rgb3, Rgb4;
extern int   Xxe, Yye, Setf, Rot;
extern int   Wfs, Wfm;
extern char  Opnmd[];
extern char  Prcmd;
extern float Flw;                           /* factor for line width */
extern char  printknd;                 /* =1:color printer, =0: gray printer */
