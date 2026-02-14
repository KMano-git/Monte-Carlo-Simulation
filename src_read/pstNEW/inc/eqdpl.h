#define MXVTX       24
#define INBUF      128
#define MAXFNM     256
#define NHEAD      640
#define NHBAS    (NHEAD + 70)
#define MAXPOLY   4000
#define NOBUF    32000
#define MAXOBUF  NOBUF + 767 - NHEAD

extern short textfont;       /* text font */
extern short textface;       /* text face */
extern short textsize;       /* text size */

extern FILE  *fpi;           /* input  file pointor */
extern FILE  *fpo;           /* output file pointor */

extern char  obuf[MAXOBUF];  /* output buffer area */
extern char  outbnm[MAXFNM]; /* output file base name */
extern char  outfnm[MAXFNM]; /* output file name */
extern short poly[MAXPOLY];  /* polygon data */
extern char  rbuf[INBUF];    /* input buffer area */
extern char  key[3];         /* keyword */
extern char  txt[INBUF];     /* string */
extern char  vtxt[MXVTX][81];
extern short vtxp[MXVTX];
extern short vtyp[MXVTX];
extern short vthg[MXVTX];

extern short qdcurrx;        /* cuurent x */
extern short qdcurry;        /* cuurent y */
extern short shei;           /* symbol height */
extern short sang;           /* symbol angle */
extern short sheio;          /* symbol height(before step) */
extern short swid;           /* simbol width */
extern float fact;           /* scaling factor */
extern short xw;             /* frame width */
extern short yw;             /* frame height */
extern short xp;             /* x-coordinate */
extern short yp;             /* y-coordinate */
extern int   rdeof;          /* EOF flag */
extern short npoly;          /* number of polygon data */
extern short nvtxt;          /* number of vertical text */
extern short xbp;            /* x coordinate of horizontal symbol */
extern short ybp;            /* y coordinate of horizontal symbol */
extern short pen;            /* pen number */
extern short pfl;            /* fill option (=1:fill) */
extern short rr;             /* radius for circle */
extern short ico;            /* color number */
extern short xf0;            /* frame lower left  x coordinate */
extern short yf0;            /* frame lower left  y coordinate */
extern short xf1;            /* frame upper right x coordinate */
extern short yf1;            /* frame upper right y coordinate */
extern short psbf;           /* buffer area for pack */
extern char  pcbf;           /* buffer area for pack */
extern int   obpnt;          /* output buffer pointer */
extern int   obgsp;          /* graph start point at output buffer */
extern int   pictsz;         /* file size */
extern short nfig;           /* PICT file number */
extern short mxvpx;          /* maximum vertical pixels */
extern short mxhpx;          /* maximum horizontal pixels */
extern float smpxs;          /* SMPXL*SMPXL */
