#define MXVTX       24
#define INBUF      128
#define MAXFNM     256
#define NHEAD      640
#define NHBAS    (NHEAD + 70)
#define MAXPOLY   4000
#define NOBUF    32000
#define MAXOBUF  NOBUF + 767 - NHEAD

short textfont;       /* text font */
short textface;       /* text face */
short textsize;       /* text size */

FILE  *fpi;           /* input  file pointor */
FILE  *fpo;           /* output file pointor */

char  obuf[MAXOBUF];  /* output buffer area */
char  outbnm[MAXFNM]; /* output file base name */
char  outfnm[MAXFNM]; /* output file name */
short poly[MAXPOLY];  /* polygon data */
char  rbuf[INBUF];    /* input buffer area */
char  key[3];         /* keyword */
char  txt[INBUF];     /* string */
char  vtxt[MXVTX][81];
short vtxp[MXVTX];
short vtyp[MXVTX];
short vthg[MXVTX];

short qdcurrx;        /* cuurent x */
short qdcurry;        /* cuurent y */
short shei;           /* symbol height */
short sang;           /* symbol angle */
short sheio;          /* symbol height(before step) */
short swid;           /* simbol width */
float fact;           /* scaling factor */
short xw;             /* frame width */
short yw;             /* frame height */
short xp;             /* x-coordinate */
short yp;             /* y-coordinate */
int   rdeof;          /* EOF flag */
short npoly;          /* number of polygon data */
short nvtxt;          /* number of vertical text */
short xbp;            /* x coordinate of horizontal symbol */
short ybp;            /* y coordinate of horizontal symbol */
short pen;            /* pen number */
short pfl;            /* fill option (=1:fill) */
short rr;             /* radius for circle */
short ico;            /* color number */
short xf0;            /* frame lower left  x coordinate */
short yf0;            /* frame lower left  y coordinate */
short xf1;            /* frame upper right x coordinate */
short yf1;            /* frame upper right y coordinate */
short psbf;           /* buffer area for pack */
char  pcbf;           /* buffer area for pack */
int   obpnt;          /* output buffer pointer */
int   obgsp;          /* graph start point at output buffer */
int   pictsz;         /* file size */
short nfig;           /* PICT file number */
short mxvpx;          /* maximum vertical pixels */
short mxhpx;          /* maximum horizontal pixels */
float smpxs;          /* SMPXL*SMPXL */
