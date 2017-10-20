/* -------------- */
/* --- nrutil --- */
/* -------------- */

/*
 * Copyright (c) 2000 - 2017, Lionel Lacassagne
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <ctype.h> 
#include <string.h>
#include <math.h> /* fabs */

#include "nrdef.h"
#include "nrutil.h"

//NR_END est maintenant defini dans nrutil.h

//#define NR_END 1
//#define FREE_ARG char*

long nr_end = NR_END;

/* ------------------------- */
void nrerror(char error_text[])
/* ------------------------- */
/* Numerical Recipes standard error handler */
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}
/* ------------------------------------------------------------------------------------- */
void generate_filename_k_ndigit(char *filename, int k, int ndigit, char *complete_filename)
/* ------------------------------------------------------------------------------------- */
{
  char *format = "%s%d.txt";

  switch(ndigit) {
    case 0 : format = "%s%d.txt";   break;
    case 1 : format = "%s%01d.txt"; break;
    case 2 : format = "%s%02d.txt"; break;
    case 3 : format = "%s%03d.txt"; break;
    case 4 : format = "%s%04d.txt"; break;
    case 5 : format = "%s%05d.txt"; break;
    case 6 : format = "%s%06d.txt"; break;
    case 7 : format = "%s%07d.txt"; break;
    case 8 : format = "%s%08d.txt"; break;
    case 9 : format = "%s%09d.txt"; break;
  }
  sprintf(complete_filename, format, filename, k);
}
/* ---------------------------------------------------------------------------------------------------------------- */
void generate_filename_k_ndigit_extension(char *filename, int k, int ndigit, char *extension, char *complete_filename)
/* ---------------------------------------------------------------------------------------------------------------- */
{
  char *format = "%s%d.%s";

  switch(ndigit) {
    case 0 : format = "%s%d.%s";   break;
    case 1 : format = "%s%01d.%s"; break;
    case 2 : format = "%s%02d.%s"; break;
    case 3 : format = "%s%03d.%s"; break;
    case 4 : format = "%s%04d.%s"; break;
    case 5 : format = "%s%05d.%s"; break;
    case 6 : format = "%s%06d.%s"; break;
    case 7 : format = "%s%07d.%s"; break;
    case 8 : format = "%s%08d.%s"; break;
    case 9 : format = "%s%09d.%s"; break;
  }
  sprintf(complete_filename, format, filename, k, extension);
}
/* -------------------------------------------------------------------------- */
void generate_path_filename(char *path, char *filename, char *complete_filename)
/* -------------------------------------------------------------------------- */
{
  sprintf(complete_filename, "%s%s", path, filename);
}
/* ----------------------------------------------------------------------------------------------------- */
void generate_path_filename_extension(char *path, char *filename, char *extension, char *complete_filename)
/* ----------------------------------------------------------------------------------------------------- */
{
  sprintf(complete_filename, "%s%s.%s", path, filename, extension);
}
/* --------------------------------------------------------------------------------------------------------------------------------- */
void generate_path_filename_k_ndigit_extension(char *path, char *filename, int k, int ndigit, char *extension, char *complete_filename)
/* --------------------------------------------------------------------------------------------------------------------------------- */
{
  char *format = "%s%s%d.%s";

  switch(ndigit) {
    case 0 : format = "%s%s%d.%s";   break;
    case 1 : format = "%s%s%01d.%s"; break;
    case 2 : format = "%s%s%02d.%s"; break;
    case 3 : format = "%s%s%03d.%s"; break;
    case 4 : format = "%s%s%04d.%s"; break;
    case 5 : format = "%s%s%05d.%s"; break;
    case 6 : format = "%s%s%06d.%s"; break;
    case 7 : format = "%s%s%07d.%s"; break;
    case 8 : format = "%s%s%08d.%s"; break;
    case 9 : format = "%s%s%09d.%s"; break;
  }
  sprintf(complete_filename, format, path, filename, k, extension);
}
/* ------------------------------------------------------------------------------------------------------------------------------------------ */
void generate_path_filename_k_ndigit_l_extension(char *path, char *filename, int k, int ndigit, int l, char *extension, char *complete_filename)
/* ------------------------------------------------------------------------------------------------------------------------------------------ */
{
  char *format = "%s%s%d_%d.%s";

  switch(ndigit) {
    case 0 : format = "%s%s%d_%d.%s";   break;
    case 1 : format = "%s%s%01d_%01d.%s"; break;
    case 2 : format = "%s%s%02d_%02d.%s"; break;
    case 3 : format = "%s%s%03d_%03d.%s"; break;
    case 4 : format = "%s%s%04d_%04d.%s"; break;
    case 5 : format = "%s%s%05d_%05d.%s"; break;
    case 6 : format = "%s%s%06d_%06d.%s"; break;
    case 7 : format = "%s%s%07d_%07d.%s"; break;
    case 8 : format = "%s%s%08d_%08d.%s"; break;
    case 9 : format = "%s%s%09d_%09d.%s"; break;
  }
  sprintf(complete_filename, format, path, filename, k, l, extension);
}
/* ------------------------------ */
sint32* si32vector(long nl, long nh)
/* ------------------------------ */
{
    sint32 *v;
    
    v=(sint32 *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(sint32)));
    if (!v) nrerror("allocation failure in si32vector()");
    if(!v) return NULL;
    return v-nl+NR_END;
}
/* ------------------------------ */
float32* f32vector(long nl, long nh)
/* ------------------------------ */
{
    float32 *v;
    
    v=(float32 *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float32)));
    if (!v) nrerror("allocation failure in f32vector()");
    if(!v) return NULL;
    return v-nl+NR_END;
}
/* ------------------------------ */
float64* f64vector(long nl, long nh)
/* ------------------------------ */
{
  float64 *v;

  v=(float64 *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float64)));
  if (!v) nrerror("allocation failure in f64vector()");
  if(!v) return NULL;
  return v-nl+NR_END;
}
/* --------------------------------------------------- */
uint8** ui8matrix(long nrl, long nrh, long ncl, long nch)
/* --------------------------------------------------- */
/* allocate an uint8 matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  uint8 **m;

  /* allocate pointers to rows */
  m=(uint8 **) malloc((size_t)((nrow+NR_END)*sizeof(uint8*)));
  if (!m) nrerror("allocation failure 1 in ui8matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(uint8 *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(uint8)));
  if (!m[nrl]) nrerror("allocation failure 2 in ui8matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}
/* --------------------------------------------------- */
sint8** si8matrix(long nrl, long nrh, long ncl, long nch)
/* --------------------------------------------------- */
/* allocate an sint8 matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  sint8 **m;

  /* allocate pointers to rows */
  m=(sint8**) malloc((size_t)((nrow+NR_END)*sizeof(sint8*)));
  if (!m) nrerror("allocation failure 1 in si8matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(sint8*) malloc((size_t)((nrow*ncol+NR_END)*sizeof(uint8)));
  if (!m[nrl]) nrerror("allocation failure 2 in si8matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}
/* --------------------------------------------------- */
uint16** ui16matrix(long nrl, long nrh, long ncl, long nch)
/* --------------------------------------------------- */
/* allocate an uint16 matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  uint16 **m;

  /* allocate pointers to rows */
  m=(uint16 **) malloc((size_t)((nrow+NR_END)*sizeof(uint16*)));
  if (!m) nrerror("allocation failure 1 in ui16matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(uint16 *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(uint16)));
  if (!m[nrl]) nrerror("allocation failure 2 in ui16matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}
/* --------------------------------------------------- */
sint16** si16matrix(long nrl, long nrh, long ncl, long nch)
/* --------------------------------------------------- */
/* allocate an sint16 matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  sint16 **m;

  /* allocate pointers to rows */
  m=(sint16 **) malloc((size_t)((nrow+NR_END)*sizeof(sint16*)));
  if (!m) nrerror("allocation failure 1 in si16matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(sint16 *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(uint16)));
  if (!m[nrl]) nrerror("allocation failure 2 in si16matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}
/* --------------------------------------------------- */
uint32** ui32matrix(long nrl, long nrh, long ncl, long nch)
/* --------------------------------------------------- */
/* allocate an uint32 matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  uint32 **m;

  /* allocate pointers to rows */
  m=(uint32 **) malloc((size_t)((nrow+NR_END)*sizeof(uint32*)));
  if (!m) nrerror("allocation failure 1 in ui16matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(uint32 *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(uint32)));
  if (!m[nrl]) nrerror("allocation failure 2 in ui16matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}
/* --------------------------------------------------- */
sint32** si32matrix(long nrl, long nrh, long ncl, long nch)
/* --------------------------------------------------- */
/* allocate an sint32 matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  sint32 **m;

  /* allocate pointers to rows */
  m=(sint32 **) malloc((size_t)((nrow+NR_END)*sizeof(sint32*)));
  if (!m) nrerror("allocation failure 1 in ui16matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(sint32 *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(sint32)));
  if (!m[nrl]) nrerror("allocation failure 2 in ui16matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}
/* ----------------------------------------------------- */
float32** f32matrix(long nrl, long nrh, long ncl, long nch)
/* ----------------------------------------------------- */
/* allocate an float32 matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  float32 **m;

  /* allocate pointers to rows */
  m=(float32 **) malloc((size_t)((nrow+NR_END)*sizeof(float32*)));
  if (!m) nrerror("allocation failure 1 in f32matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(float32 *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float32)));
  if (!m[nrl]) nrerror("allocation failure 2 in f32matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}
/* ----------------------------------------------------- */
float64** f64matrix(long nrl, long nrh, long ncl, long nch)
/* ----------------------------------------------------- */
/* allocate an float32 matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    float64 **m;
    
    /* allocate pointers to rows */
    m=(float64 **) malloc((size_t)((nrow+NR_END)*sizeof(float64*)));
    if (!m) nrerror("allocation failure 1 in f64matrix()");
    m += NR_END;
    m -= nrl;
    
    /* allocate rows and set pointers to them */
    m[nrl]=(float64 *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float64)));
    if (!m[nrl]) nrerror("allocation failure 2 in f32matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;
    
    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
    
    /* return pointer to array of pointers to rows */
    return m;
}
/* --------------------------------------------------- */
rgb8** rgb8matrix(long nrl, long nrh, long ncl, long nch)
/* -------------------------------------------------- */
/* allocate an rgb8 matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  rgb8 **m;

  /* allocate pointers to rows */
  m=(rgb8**) malloc((size_t)((nrow+NR_END)*sizeof(rgb8*)));
  if (!m) nrerror("allocation failure 1 in rgb8matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(rgb8*) malloc((size_t)((nrow*ncol+NR_END)*sizeof(rgb8)));
  if (!m[nrl]) nrerror("allocation failure 2 in rgb8matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}
/* ------------------------------------------- */
void free_si32vector(sint32 *v, long nl, long nh)
/* ------------------------------------------- */
/* free a sint32 vector allocated with si32vector() */
{
    free((FREE_ARG) (v+nl-NR_END));
}
/* ------------------------------------------- */
void free_f32vector(float32 *v, long nl, long nh)
/* ------------------------------------------- */
/* free a f32 vector allocated with f32vector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}
/* ------------------------------------------- */
void free_f64vector(float64 *v, long nl, long nh)
/* ------------------------------------------- */
/* free a double vector allocated with f64vector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}
/* ---------------------------------------------------------------- */
void free_ui8matrix(uint8 **m, long nrl, long nrh, long ncl, long nch)
/* ---------------------------------------------------------------- */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}
/* ---------------------------------------------------------------- */
void free_si8matrix(sint8 **m, long nrl, long nrh, long ncl, long nch)
/* ---------------------------------------------------------------- */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}
/* ------------------------------------------------------------------ */
void free_ui16matrix(uint16 **m, long nrl, long nrh, long ncl, long nch)
/* ------------------------------------------------------------------ */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}
/* ------------------------------------------------------------------ */
void free_si16matrix(sint16 **m, long nrl, long nrh, long ncl, long nch)
/* ------------------------------------------------------------------ */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}
/* ------------------------------------------------------------------ */
void free_ui32matrix(uint32 **m, long nrl, long nrh, long ncl, long nch)
/* ------------------------------------------------------------------ */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}
/* ------------------------------------------------------------------ */
void free_si32matrix(sint32 **m, long nrl, long nrh, long ncl, long nch)
/* ------------------------------------------------------------------ */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}
/* ------------------------------------------------------------------ */
void free_f32matrix(float32 **m, long nrl, long nrh, long ncl, long nch)
/* ------------------------------------------------------------------ */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}
/* ------------------------------------------------------------------ */
void free_f64matrix(float64 **m, long nrl, long nrh, long ncl, long nch)
/* ------------------------------------------------------------------ */
{
    free((FREE_ARG) (m[nrl]+ncl-NR_END));
    free((FREE_ARG) (m+nrl-NR_END));
}
/* ---------------------------------------------------------------- */
void free_rgb8matrix(rgb8 **m, long nrl, long nrh, long ncl, long nch)
/* ---------------------------------------------------------------- */
/* free an uchar matrix allocated by rgb8matrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}
/* ------------ */
/* --- zero --- */
/* ------------ */

/* ----------------------------------------- */
void zero_si32vector(sint32* v, int j0, int j1)
/* ----------------------------------------- */
{
    int j;
    for(j=j0; j<=j1; j++) {
        v[j] = 0;
    }
}
/* ----------------------------------------- */
void zero_f32vector(float32* v, int j0, int j1)
/* ----------------------------------------- */
{
    int j;
    for(j=j0; j<=j1; j++) {
        v[j] = 0;
    }
}
/* ----------------------------------------- */
void zero_f64vector(float64* v, int j0, int j1)
/* ----------------------------------------- */
{
    int j;
    for(j=j0; j<=j1; j++) {
        v[j] = 0;
    }
}
/* ---------------------------------------------------------- */
void zero_si32matrix(sint32** m, int i0, int i1, int j0, int j1)
/* ---------------------------------------------------------- */
{
    int i;
    
    for(i=i0; i<=i1; i++) {
        zero_si32vector(m[i], j0, j1);
    }
}
/* ---------------------------------------------------------- */
void zero_f32matrix(float32** m, int i0, int i1, int j0, int j1)
/* ---------------------------------------------------------- */
{
    int i;
    
    for(i=i0; i<=i1; i++) {
        zero_f32vector(m[i], j0, j1);
    }
}
/* ---------------------------------------------------------- */
void zero_f64matrix(float64** m, int i0, int i1, int j0, int j1)
/* ---------------------------------------------------------- */
{
    int i;
    
    for(i=i0; i<=i1; i++) {
        zero_f64vector(m[i], j0, j1);
    }
}

/* ----------------- */
/* --- set_param --- */
/* ----------------- */

/* ---------------------------------------------------------------------- */
void set_si32vector_param(sint32 *v, int j0, int j1, sint32 x, sint32 xstep)
/* ---------------------------------------------------------------------- */
{
    int j;   
    for(j=j0; j<=j1; j++) {
        v[j] = x;
        x += xstep;
    }
}
/* ------------------------------------------------------------------------ */
void set_f32vector_param(float32 *v, int j0, int j1, float32 x, float32 xstep)
/* ------------------------------------------------------------------------ */
{
    int j;   
    for(j=j0; j<=j1; j++) {
        v[j] = x;
        x += xstep;
    }
}
/* ------------------------------------------------------------------------ */
void set_f64vector_param(float64 *v, int j0, int j1, float64 x, float64 xstep)
/* ------------------------------------------------------------------------ */
{
    int j;   
    for(j=j0; j<=j1; j++) {
        v[j] = x;
        x += xstep;
    }
}
/* ----------------------------------------------------------------------------------------------------- */
void set_si32matrix_param(sint32 **m, int i0, int i1, int j0, int j1, sint32 x, sint32 xstep, sint32 ystep)
/* ----------------------------------------------------------------------------------------------------- */
{
    int i;
    for(i=i0; i<=i1; i++) {
        set_si32vector_param(m[i], j0, j1, x, xstep);
        x += ystep;
    }
}
/* -------------------------------------------------------------------------------------------------------- */
void set_f32matrix_param(float32 **m, int i0, int i1, int j0, int j1, float32 x, float32 xstep, float32 ystep)
/* -------------------------------------------------------------------------------------------------------- */
{
    int i;
    for(i=i0; i<=i1; i++) {
        set_f32vector_param(m[i], j0, j1, x, xstep);
        x += ystep;
    }
}
/* -------------------------------------------------------------------------------------------------------- */
void set_f64matrix_param(float64 **m, int i0, int i1, int j0, int j1, float64 x, float64 xstep, float64 ystep)
/* -------------------------------------------------------------------------------------------------------- */
{
    int i;
    for(i=i0; i<=i1; i++) {
        set_f64vector_param(m[i], j0, j1, x, xstep);
        x += ystep;
    }
}
/* ------------ */
/* --- copy --- */
/* ------------ */

/* ----------------------------------------------------------------------------- */
void copy_ui8matrix_ui8matrix(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
/* ----------------------------------------------------------------------------- */
{
    int i, j;

    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            Y[i][j] = X[i][j];
        }
    }
}
/* ------------------------------------------------------------------------------- */
void copy_ui8matrix_ui16matrix(uint8 **X, int i0, int i1, int j0, int j1, uint16 **Y)
/* ------------------------------------------------------------------------------- */
{
    int i, j;

    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            Y[i][j] = X[i][j];
        }
    }
}
/* ------------------------------------------------------------------------------- */
void copy_ui16matrix_ui8matrix(uint16 **X, int i0, int i1, int j0, int j1, uint8 **Y)
/* ------------------------------------------------------------------------------- */
{
    int i, j;

    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            Y[i][j] = (uint8) X[i][j];
        }
    }
}
/* ------------------------------------------------------------------------------- */
void copy_ui8matrix_f32matrix(uint8 **X, int i0, int i1, int j0, int j1, float32 **Y)
/* ------------------------------------------------------------------------------- */
{
    int i, j;

    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            Y[i][j] = X[i][j];
        }
    }
}
/* ------------------------------------------------------------------------------- */
void copy_f32matrix_ui8matrix(float32 **X, int i0, int i1, int j0, int j1, uint8 **Y)
/* ------------------------------------------------------------------------------- */
{
    int i, j;

    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            Y[i][j] = (uint8) X[i][j];
        }
    }
}

/* -------------------------------------------------------------------- */
void display_si8vector(sint8 *v,long nl,long nh, char *format, char *name)
/* -------------------------------------------------------------------- */
{
  long i;

  if(name != NULL) printf("%s", name);

  for(i=nl; i<=nh; i++) {
    printf(format, v[i]);
  }
  putchar('\n');
}
/* -------------------------------------------------------------------- */
void display_ui8vector(uint8 *v,long nl,long nh, char *format, char *name)
/* -------------------------------------------------------------------- */
{
  long i;

  if(name != NULL) printf("%s", name);

  for(i=nl; i<=nh; i++) {
    printf(format, v[i]);
  }
  putchar('\n');
}
/* ---------------------------------------------------------------------- */
void display_si16vector(sint16 *v,long nl,long nh, char *format, char *name)
/* ---------------------------------------------------------------------- */
{
  long i;

  if(name != NULL) printf("%s", name);

  for(i=nl; i<=nh; i++) {
    printf(format, v[i]);
  }
  putchar('\n');
}
/* ---------------------------------------------------------------------- */
void display_ui16vector(uint16 *v,long nl,long nh, char *format, char *name)
/* ---------------------------------------------------------------------- */
{
  long i;

  if(name != NULL) printf("%s", name);

  for(i=nl; i<=nh; i++) {
    printf(format, v[i]);
  }
  putchar('\n');
}
/* ---------------------------------------------------------------------- */
void display_si32vector(sint32 *v,long nl,long nh, char *format, char *name)
/* ---------------------------------------------------------------------- */
{
  long i;

  if(name != NULL) printf("%s", name);

  for(i=nl; i<=nh; i++) {
    printf(format, v[i]);
  }
  putchar('\n');
}
/* ---------------------------------------------------------------------- */
void display_ui32vector(uint32 *v,long nl,long nh, char *format, char *name)
/* ---------------------------------------------------------------------- */
{
  long i;

  if(name != NULL) printf("%s", name);

  for(i=nl; i<=nh; i++) {
    printf(format, v[i]);
  }
  putchar('\n');
}
/* ---------------------------------------------------------------------- */
void display_f32vector(float32 *v,long nl,long nh, char *format, char *name)
/* ---------------------------------------------------------------------- */
{
  long i;

  if(name != NULL) printf("%s", name);

  for(i=nl; i<=nh; i++) {
    printf(format, v[i]);
  }
  putchar('\n');
}
/* ---------------------------------------------------------------------- */
void display_f64vector(float64 *v,long nl,long nh, char *format, char *name)
/* ---------------------------------------------------------------------- */
{
    long i;
    
    if(name != NULL) printf("%s", name);
    
    for(i=nl; i<=nh; i++) {
        printf(format, v[i]);
    }
    putchar('\n');
}
/* ------------------------------------------------------------------------------------------ */
void display_si8matrix(sint8 **m,long nrl,long nrh,long ncl, long nch, char *format, char *name)
/* ------------------------------------------------------------------------------------------ */
{
  long i,j;

  if(name != NULL) puts(name);

  for(i=nrl; i<=nrh; i++) {
    for(j=ncl; j<=nch; j++) {
      printf(format, m[i][j]);
    }
    putchar('\n');
  }
}
/* ------------------------------------------------------------------------------------------ */
void display_ui8matrix(uint8 **m,long nrl,long nrh,long ncl, long nch, char *format, char *name)
/* ------------------------------------------------------------------------------------------ */
{
  long i,j;

  if(name != NULL) puts(name);

  for(i=nrl; i<=nrh; i++) {
    for(j=ncl; j<=nch; j++) {
      printf(format, m[i][j]);
    }
    putchar('\n');
  }
}
/* -------------------------------------------------------------------------------------------- */
void display_si16matrix(sint16 **m,long nrl,long nrh,long ncl, long nch, char *format, char *name)
/* -------------------------------------------------------------------------------------------- */
{
  long i,j;

  if(name != NULL) puts(name);

  for(i=nrl; i<=nrh; i++) {
    for(j=ncl; j<=nch; j++) {
      printf(format, m[i][j]);
    }
    putchar('\n');
  }
}
/* -------------------------------------------------------------------------------------------- */
void display_ui16matrix(uint16 **m,long nrl,long nrh,long ncl, long nch, char *format, char *name)
/* -------------------------------------------------------------------------------------------- */
{
  long i,j;

  if(name != NULL) puts(name);

  for(i=nrl; i<=nrh; i++) {
    for(j=ncl; j<=nch; j++) {
      printf(format, m[i][j]);
    }
    putchar('\n');
  }
}
/* -------------------------------------------------------------------------------------------- */
void display_si32matrix(sint32 **m,long nrl,long nrh,long ncl, long nch, char *format, char *name)
/* -------------------------------------------------------------------------------------------- */
{
  long i,j;

  if(name != NULL) puts(name);

  for(i=nrl; i<=nrh; i++) {
    for(j=ncl; j<=nch; j++) {
      printf(format, m[i][j]);
    }
    putchar('\n');
  }
}
/* -------------------------------------------------------------------------------------------- */
void display_ui32matrix(uint32 **m,long nrl,long nrh,long ncl, long nch, char *format, char *name)
/* -------------------------------------------------------------------------------------------- */
{
  long i,j;

  if(name != NULL) puts(name);

  for(i=nrl; i<=nrh; i++) {
    for(j=ncl; j<=nch; j++) {
      printf(format, m[i][j]);
    }
    putchar('\n');
  }
}
/* -------------------------------------------------------------------------------------------- */
void display_f32matrix(float32 **m,long nrl,long nrh,long ncl, long nch, char *format, char *name)
/* -------------------------------------------------------------------------------------------- */
{
  long i,j;

  if(name != NULL) puts(name);

  for(i=nrl; i<=nrh; i++) {
    for(j=ncl; j<=nch; j++) {
      printf(format, m[i][j]);
    }
    putchar('\n');
  }
}
/* -------------------------------------------------------------------------------------------- */
void display_f64matrix(float64 **m,long nrl,long nrh,long ncl, long nch, char *format, char *name)
/* -------------------------------------------------------------------------------------------- */
{
    long i,j;
    
    if(name != NULL) puts(name);
    
    for(i=nrl; i<=nrh; i++) {
        for(j=ncl; j<=nch; j++) {
            printf(format, m[i][j]);
        }
        putchar('\n');
    }
}
/* ------------------------ */
/* -- PGM IO for bmatrix -- */
/* ------------------------ */

char *readitem   (FILE *file, char *buffer);
void  ReadPGMrow (FILE *file, long width, uint8  *line);
void  WritePGMrow(uint8 *line, long width, FILE  *file);

/* --------------------------------- */
char *readitem(FILE *file,char *buffer)
/* --------------------------------- */
/* lecture d'un mot */
{
  char *aux;
  int k;

  k=0;
  aux=buffer;
  while (!feof(file))
    {
      *aux=fgetc(file);
      switch(k)
        {
        case 0:
          if (*aux=='#') k=1;
          if (isalnum(*aux)) k=2,aux++;
          break;
        case 1:
          if (*aux==0xA) k=0;
          break;
        case 2:
          if (!isalnum(*aux))
            {
              *aux=0;
              return buffer;
            }
          aux++;
          break;
        }
    }
  *aux=0;
  return buffer;
}
/* ---------------------------------------------- */
void ReadPGMrow(FILE *file, long width, uint8  *line)
/* ---------------------------------------------- */
{
    /* Le fichier est ouvert (en lecture) et ne sera pas ferme a la fin */
     fread(&(line[0]), sizeof(uint8), width, file);
}
/* ----------------------------------------------- */
void WritePGMrow(uint8 *line, long width, FILE  *file)
/* ----------------------------------------------- */
{
/* Le fichier est deja ouvert et ne sera pas ferme a la fin */

   fwrite(&(line[0]), sizeof(uint8), width, file);
}
/* ------------------------------------------------------------------------------ */
uint8** LoadPGM_ui8matrix(char *filename, long *nrl, long *nrh, long *ncl, long *nch)
/* ------------------------------------------------------------------------------ */
{
  /* cette version ne lit plus que le type P5 */

  long height, width, gris;
  uint8 **m;
  FILE *file;
  /*int   format;/**/

  char *buffer;
  /*char  c;/**/
  int i;
  
  buffer = (char*) calloc(80, sizeof(char));
  /* ouverture du fichier */
  file = fopen(filename,"rb");
  if (file==NULL)
    nrerror("ouverture du fichier impossible\n");
    //nrerror("ouverture du fichier %s impossible\n", filename);

  /* lecture de l'entete du fichier pgm */
  readitem(file, buffer);
  /*fscanf(fichier, "%s", buffer);*/
  if(strcmp(buffer, "P5") != 0)
    nrerror("entete du fichier %s invalide\n");
    //nrerror("entete du fichier %s invalide\n", filename);

  width  = atoi(readitem(file, buffer));
  height = atoi(readitem(file, buffer));
  gris   = atoi(readitem(file, buffer));

  *nrl = 0;
  *nrh = height - 1;
  *ncl = 0;
  *nch = width - 1;
  m = ui8matrix(*nrl, *nrh, *ncl, *nch);
  
  for(i=0; i<height; i++) {
    ReadPGMrow(file, width, m[i]);
  }

  fclose(file);
  free(buffer);

  return m;
}
/* ------------------------------------------------------------------------------- */
void MLoadPGM_ui8matrix(char *filename, int nrl, int nrh, int ncl, int nch, uint8 **m)
/* ------------------------------------------------------------------------------- */
{
    /* cette version ne lit plus que le type P5 */
    
    int height, width, gris;
    FILE *file;
    
    char *buffer;
    int i;
    
    buffer = (char*) calloc(80, sizeof(char));
    /* ouverture du fichier */
    file = fopen(filename,"rb");
    if (file==NULL)
        nrerror("ouverture du fichier impossible\n");
    //nrerror("ouverture du fichier %s impossible\n", filename);
    
    /* lecture de l'entete du fichier pgm */
    readitem(file, buffer);
    /*fscanf(fichier, "%s", buffer);*/
    if(strcmp(buffer, "P5") != 0)
        nrerror("entete du fichier %s invalide\n");
    //nrerror("entete du fichier %s invalide\n", filename);
    
    width  = atoi(readitem(file, buffer));
    height = atoi(readitem(file, buffer));
    gris   = atoi(readitem(file, buffer));
    
    for(i=0; i<height; i++) {
        ReadPGMrow(file, width, m[i]);
    }
    
    fclose(file);
    free(buffer);
}
/* ----------------------------------------------------------------------------------- */
void SavePGM_ui8matrix(uint8 **m, long nrl, long nrh, long ncl, long nch, char *filename)
/* ----------------------------------------------------------------------------------- */
{
  long nrow = nrh-nrl+1;
  long ncol = nch-ncl+1;

  char buffer[80];
  
  FILE *file;
  int  i;

  file = fopen(filename, "wb");
  if (file == NULL)
    //nrerror("ouverture du fichier %s impossible dans SavePGM_bmatrix\n", filename);
    nrerror("ouverture du fichier %s impossible dans SavePGM_ui8matrix\n");

  /* enregistrement de l'image au format rpgm */

  sprintf(buffer,"P5\n%d %d\n255\n",ncol, nrow);
  fwrite(buffer,strlen(buffer),1,file);
  for(i=nrl; i<=nrh; i++)
    WritePGMrow(m[i], ncol, file);

  /* fermeture du fichier */
  fclose(file);
}
/* --------------------------- */
/* -- PNM IO for rgb8matrix -- */
/* --------------------------- */

/* ----------------------------------------------- */
void ReadPNMrow(FILE  *file, long width, uint8  *line)
/* ----------------------------------------------- */
{
    /* Le fichier est ouvert (en lecture) et ne sera pas ferme a la fin */
     fread(&(line[0]), sizeof(uint8), 3*sizeof(uint8)*width, file);
}
/* ------------------------------------------------ */
void WritePNMrow(uint8  *line, long width, FILE  *file)
/* ------------------------------------------------ */
{
/* Le fichier est deja ouvert et ne sera pas ferme a la fin */

   fwrite(&(line[0]), sizeof(uint8), 3*sizeof(uint8)*width, file);
}
/* ------------------------------------------------------------------------------- */
rgb8** LoadPPM_rgb8matrix(char *filename, long *nrl, long *nrh, long *ncl, long *nch)
/* ------------------------------------------------------------------------------- */
{
  /* cette version ne lit plus que le type P6 */

  long height, width, gris;
  rgb8 **m;
  FILE *file;
  /*int   format;/**/

  char *buffer;
  /*char  c;/**/
  int i;
  
  buffer = (char*) calloc(80, sizeof(char));
  /* ouverture du fichier */
  file = fopen(filename,"rb");
  if (file==NULL)
    nrerror("ouverture du fichier impossible\n");
    //nrerror("ouverture du fichier %s impossible\n", filename);

  /* lecture de l'entete du fichier pgm */
  readitem(file, buffer);
  /*fscanf(fichier, "%s", buffer);*/
  if(strcmp(buffer, "P6") != 0)
    nrerror("entete du fichier %s invalide\n");
    //nrerror("entete du fichier %s invalide\n", filename);

  width  = atoi(readitem(file, buffer));
  height = atoi(readitem(file, buffer));
  gris   = atoi(readitem(file, buffer));

  *nrl = 0;
  *nrh = height - 1;
  *ncl = 0;
  *nch = width - 1;
  m = rgb8matrix(*nrl, *nrh, *ncl, *nch);
  
  for(i=0; i<height; i++) {
    ReadPNMrow(file, width, (uint8*)m[i]);
  }

  fclose(file);
  free(buffer);

  return m;
}
/* ----------------------------------------------------------------------------------- */
void SavePPM_rgb8matrix(rgb8 **m, long nrl, long nrh, long ncl, long nch, char *filename)
/* ----------------------------------------------------------------------------------- */
{
  long nrow = nrh-nrl+1;
  long ncol = nch-ncl+1;

  char buffer[80];
  
  FILE *file;
  int  i;

  file = fopen(filename, "wb");
  if (file == NULL)
    //nrerror("ouverture du fichier %s impossible dans SavePGM_bmatrix\n", filename);
    nrerror("ouverture du fichier %s impossible dans SavePPM_bmatrix\n");

  /* enregistrement de l'image au format rpgm */

  sprintf(buffer,"P6\n%d %d\n255\n",ncol, nrow);
  fwrite(buffer,strlen(buffer),1,file);
  for(i=nrl; i<=nrh; i++)
    WritePNMrow((uint8*)m[i], ncol, file);

  /* fermeture du fichier */
  fclose(file);
}
