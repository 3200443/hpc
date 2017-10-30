/* ----------------- */
/* --- vnrutil.c --- */
/* ----------------- */


/*
 * Copyright (c) 2004-2006 Lionel Lacassagne
 */

#include <stdio.h>
#include <stdlib.h>

#include <xmmintrin.h>
#include <emmintrin.h>
#include <mmintrin.h>

#include "nrdef.h"
#include "vnrdef.h"
#include "vnrutil.h"

#define LINEAR_ALLOCATION
/* -------------------------- */
void vnrerror(char error_text[])
/* -------------------------- */
/* Numerical Recipes standard error handler */
{
    fprintf(stderr,"Numerical Recipes run-time error...\n");
    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"...now exiting to system...\n");
    exit(1);
}
/* ------------------------------------------------- */
void s2v1D(int sj0, int sj1, int c, int *vj0, int *vj1)
/* ------------------------------------------------- */
{
    if(sj0%c) {
        *vj0 = (sj0-c)/c;
    } else {
        *vj0 = sj0/c;
    }
    *vj1 = sj1/c; 
}
/* ------------------------------------------------- */
void v2m1D(int vj0, int vj1, int c, int *mj0, int *mj1)
/* ------------------------------------------------- */
{
    *mj0 = c * vj0;
    *mj1 = c * (vj1+1)-1;
}
/* ------------------------------------------------------------------------------------- */
void s2v(int si0, int si1, int sj0, int sj1, int c, int *vi0, int *vi1, int *vj0, int *vj1)
/* ------------------------------------------------------------------------------------- */
{
    // conversion indices scalaires en indices vectoriels
    *vi0 = si0;
    *vi1 = si1;
    
    s2v1D(sj0, sj1, c, vj0, vj1);
    
    /*if(sj0%c) {
        *vj0 = (sj0-c)/c;
    } else {
        *vj0 = sj0/c;
    }
    *vj1 = sj1/c;*/
}
/* ------------------------------------------------------------------------------------- */
void v2m(int vi0, int vi1, int vj0, int vj1, int c, int *mi0, int *mi1, int *mj0, int *mj1)
/* ------------------------------------------------------------------------------------- */
{
    // conversion indices vectoriels en indices scalaires pour malloc alignés
    *mi0 = vi0;
    *mi1 = vi1;
    
    /**mj0 = c * vj0;
    *mj1 = c * (vj1+1)-1;*/
    
    v2m1D(vj0, vj1, c, mj0, mj1);
}
/* --------------- */
int card_vuint8(void)
/* --------------- */
{
    return 16;
}
/* --------------- */
int card_vsint8(void)
/* --------------- */
{
    return 16;
}
/* ---------------- */
int card_vuint16(void)
/* ---------------- */
{
    return 8;
}
/* ---------------- */
int card_vsint16(void)
/* ---------------- */
{
    return 8;
}
/* ---------------- */
int card_vuint32(void)
/* ---------------- */
{
    return 4;
}
/* ---------------- */
int card_vsint32(void)
/* ---------------- */
{
    return 4;
}
/* ----------------- */
int card_vfloat32(void)
/* ----------------- */
{
    return 4;
}
/* ---------------------------- */
vuint8* vui8vector(int nl, int nh)
/* ---------------------------- */
{
    vuint8 *v;

    v=(vuint8 *)_mm_malloc ((size_t) ((nh-nl+1)*sizeof(vuint8)), 16);
    if (!v) vnrerror("allocation failure in vui8vector()");
    if(!v) return NULL;
    return v-nl;
}
/* ---------------------------- */
vsint8* vsi8vector(int nl, int nh)
/* ---------------------------- */
{
    vsint8 *v;
    
    v=(vsint8 *)_mm_malloc ((size_t) ((nh-nl+1)*sizeof(vsint8)), 16);
    if (!v) vnrerror("allocation failure in vsi8vector()");
    if(!v) return NULL;
    return v-nl;
}
/* ------------------------------ */
vuint16* vui16vector(int nl, int nh)
/* ------------------------------ */
{
    vuint16 *v;
    
    v=(vuint16*)_mm_malloc ((size_t) ((nh-nl+1)*sizeof(vuint16)), 16);
    if (!v) vnrerror("allocation failure in vui16vector()");
    if(!v) return NULL;
    return v-nl;
}
/* ------------------------------ */
vsint16* vsi16vector(int nl, int nh)
/* ------------------------------ */
{
    vsint16 *v;
    
    v=(vsint16*)_mm_malloc ((size_t) ((nh-nl+1)*sizeof(vsint16)), 16);
    if (!v) vnrerror("allocation failure in vsi16vector()");
    if(!v) return NULL;
    return v-nl;
}
/* ------------------------------ */
vuint32* vui32vector(int nl, int nh)
/* ------------------------------ */
{
    vuint32 *v;
    
    v=(vuint32*)_mm_malloc ((size_t) ((nh-nl+1)*sizeof(vuint32)), 16);
    if (!v) vnrerror("allocation failure in vui32vector()");
    if(!v) return NULL;
    return v-nl;
}
/* ------------------------------ */
vsint32* vsi32vector(int nl, int nh)
/* ------------------------------ */
{
    vsint32 *v;
    
    v=(vsint32*)_mm_malloc ((size_t) ((nh-nl+1)*sizeof(vsint32)), 16);
    if (!v) vnrerror("allocation failure in vsi32vector()");
    if(!v) return NULL;
    return v-nl;
}
/* -------------------------- */
vfloat* vfvector(int nl, int nh)
/* -------------------------- */
{
    vfloat *v;
    
    v=(vfloat*)_mm_malloc ((size_t) ((nh-nl+1)*sizeof(vfloat)), 16);
    if (!v) vnrerror("allocation failure in vfvector()");
    if(!v) return NULL;
    return v-nl;
}
/* ------------------------------ */
vfloat32* vf32vector(int nl, int nh)
/* ------------------------------ */
{
    vfloat32 *v;
    
    v=(vfloat32*)_mm_malloc ((size_t) ((nh-nl+1)*sizeof(vfloat32)), 16);
    if (!v) vnrerror("allocation failure in vf32vector()");
    if(!v) return NULL;
    v -= nl;
    //return v-nl;
    return v;
}
/* ----------------------------------------- */
void free_vui8vector(vuint8 *v, int nl, int nh)
/* ----------------------------------------- */
{
    _mm_free(v+nl);
}
/* ----------------------------------------- */
void free_vsi8vector(vsint8 *v, int nl, int nh)
/* ----------------------------------------- */
{
    _mm_free(v+nl);
}
/* ------------------------------------------- */
void free_vui16vector(vuint16 *v, int nl, int nh)
/* ------------------------------------------- */
{
    _mm_free(v+nl);
}
/* ------------------------------------------- */
void free_vsi16vector(vsint16 *v, int nl, int nh)
/* ------------------------------------------- */
{
    _mm_free(v+nl);
}
/* ------------------------------------------- */
void free_vui32vector(vuint32 *v, int nl, int nh)
/* ------------------------------------------- */
{
    _mm_free(v+nl);
}
/* ------------------------------------------- */
void free_vsi32vector(vsint32 *v, int nl, int nh)
/* ------------------------------------------- */
{
    _mm_free(v+nl);
}
/* --------------------------------------- */
void free_vfvector(vfloat *v, int nl, int nh)
/* --------------------------------------- */
{
    _mm_free(v+nl);
}
/* ------------------------------------------- */
void free_vf32vector(vfloat32 *v, int nl, int nh)
/* ------------------------------------------- */
{
    _mm_free(v+nl);
}
/* ---------- */
/* --- 2D --- */
/* ---------- */
/* ------------------------------------------------- */
vuint8 **vui8matrix(int nrl, int nrh, int ncl, int nch)
/* ------------------------------------------------- */
/* allocate a vuint8 matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    vuint8 **m;

    /* allocate pointers to rows */
    m=(vuint8 **) _mm_malloc ((size_t)((nrow)*sizeof(vuint8*)), 16);
    if (!m) vnrerror("allocation failure 1 in vui8matrix()");
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl]=(vuint8 *) _mm_malloc ((size_t)((nrow*ncol)*sizeof(vuint8)), 16);
    if (!m[nrl]) vnrerror("allocation failure 2 in vui8matrix()");
    m[nrl] -= ncl;
    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
    /* return pointer to array of pointers to rows */
    return m;
}
/* ------------------------------------------------- */
vsint8 **vsi8matrix(int nrl, int nrh, int ncl, int nch)
/* ------------------------------------------------- */
/* allocate a vsint8 matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    vsint8 **m;
    
    /* allocate pointers to rows */
    m=(vsint8 **) _mm_malloc ((size_t)((nrow)*sizeof(vsint8*)), 16);
    if (!m) vnrerror("allocation failure 1 in vsi8matrix()");
    m -= nrl;
    
    /* allocate rows and set pointers to them */
    m[nrl]=(vsint8 *) _mm_malloc ((size_t)((nrow*ncol)*sizeof(vsint8)), 16);
    if (!m[nrl]) vnrerror("allocation failure 2 in vui8matrix()");
    m[nrl] -= ncl;
    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
    /* return pointer to array of pointers to rows */
    return m;
}
/* --------------------------------------------------- */
vuint16 **vui16matrix(int nrl, int nrh, int ncl, int nch)
/* --------------------------------------------------- */
/* allocate a vuint16 matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    vuint16 **m;

    /* allocate pointers to rows */
    m=(vuint16 **) _mm_malloc ((size_t)((nrow)*sizeof(vuint16*)), 16);
    if (!m) vnrerror("allocation failure 1 in vui16matrix()");
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl]=(vuint16 *) _mm_malloc ((size_t)((nrow*ncol)*sizeof(vuint16)), 16);
    if (!m[nrl]) vnrerror("allocation failure 2 in vui16matrix()");
    m[nrl] -= ncl;
    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
    /* return pointer to array of pointers to rows */
    return m;
}
/* --------------------------------------------------- */
vsint16 **vsi16matrix(int nrl, int nrh, int ncl, int nch)
/* --------------------------------------------------- */
/* allocate a vsint16 matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    vsint16 **m;
    
    /* allocate pointers to rows */
    m=(vsint16 **) _mm_malloc ((size_t)((nrow)*sizeof(vsint16*)), 16);
    if (!m) vnrerror("allocation failure 1 in vsi16matrix()");
    m -= nrl;
    
    /* allocate rows and set pointers to them */
    m[nrl]=(vsint16 *) _mm_malloc ((size_t)((nrow*ncol)*sizeof(vsint16)), 16);
    if (!m[nrl]) vnrerror("allocation failure 2 in vsi16matrix()");
    m[nrl] -= ncl;
    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
    /* return pointer to array of pointers to rows */
    return m;
}
/* --------------------------------------------------- */
vuint32 **vui32matrix(int nrl, int nrh, int ncl, int nch)
/* --------------------------------------------------- */
/* allocate a vuint32 matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    vuint32 **m;
    
    /* allocate pointers to rows */
    m=(vuint32 **) _mm_malloc ((size_t)((nrow)*sizeof(vuint32*)), 16);
    if (!m) vnrerror("allocation failure 1 in vui32matrix()");
    m -= nrl;
    
    /* allocate rows and set pointers to them */
    m[nrl]=(vuint32 *) _mm_malloc ((size_t)((nrow*ncol)*sizeof(vuint32)), 16);
    if (!m[nrl]) vnrerror("allocation failure 2 in vui32matrix()");
    m[nrl] -= ncl;
    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
    /* return pointer to array of pointers to rows */
    return m;
}
/* --------------------------------------------------- */
vsint32 **vsi32matrix(int nrl, int nrh, int ncl, int nch)
/* --------------------------------------------------- */
/* allocate a vsint32 matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    vsint32 **m;
    
    /* allocate pointers to rows */
    m=(vsint32 **) _mm_malloc ((size_t)((nrow)*sizeof(vsint32*)), 16);
    if (!m) vnrerror("allocation failure 1 in vsi32matrix()");
    m -= nrl;
    
    /* allocate rows and set pointers to them */
    m[nrl]=(vsint32 *) _mm_malloc ((size_t)((nrow*ncol)*sizeof(vsint32)), 16);
    if (!m[nrl]) vnrerror("allocation failure 2 in vsi32matrix()");
    m[nrl] -= ncl;
    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
    /* return pointer to array of pointers to rows */
    return m;
}
/* --------------------------------------------------- */
vfloat32 **vf32matrix(int nrl, int nrh, int ncl, int nch)
/* --------------------------------------------------- */
/* allocate a vfloat matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    vfloat32 **m;
    
    /* allocate pointers to rows */
    m=(vfloat32 **) _mm_malloc ((size_t)((nrow)*sizeof(vfloat32*)), 16);
    if (!m) vnrerror("allocation failure 1 in vf32matrix()");
    m -= nrl;
    
    /* allocate rows and set pointers to them */
    m[nrl]=(vfloat32 *) _mm_malloc ((size_t)((nrow*ncol)*sizeof(vfloat32)), 16);
    if (!m[nrl]) vnrerror("allocation failure 2 in vf32matrix()");
    m[nrl] -= ncl;
    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
    /* return pointer to array of pointers to rows */
    return m;
}
/* -------------------------------------------------------------- */
void free_vui8matrix(vuint8 **m, int nrl, int nrh, int ncl, int nch)
/* -------------------------------------------------------------- */
{
    _mm_free(m[nrl]+ncl);
    _mm_free(m+nrl);
}
/* -------------------------------------------------------------- */
void free_vsi8matrix(vsint8 **m, int nrl, int nrh, int ncl, int nch)
/* -------------------------------------------------------------- */
{
    _mm_free(m[nrl]+ncl);
    _mm_free(m+nrl);
}
/* ---------------------------------------------------------------- */
void free_vui16matrix(vuint16 **m, int nrl, int nrh, int ncl, int nch)
/* ---------------------------------------------------------------- */
{
    _mm_free(m[nrl]+ncl);
    _mm_free(m+nrl);
}
/* ---------------------------------------------------------------- */
void free_vsi16matrix(vsint16 **m, int nrl, int nrh, int ncl, int nch)
/* ---------------------------------------------------------------- */
{
    _mm_free(m[nrl]+ncl);
    _mm_free(m+nrl);
}
/* ---------------------------------------------------------------- */
void free_vui32matrix(vuint32 **m, int nrl, int nrh, int ncl, int nch)
/* ---------------------------------------------------------------- */
{
    _mm_free(m[nrl]+ncl);
    _mm_free(m+nrl);
}
/* --------------------------------------------------------------- */
void free_vsi32matrix(vsint32 **m, int nrl, int nrh, int ncl, int nch)
/* --------------------------------------------------------------- */
{
    _mm_free(m[nrl]+ncl);
    _mm_free(m+nrl);
}
/* ---------------------------------------------------------------- */
void free_vf32matrix(vfloat32 **m, int nrl, int nrh, int ncl, int nch)
/* ---------------------------------------------------------------- */
{
    _mm_free(m[nrl]+ncl);
    _mm_free(m+nrl);
}
/* ----------------------------------------------- */
vuint8 **vui8matrix_s(int i0, int i1, int j0, int j1)
/* ----------------------------------------------- */
{
    int c = 16;
    int vi0 = i0;
    int vi1 = i1;
    int vj0, vj1;
    s2v(i0, i1, j0, j1, c, &vi0, &vi1, &vj0, &vj1);
    return vui8matrix(vi0, vi1, vj0, vj1);
}
/* ----------------------------------------------- */
vsint8 **vsi8matrix_s(int i0, int i1, int j0, int j1)
/* ----------------------------------------------- */
{
    int c = 16;
    int vi0 = i0;
    int vi1 = i1;
    int vj0, vj1;
    s2v(i0, i1, j0, j1, c, &vi0, &vi1, &vj0, &vj1);
    return vsi8matrix(vi0, vi1, vj0, vj1);
}
/* ------------------------------------------------- */
vuint16 **vui16matrix_s(int i0, int i1, int j0, int j1)
/* ------------------------------------------------- */
{
    int c = 8;
    int vi0 = i0;
    int vi1 = i1;
    int vj0, vj1;
    s2v(i0, i1, j0, j1, c, &vi0, &vi1, &vj0, &vj1);
    return vui16matrix(vi0, vi1, vj0, vj1);
}
/* ------------------------------------------------- */
vsint16 **vsi16matrix_s(int i0, int i1, int j0, int j1)
/* ------------------------------------------------- */
{
    int c = 8;
    int vi0 = i0;
    int vi1 = i1;
    int vj0, vj1;
    s2v(i0, i1, j0, j1, c, &vi0, &vi1, &vj0, &vj1);
    return vsi16matrix(vi0, vi1, vj0, vj1);
}
/* ----------------------------------------------- */
vuint32 **vui32matrix_s(int i0, int i1, int j0, int j1)
/* ----------------------------------------------- */
{
    int c = 4;
    int vi0 = i0;
    int vi1 = i1;
    int vj0, vj1;
    s2v(i0, i1, j0, j1, c, &vi0, &vi1, &vj0, &vj1);
    return vui32matrix(vi0, vi1, vj0, vj1);
}
/* ----------------------------------------------- */
vsint32 **vsi32matrix_s(int i0, int i1, int j0, int j1)
/* ----------------------------------------------- */
{
    int c = 4;
    int vi0 = i0;
    int vi1 = i1;
    int vj0, vj1;
    s2v(i0, i1, j0, j1, c, &vi0, &vi1, &vj0, &vj1);
    return vsi32matrix(vi0, vi1, vj0, vj1);
}
/* ------------------------------------------------- */
vfloat32 **vf32matrix_s(int i0, int i1, int j0, int j1)
/* ------------------------------------------------- */
{
    int c = 4;
    int vi0 = i0;
    int vi1 = i1;
    int vj0, vj1;
    s2v(i0, i1, j0, j1, c, &vi0, &vi1, &vj0, &vj1);
    return vf32matrix(vi0, vi1, vj0, vj1);
}

/* --------------------- */
vuint8 init_vuint8(uint8 x)
/* --------------------- */
{
    vuint8 T[1];
    uint8 *p = (uint8*) T;
    
    p[0]  = p[1]  = p[2]  = p[3]  = x;
    p[4]  = p[5]  = p[6]  = p[7]  = x;
    p[8]  = p[9]  = p[10] = p[11] = x;
    p[12] = p[13] = p[14] = p[15] = x;
    
    return T[0];
}
/* --------------------- */
vsint8 init_vsint8(sint8 x)
/* --------------------- */
{
    vsint8 T[1];
    sint8 *p = (sint8*) T;
    
    p[0]  = p[1]  = p[2]  = p[3]  = x;
    p[4]  = p[5]  = p[6]  = p[7]  = x;
    p[8]  = p[9]  = p[10] = p[11] = x;
    p[12] = p[13] = p[14] = p[15] = x;
    
    return T[0];
}
/* ------------------------ */
vuint16 init_vuint16(uint16 x)
/* ------------------------ */
{
    vuint16 T[1];
    uint16 *p = (uint16*) T;
    
    p[0] = p[1] = p[2] = p[3] = x;
    p[4] = p[5] = p[6] = p[7] = x;
    
    return T[0];
}
/* ------------------------ */
vsint16 init_vsint16(sint16 x)
/* ------------------------ */
{
    vsint16 T[1];
    sint16 *p = (sint16*) T;
    
    p[0] = p[1] = p[2] = p[3] = x;
    p[4] = p[5] = p[6] = p[7] = x;
    
    return T[0];
}
/* ------------------------ */
vuint32 init_vuint32(uint32 x)
/* ------------------------ */
{
    vuint32 T[1];
    uint32 *p = (uint32*) T;
    
    p[0] = p[1] = p[2] = p[3] = x;
    
    return T[0];
}
/* ------------------------ */
vsint32 init_vsint32(sint32 x)
/* ------------------------ */
{
    vsint32 T[1];
    sint32 *p = (sint32*) T;
    
    p[0] = p[1] = p[2] = p[3] = x;
    
    return T[0];
}
/* --------------------------- */
vfloat32 init_vfloat32(float32 x)
/* --------------------------- */
{
    vfloat32 T[1];
    float32 *p = (float32*) T;
    
    p[0] = p[1] = p[2] = p[3] = x;
    
    return T[0];
}
/* ---------------------------------------- */
vuint8 init_vuint8_param(uint8 x0, uint8 step)
/* ---------------------------------------- */
{
    int i;
    vuint8 T[1];
    uint8 *p = (uint8*) T;
    
    for(i=0; i<16; i++) {
        p[i] = x0;
        x0 += step;
    }
    return T[0];
}
/* ---------------------------------------- */
vsint8 init_vsint8_param(sint8 x0, sint8 step)
/* ---------------------------------------- */
{
    int i;
    vsint8 T[1];
    sint8 *p = (sint8*) T;
    
    for(i=0; i<16; i++) {
        p[i] = x0;
        x0 += step;
    }
    return T[0];
}
/* -------------------------------------------- */
vuint16 init_vuint16_param(uint16 x0, uint16 step)
/* -------------------------------------------- */
{
    int i;
    vuint16 T[1];
    uint16 *p = (uint16*) T;
    
    for(i=0; i<8; i++) {
        p[i] = x0;
        x0 += step;
    }
    return T[0];
}
/* -------------------------------------------- */
vsint16 init_vsint16_param(sint16 x0, sint16 step)
/* -------------------------------------------- */
{
    int i;
    vsint16 T[1];
    sint16 *p = (sint16*) T;
    
    for(i=0; i<8; i++) {
        p[i] = x0;
        x0 += step;
    }
    return T[0];
}
/* -------------------------------------------- */
vuint32 init_vuint32_param(uint32 x0, uint32 step)
/* -------------------------------------------- */
{
    int i;
    vuint32 T[1];
    uint32 *p = (uint32*) T;
    
    for(i=0; i<4; i++) {
        p[i] = x0;
        x0 += step;
    }
    return T[0];
}
/* -------------------------------------------- */
vsint32 init_vsint32_param(sint32 x0, sint32 step)
/* -------------------------------------------- */
{
    int i;
    vsint32 T[1];
    sint32 *p = (sint32*) T;
    
    for(i=0; i<4; i++) {
        p[i] = x0;
        x0 += step;
    }
    return T[0];
}
/* ------------------------------------------------ */
vfloat32 init_vfloat32_param(float32 x0, float32 step)
/* ------------------------------------------------ */
{
    int i;
    vfloat32 T[1];
    float32 *p = (float32*) T;
    
    for(i=0; i<4; i++) {
        p[i] = x0;
        x0 += step;
    }
    return T[0];
}
/* -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- */
vuint8 init_vuint8_all(uint8 x0, uint8 x1, uint8 x2, uint8 x3, uint8 x4, uint8 x5, uint8 x6, uint8 x7, uint8 x8, uint8 x9, uint8 x10, uint8 x11, uint8 x12, uint8 x13, uint8 x14, uint8 x15)
/* -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- */
{
    vuint8 T[1];
    uint8 *p = (uint8*) T;
    
    p[0] = x0;
    p[1] = x1;
    p[2] = x2;
    p[3] = x3;
    p[4] = x4;
    p[5] = x5;
    p[6] = x6;
    p[7] = x7;
    p[8] = x8;
    p[9] = x9;
    p[10] = x10;
    p[11] = x11;
    p[12] = x12;
    p[13] = x13;
    p[14] = x14;
    p[15] = x15;
    
    return T[0];
}
/* -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- */
vsint8 init_vsint8_all(sint8 x0, sint8 x1, sint8 x2, sint8 x3, sint8 x4, sint8 x5, sint8 x6, sint8 x7, sint8 x8, sint8 x9, sint8 x10, sint8 x11, sint8 x12, sint8 x13, sint8 x14, sint8 x15)
/* -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- */
{
    vsint8 T[1];
    sint8 *p = (sint8*) T;
    
    p[0] = x0;
    p[1] = x1;
    p[2] = x2;
    p[3] = x3;
    p[4] = x4;
    p[5] = x5;
    p[6] = x6;
    p[7] = x7;
    p[8] = x8;
    p[9] = x9;
    p[10] = x10;
    p[11] = x11;
    p[12] = x12;
    p[13] = x13;
    p[14] = x14;
    p[15] = x15;
    
    return T[0];
}
/* ---------------------------------------------------------------------------------------------------------- */
vuint16 init_vuint16_all(uint16 x0, uint16 x1, uint16 x2, uint16 x3, uint16 x4, uint16 x5, uint16 x6, uint16 x7)
/* ---------------------------------------------------------------------------------------------------------- */
{
    vuint16 T[1];
    uint16 *p = (uint16*) T;
    
    p[0] = x0;
    p[1] = x1;
    p[2] = x2;
    p[3] = x3;
    p[4] = x4;
    p[5] = x5;
    p[6] = x6;
    p[7] = x7;
    
    return T[0];
}
/* ---------------------------------------------------------------------------------------------------------- */
vsint16 init_vsint16_all(sint16 x0, sint16 x1, sint16 x2, sint16 x3, sint16 x4, sint16 x5, sint16 x6, sint16 x7)
/* ---------------------------------------------------------------------------------------------------------- */
{
    vsint16 T[1];
    sint16 *p = (sint16*) T;
    
    p[0] = x0;
    p[1] = x1;
    p[2] = x2;
    p[3] = x3;
    p[4] = x4;
    p[5] = x5;
    p[6] = x6;
    p[7] = x7;
    
    return T[0];
}
/* -------------------------------------------------------------- */
vuint32 init_vuint32_all(uint32 x0, uint32 x1, uint32 x2, uint32 x3)
/* -------------------------------------------------------------- */
{
    vuint32 T[1];
    uint32 *p = (uint32*) T;
    
    p[0] = x0;
    p[1] = x1;
    p[2] = x2;
    p[3] = x3;
    
    return T[0];
}
/* -------------------------------------------------------------- */
vsint32 init_vsint32_all(sint32 x0, sint32 x1, sint32 x2, sint32 x3)
/* -------------------------------------------------------------- */
{
    vsint32 T[1];
    sint32 *p = (sint32*) T;
    
    p[0] = x0;
    p[1] = x1;
    p[2] = x2;
    p[3] = x3;
    
    return T[0];
}
/* -------------------------------------------------------------------- */
vfloat32 init_vfloat32_all(float32 x0, float32 x1, float32 x2, float32 x3)
/* -------------------------------------------------------------------- */
{
    vfloat32 T[1];
    float32 *p = (float32*) T;
    
    p[0] = x0;
    p[1] = x1;
    p[2] = x2;
    p[3] = x3;
    
    return T[0];
}
/* ------------------ */
/* --- 1D Altivec --- */
/* ------------------ */

/* ----------------------------------------- */
void zero_vsi8vector(vsint8 *v, int j0, int j1)
/* ----------------------------------------- */
{
    int j;
    //vsint8 z = vec_splat_s8(0); // AV
    vsint8 z = init_vsint8(0);
    
    for(j=j0; j<=j1; j++) {
         //vec_st(z, 0, &v[j]);
        _mm_store_si128(&v[j], z);
    }
}
/* ----------------------------------------- */
void zero_vui8vector(vuint8 *v, int j0, int j1)
/* ----------------------------------------- */
{
    int j;
    //vuint8 z = vec_splat_u8(0); // AV
    vuint8 z = init_vuint8(0);

    for(j=j0; j<=j1; j++) {
        //vec_st(z, 0, &v[j]);
        _mm_store_si128(&v[j], z);
    }
}
/* ------------------------------------------- */
void zero_vsi16vector(vsint16 *v, int j0, int j1)
/* ------------------------------------------- */
{
    int j;
    //vsint16 z = vec_splat_s16(0); // AV
    vsint16 z = init_vsint16(0);

    for(j=j0; j<=j1; j++) {
        //vec_st(z, 0, &v[j]);
        _mm_store_si128(&v[j], z);
    }
}
/* ------------------------------------------- */
void zero_vui16vector(vuint16 *v, int j0, int j1)
/* ------------------------------------------- */
{
    int j;
    //vuint16 z = vec_splat_u16(0); //AV
    vuint16 z = init_vuint16(0);
    
    for(j=j0; j<=j1; j++) {
        //vec_st(z, 0, &v[j]);
        _mm_store_si128(&v[j], z);
    }
}
/* ------------------------------------------- */
void zero_vsi32vector(vsint32 *v, int j0, int j1)
/* ------------------------------------------- */
{
    int j;
    //vsint32 z = vec_splat_s32(0); // AV
    vsint32 z = init_vsint32(0);
    
    for(j=j0; j<=j1; j++) {
        //vec_st(z, 0, &v[j]);
        _mm_store_si128(&v[j], z);
    }
}
/* -------------------------------------------- */
void zero_vui32vectorx(vuint32 *v, int j0, int j1)
/* -------------------------------------------- */
{
    int j;
    //vuint32 z = vec_splat_u32(0); // AV
    vuint32 z = init_vuint32(0);

    for(j=j0; j<=j1; j++) {
        //vec_st(z, 0, &v[j]);
        _mm_store_si128(&v[j], z);
    }
}
/* ------------------------------------------- */
void zero_vf32vector(vfloat32 *v, int j0, int j1)
/* ------------------------------------------- */
{
    int j;
    //vfloat32 z = (vfloat32) (0,0,0,0); // AV
    vfloat32 z = init_vfloat32(0);

    for(j=j0; j<=j1; j++) {
        //vec_st(z, 0, &v[j]);
        _mm_store_ps((float*)&v[j], z);
    }
}
/* --------------------------------------------------- */
void dup_vui8vector(vuint8 *X, int j0, int j1, vuint8 *Y)
/* --------------------------------------------------- */
{
    int j;
    vuint8 x;
    
    for(j=j0; j<=j1; j++) {
        //x = vec_ld(0, &X[j]);
        //vec_st(x, 0, &Y[j]);
        x = _mm_load_si128(&X[j]);
        _mm_store_si128(&Y[j], x);
    }
}
/* --------------------------------------------------- */
void dup_vsi8vector(vsint8 *X, int j0, int j1, vsint8 *Y)
/* --------------------------------------------------- */
{
    int j;
    vsint8 x;
    
    for(j=j0; j<=j1; j++) {
        //x = vec_ld(0, &X[j]);
        //vec_st(x, 0, &Y[j]);
         x = _mm_load_si128(&X[j]);
        _mm_store_si128(&Y[j], x);
    }
}
/* ------------------------------------------------------ */
void dup_vui16vector(vuint16 *X, int j0, int j1, vuint16 *Y)
/* ------------------------------------------------------ */
{
    int j;
    vuint16 x;
    
    for(j=j0; j<=j1; j++) {
        //x = vec_ld(0, &X[j]);
        //vec_st(x, 0, &Y[j]);
         x = _mm_load_si128(&X[j]);
        _mm_store_si128(&Y[j], x);
    }
}
/* ------------------------------------------------------ */
void dup_vsi16vector(vsint16 *X, int j0, int j1, vsint16 *Y)
/* ------------------------------------------------------ */
{
    int j;
    vsint16 x;
    
    for(j=j0; j<=j1; j++) {
        //x = vec_ld(0, &X[j]);
        //vec_st(x, 0, &Y[j]);
         x = _mm_load_si128(&X[j]);
        _mm_store_si128(&Y[j], x);
    }
}
/* ------------------------------------------------------ */
void dup_vui32vector(vuint32 *X, int j0, int j1, vuint32 *Y)
/* ------------------------------------------------------ */
{
    int j;
    vuint32 x;
    
    for(j=j0; j<=j1; j++) {
        //x = vec_ld(0, &X[j]);
        //vec_st(x, 0, &Y[j]);
        x = _mm_load_si128(&X[j]);
        _mm_store_si128(&Y[j], x);
    }
}
/* ------------------------------------------------------ */
void dup_vsi32vector(vsint32 *X, int j0, int j1, vsint32 *Y)
/* ------------------------------------------------------ */
{
    int j;
    vsint32 x;
    
    for(j=j0; j<=j1; j++) {
        //x = vec_ld(0, &X[j]);
        //vec_st(x, 0, &Y[j]);
         x = _mm_load_si128(&X[j]);
        _mm_store_si128(&Y[j], x);
    }
}
/* ------------------------------------------------------- */
void dup_vf32vector(vfloat32 *X, int j0, int j1, vfloat32 *Y)
/* ------------------------------------------------------- */
{
    int j;
    vfloat32 x;
    
    for(j=j0; j<=j1; j++) {
        //x = vec_ld(0, &X[j]);
        //vec_st(x, 0, &Y[j]);
         x = _mm_load_ps((float*) &X[j]);
        _mm_store_ps((float*) &Y[j], x);
    }
}
/* ------------------ */
/* --- 2D Altivec --- */
/* ------------------ */
/* --------------------------------------------------------------------- */
void dup_vui8matrix(vuint8 **X, int i0, int i1, int j0, int j1, vuint8 **Y)
/* --------------------------------------------------------------------- */
{
    int i, j;
    
    vuint8 x;
    
    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            
            //x = vec_ld(0, &X[i][j]);
            //vec_st(x, 0, &Y[i][j]);
            x = _mm_load_si128(&X[i][j]);
            _mm_store_si128(&Y[i][j], x);
        }
    }
}
/* --------------------------------------------------------------------- */
void dup_vsi8matrix(vsint8 **X, int i0, int i1, int j0, int j1, vsint8 **Y)
/* --------------------------------------------------------------------- */
{
    int i, j;
    
    vsint8 x;
    
    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            
            //x = vec_ld(0, &X[i][j]);
            //vec_st(x, 0, &Y[i][j]);
             x = _mm_load_si128(&X[i][j]);
            _mm_store_si128(&Y[i][j], x);
        }
    }
}
/* ------------------------------------------------------------------------ */
void dup_vui16matrix(vuint16 **X, int i0, int i1, int j0, int j1, vuint16 **Y)
/* ------------------------------------------------------------------------ */
{
    int i, j;
    
    vuint16 x;
    
    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            
            //x = vec_ld(0, &X[i][j]);
            //vec_st(x, 0, &Y[i][j]);
             x = _mm_load_si128(&X[i][j]);
            _mm_store_si128(&Y[i][j], x);
        }
    }
}
/* ------------------------------------------------------------------------ */
void dup_vsi16matrix(vsint16 **X, int i0, int i1, int j0, int j1, vsint16 **Y)
/* ------------------------------------------------------------------------ */
{
    int i, j;
    
    vsint16 x;
    
    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            
            //x = vec_ld(0, &X[i][j]);
            //vec_st(x, 0, &Y[i][j]);
            x = _mm_load_si128(&X[i][j]);
            _mm_store_si128(&Y[i][j], x);
        }
    }
}
/* ------------------------------------------------------------------------ */
void dup_vui32matrix(vuint32 **X, int i0, int i1, int j0, int j1, vuint32 **Y)
/* ------------------------------------------------------------------------ */
{
    int i, j;
    
    vuint32 x;
    
    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            
            //x = vec_ld(0, &X[i][j]);
            //vec_st(x, 0, &Y[i][j]);
             x = _mm_load_si128(&X[i][j]);
            _mm_store_si128(&Y[i][j], x);
        }
    }
}
/* ------------------------------------------------------------------------ */
void dup_vsi32matrix(vsint32 **X, int i0, int i1, int j0, int j1, vsint32 **Y)
/* ------------------------------------------------------------------------ */
{
    int i, j;
    
    vsint32 x;
    
    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            
            //x = vec_ld(0, &X[i][j]);
            //vec_st(x, 0, &Y[i][j]);
             x = _mm_load_si128(&X[i][j]);
            _mm_store_si128(&Y[i][j], x);
        }
    }
}
/* ------------------------------------------------------------------------- */
void dup_vf32matrix(vfloat32 **X, int i0, int i1, int j0, int j1, vfloat32 **Y)
/* ------------------------------------------------------------------------- */
{
    int i, j;
    
    vfloat32 x;
    
    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            
            //x = vec_ld(0, &X[i][j]);
            //vec_st(x, 0, &Y[i][j]);
             x = _mm_load_ps((float*)&X[i][j]);
            _mm_store_ps((float*) &Y[i][j], x);
        }
    }
}
/* ---------------------------------------------------------------------------- */
void init_vf32vector_param(vfloat32 *v, int j0, int j1, float32 x0, float32 xstep)
/* ---------------------------------------------------------------------------- */
{
    int j;
    //vfloat32 vx = init_vfloat32_all(x0, x0+xstep, x0+2*xstep, x0+3*xstep);
    vfloat32 vx = init_vfloat32_param(x0, xstep);
    vfloat32 vxstep = init_vfloat32(4*xstep);
    
    for(j=j0; j<=j1; j++) {
        //vec_st(vx, 0, &v[j]);
        //vx = vec_add(vx, vxstep);
        _mm_store_ps((float*) &v[j], vx);
        vx = _mm_add_ps(vx, vxstep);
    }
}
/* ---------------------------------------------------------- */
void init_vui8matrix(vuint8 **m, int i0, int i1, int j0, int j1)
/* ---------------------------------------------------------- */
{
    int i, j;
    //vuint8 x, x0 = (vuint8) (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);
    vuint8 x, x0 = init_vuint8_param(1, 1);
    //vuint8 xstep = (vuint8) (16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16);
    vuint8 xstep = init_vuint8(16);
    vuint8 ystep = init_vuint8(10);
    
    for(i=i0; i<=i1; i++) {
        x = x0;
        for(j=j0; j<=j1; j++) {
            //vec_st(x, 0, &m[i][j]);
            //x = vec_add(x, xstep);
            _mm_store_si128(&m[i][j], x);
            x = _mm_add_epi8(x, xstep);
        }
        //x0 = vec_add(x0, ystep);
        x0 = _mm_add_epi8(x0, ystep);
    }
}
/* ------------------------------------------------------------ */
void init_vui16matrix(vuint16 **m, int i0, int i1, int j0, int j1)
/* ------------------------------------------------------------ */
{
    int i, j;
    //vuint16 x, x0 = (vuint16) (1,2,3,4,5,6,7,8);
    vuint16 x, x0 = init_vuint16_param(1, 1);
    vuint16 xstep = init_vuint16(8);
    vuint16 ystep = init_vuint16(10);
    
    for(i=i0; i<=i1; i++) {
        x = x0;
        for(j=j0; j<=j1; j++) {
            //vec_st(x, 0, &m[i][j]);
            //x = vec_add(x, xstep);
            _mm_store_si128(&m[i][j], x);
            x = _mm_add_epi16(x, xstep);
        }
        //x0 = vec_add(x0, ystep);
        x0 = _mm_add_epi16(x0, ystep);
    }
}
/* ----------------------------------------------------------- */
void init_vui32marix(vuint32 **m, int i0, int i1, int j0, int j1)
/* ----------------------------------------------------------- */
{
    int i, j;
    //vuint32 x, x0 = (vuint32) (1,2,3,4);
    vuint32 x, x0 = init_vuint32_param(1, 1); // (1, 2, 3, 4);
    vuint32 xstep = init_vuint32(4);
    vuint32 ystep = init_vuint32(8);

    //ystep = vec_add(ystep, ystep);
    //ystep = vec_add(ystep, ystep); // y = 32
    
    for(i=i0; i<=i1; i++) {
        x = x0;
        for(j=j0; j<=j1; j++) {
            //a = &m[i][j];
            //vec_st(x, 0, &m[i][j]);
            //x = vec_add(x, xstep);
            _mm_store_si128(&m[i][j], x);
            x = _mm_add_epi32(x, xstep);
        }
        //x0 = vec_add(x0, ystep);
        x0 = _mm_add_epi32(x0, ystep);
    }
}
/* ----------------------------------------------------------- */
void init_vsi32marix(vsint32 **m, int i0, int i1, int j0, int j1)
/* ----------------------------------------------------------- */
{
    int i, j;
    //vuint32 x, x0 = (vuint32) (1,2,3,4);
    vsint32 x, x0 =  init_vsint32_param(1, 1); // (1, 2, 3, 4);
    vsint32 xstep = init_vsint32(4);
    vsint32 ystep = init_vsint32(8);
    
    for(i=i0; i<=i1; i++) {
        x = x0;
        for(j=j0; j<=j1; j++) {
            //vec_st(x, 0, &m[i][j]);
            //x = vec_add(x, xstep);
            _mm_store_si128(&m[i][j], x);
            x = _mm_add_epi32(x, xstep);
        }
        //x0 = vec_add(x0, ystep);
        x0 = _mm_add_epi32(x0, ystep);
    }
}
/* ----------------------------------------------------------- */
void init_vf32marix(vfloat32 **m, int i0, int i1, int j0, int j1)
/* ----------------------------------------------------------- */
{
    int i, j;
    //vfloat32 x, x0 = (vfloat32) (1,2,3,4);
    //vfloat32 xstep = (vfloat32) (4,4,4,4);
    //vfloat32 ystep = (vfloat32) (8,8,8,8);
    
    vfloat32 x, x0 = init_vfloat32_param(1, 1);
    vfloat32 xstep = init_vfloat32(4);
    vfloat32 ystep = init_vfloat32(8);
    
    for(i=i0; i<=i1; i++) {
        x = x0;
        for(j=j0; j<=j1; j++) {
            //vec_st(x, 0, &m[i][j]);
            //x = vec_add(x, xstep);
            _mm_store_ps((float*) &m[i][j], x);
            x = _mm_add_ps(x, xstep);
        }
        //x0 = vec_add(x0, ystep);
        x0 = _mm_add_ps(x0, ystep);
    }
}
/* ---------------------------------------------------------------------------------------------------- */
void init_vui8matrix_param(vuint8 **m, int i0, int i1, int j0, int j1, uint8 x0, uint8 xstep, uint8 ystep)
/* ---------------------------------------------------------------------------------------------------- */
{
    int i, j;
    vuint8 vxstep = init_vuint8(16*xstep);
    vuint8 vystep = init_vuint8(ystep);
    //vuint8 vx, vx0 = init_vuint8_all(x0, x0+xstep, x0+2*xstep, x0+3*xstep, x0+4*xstep,  x0+5*xstep,  x0+6*xstep,  x0+7*xstep, x0+8*xstep, x0+9*xstep, x0+10*xstep, x0+11*xstep, x0+12*xstep, x0+13*xstep, x0+14*xstep, x0+15*xstep);
    vuint8 vx, vx0 = init_vuint8_param(x0, xstep);
    
    for(i=i0; i<=i1; i++) {
        vx = vx0;
        for(j=j0; j<=j1; j++) {
            //vec_st(vx, 0, &m[i][j]);
            //vx = vec_add(vx, vxstep);
            _mm_store_si128(&m[i][j], vx);
            vx = _mm_add_epi8(vx, vxstep);
        }
        //vx0 = vec_add(vx0, vystep);
        vx0 = _mm_add_epi8(vx0, vystep);
    }
}
/* --------------------------------------------------------------------------------------------------------- */
void init_vui16matrix_param(vuint16 **m, int i0, int i1, int j0, int j1, uint16 x0, uint16 xstep, uint16 ystep)
/* --------------------------------------------------------------------------------------------------------- */
{
    int i, j;
    vuint16 vx, vx0 = init_vuint16_param(x0, xstep);
    vuint16 vxstep = init_vuint16(8*xstep);
    vuint16 vystep = init_vuint16(ystep);
    
    for(i=i0; i<=i1; i++) {
        vx = vx0;
        for(j=j0; j<=j1; j++) {
            //vec_st(vx, 0, &m[i][j]);
            //vx = vec_add(vx, vxstep);
            _mm_store_si128(&m[i][j], vx);
            vx = _mm_add_epi16(vx, vxstep);
        }
        //vx0 = vec_add(vx0, vystep);
        vx0 = _mm_add_epi16(vx0, vystep);
    }
}
/* --------------------------------------------------------------------------------------------------------- */
void init_vsi16matrix_param(vsint16 **m, int i0, int i1, int j0, int j1, sint16 x0, sint16 xstep, sint16 ystep)
/* --------------------------------------------------------------------------------------------------------- */
{
    int i, j;
    vsint16 vx, vx0 = init_vsint16_param(x0, xstep);
    vsint16 vxstep  = init_vsint16(8*xstep);
    vsint16 vystep  = init_vsint16(ystep);
    
    for(i=i0; i<=i1; i++) {
        vx = vx0;
        for(j=j0; j<=j1; j++) {
            //vec_st(vx, 0, &m[i][j]);
            //vx = vec_add(vx, vxstep);
            _mm_store_si128(&m[i][j], vx);
            vx = _mm_add_epi16(vx, vxstep);
        }
        //vx0 = vec_add(vx0, vystep);
        vx0 = _mm_add_epi16(vx0, vystep);
    }
}
/* ------------------------------------------------------------------------------------------------------------ */
void init_vui32matrix_param(vuint32 **m, int i0, int i1, int j0, int j1, uint32  x0, uint32  xstep, uint32  ystep)
/* ------------------------------------------------------------------------------------------------------------ */
{
    int i, j;
    vuint32 vx, vx0 = init_vuint32_param(x0, xstep);
    vuint32 vxstep = init_vuint32(4*xstep);
    vuint32 vystep = init_vuint32(ystep);
    
    for(i=i0; i<=i1; i++) {
        vx = vx0;
        for(j=j0; j<=j1; j++) {
            //vec_st(vx, 0, &m[i][j]);
            //vx = vec_add(vx, vxstep);
            _mm_store_si128(&m[i][j], vx);
            vx = _mm_add_epi32(vx, vxstep);
        }
        //vx0 = vec_add(vx0, vystep);
        vx0 = _mm_add_epi32(vx0, vystep);
    }
}
/* ------------------------------------------------------------------------------------------------------------ */
void init_vsi32matrix_param(vsint32 **m, int i0, int i1, int j0, int j1, sint32  x0, sint32  xstep, sint32  ystep)
/* ------------------------------------------------------------------------------------------------------------ */
{
    int i, j;
    vsint32 vx, vx0 = init_vsint32_param(x0, xstep);
    vsint32 vxstep  = init_vsint32(4*xstep);
    vsint32 vystep  = init_vsint32(ystep);
    
    for(i=i0; i<=i1; i++) {
        vx = vx0;
        for(j=j0; j<=j1; j++) {
            //vec_st(vx, 0, &m[i][j]);
            //vx = vec_add(vx, vxstep);
            _mm_store_si128(&m[i][j], vx);
            vx = _mm_add_epi32(vx, vxstep);
        }
        //vx0 = vec_add(vx0, vystep);
        vx0 = _mm_add_epi32(vx0, vystep);
    }
}
/* ------------------------------------------------------------------------------------------------------------ */
void init_vf32matrix_param(vfloat32 **m, int i0, int i1, int j0, int j1, float32 x0, float32 xstep, float32 ystep)
/* ------------------------------------------------------------------------------------------------------------ */
{
    int i, j;
    vfloat32 vx, vx0 = init_vfloat32_param(x0, xstep);
    vfloat32 vxstep = init_vfloat32(4*xstep);
    vfloat32 vystep = init_vfloat32(ystep);
    
    for(i=i0; i<=i1; i++) {
        vx = vx0;
        for(j=j0; j<=j1; j++) {
            //vec_st(vx, 0, &m[i][j]);
            //vx = vec_add(vx, vxstep);
            _mm_store_ps((float*)&m[i][j], vx);
            vx = _mm_add_ps(vx, vxstep);
        }
        //vx0 = vec_add(vx0, vystep);
        vx0 = _mm_add_ps(vx0, vystep);
    }
}
/* -------------------------------------------------------- */
void init_ui8matrix(uint8 **m, int i0, int i1, int j0, int j1)
/* -------------------------------------------------------- */
{
    int i, j;
    uint8 x, x0 = 1;
    
    uint8 xstep = 1;
    uint8 ystep = 16;
    
    for(i=i0; i<=i1; i++) {
        x = x0;
        for(j=j0; j<=j1; j++) {
            m[i][j] =x;
            x += xstep;
        }
        x0 += ystep;
    }
}
/* ---------------------------------------------------------- */
void init_ui32matrix(uint32 **m, int i0, int i1, int j0, int j1)
/* ---------------------------------------------------------- */
{
    int i, j;
    uint32 x, x0 = 1;
    uint32 *a;
    
    uint32 xstep = 1;
    uint32 ystep = 8;
    
    for(i=i0; i<=i1; i++) {
        x = x0;
        for(j=j0; j<=j1; j++) {
            a = &m[i][j];
            m[i][j] = x;
            x += xstep;
        }
        x0 += ystep;
    }
}
/* ---------------------------------------------------------- */
void zero_vsi8matrix(vsint8 **m, int i0, int i1, int j0, int j1)
/* ---------------------------------------------------------- */
{
    int i, j;
    vsint8 z = init_vsint8(0);
    
    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            //vec_st(z, 0, &m[i][j]);
            _mm_store_si128(&m[i][j], z);
        }
    }
}
/* ---------------------------------------------------------- */
void zero_vui8matrix(vuint8 **m, int i0, int i1, int j0, int j1)
/* ---------------------------------------------------------- */
{
    int i, j;
    vuint8 z = init_vuint8(0);
    
    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            //vec_st(z, 0, &m[i][j]);
            _mm_store_si128(&m[i][j], z);
        }
    }
}
/* ------------------------------------------------------------ */
void zero_vsi16matrix(vsint16 **m, int i0, int i1, int j0, int j1)
/* ------------------------------------------------------------ */
{
    int i, j;
    vsint16 z = init_vsint16(0);
    
    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            //vec_st(z, 0, &m[i][j]);
            _mm_store_si128(&m[i][j], z);
        }
    }
}
/* ------------------------------------------------------------ */
void zero_vui16matrix(vuint16 **m, int i0, int i1, int j0, int j1)
/* ------------------------------------------------------------ */
{
    int i, j;
    vuint16 z = init_vuint16(0);
    
    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            //vec_st(z, 0, &m[i][j]);
            _mm_store_si128(&m[i][j], z);
        }
    }
}
/* ------------------------------------------------------------ */
void zero_vsi32matrix(vsint32 **m, int i0, int i1, int j0, int j1)
/* ------------------------------------------------------------ */
{
    int i, j;
    vsint32 z = init_vsint32(0);
    
    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            //vec_st(z, 0, &m[i][j]);
            _mm_store_si128(&m[i][j], z);
        }
    }
}
/* ------------------------------------------------------------ */
void zero_vui32matrix(vuint32 **m, int i0, int i1, int j0, int j1)
/* ------------------------------------------------------------ */
{
    int i, j;
    vuint32 z = init_vuint32(0);
    
    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            //vec_st(z, 0, &m[i][j]);
            _mm_store_si128(&m[i][j], z);
        }
    }
}
/* ------------------------------------------------------------- */
void zero_vf32matrix (vfloat32 **m, int i0, int i1, int j0, int j1)
/* ------------------------------------------------------------- */
{
    int i, j;
    vfloat32 z = init_vfloat32(0);

    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            //vec_st(z, 0, &m[i][j]);
            _mm_store_ps((float*) &m[i][j], z);
        }
    }
}
/* ------------------------------------------------- */
void display_vuint8(vuint8 x, char *format, char *name)
/* ------------------------------------------------- */
{
    int i;
    vuint8 T[1];
    uint8 *p = (uint8*) T;
    
    //vec_st(x, 0, T);
    _mm_store_si128(T, x);
    
    if(name != NULL) printf("%s", name);
    
    for(i=0; i<16; i++)
        printf(format, p[i]);
}
/* ------------------------------------------------- */
void display_vsint8(vsint8 x, char *format, char *name)
/* ------------------------------------------------- */
{
    int i;
    vsint8 T[1];
    sint8 *p = (sint8*) T;
    
    //vec_st(x, 0, T);
    _mm_store_si128(T, x);
    
    if(name != NULL) printf("%s", name);
    
    for(i=0; i<16; i++)
        printf(format, p[i]);
}
/* --------------------------------------------------- */
void display_vuint16(vuint16 x, char *format, char *name)
/* --------------------------------------------------- */
{
    int i;
    vuint16 T[1];
    uint16 *p = (uint16*) T;
    
    //vec_st(x, 0, T);
    _mm_store_si128(T, x);
    
    if(name != NULL) printf("%s", name);
    
    for(i=0; i<8; i++)
        printf(format, p[i]);
}
/* --------------------------------------------------- */
void display_vsint16(vsint16 x, char *format, char *name)
/* --------------------------------------------------- */
{
    int i;
    vsint16 T[1];
    sint16 *p = (sint16*) T;
    
    //vec_st(x, 0, T);
    _mm_store_si128(T, x);
    
    if(name != NULL) printf("%s", name);
    
    for(i=0; i<8; i++)
        printf(format, p[i]);
    
    printf("\n");
}
/* --------------------------------------------------- */
void display_vuint32(vuint32 x, char *format, char *name)
/* --------------------------------------------------- */
{
    int i;
    vuint32 T[1];
    uint32 *p = (uint32*) T;
    
    //vec_st(x, 0, T);
    _mm_store_si128(T, x);
    
    if(name != NULL) printf("%s", name);
    
    for(i=0; i<4; i++)
        printf(format, p[i]);
}
/* --------------------------------------------------- */
void display_vsint32(vsint32 x, char *format, char *name)
/* --------------------------------------------------- */
{
    int i;
    vsint32 T[1];
    sint32 *p = (sint32*) T;
    
    //vec_st(x, 0, T);
    _mm_store_si128(T, x);
    
    if(name != NULL) printf("%s", name);
    
    for(i=0; i<4; i++)
        printf(format, p[i]);
}
/* ----------------------------------------------------- */
void display_vfloat32(vfloat32 x, char *format, char *name)
/* ----------------------------------------------------- */
{
    int i;
    vfloat32 T[1];
    float32 *p = (float32*) T;
    
    //vec_st(x, 0, T);
    _mm_store_ps((float*)T, x);
    
    if(name != NULL) printf("%s", name);
    
    for(i=0; i<4; i++)
        printf(format, p[i]);
}
/* ---------------------------------------------------------- */
void display_vuint8_k(vuint8 x, char *format, char *name, int k)
/* ---------------------------------------------------------- */
{
    int i;
    vuint8 T[1];
    uint8 *p = (uint8*) T;
    
    //vec_st(x, 0, T);
    _mm_store_si128(T, x);
    
    if(name != NULL) printf("%s%d", name, k);
    
    for(i=0; i<16; i++)
        printf(format, p[i]);
}
/* ---------------------------------------------------------- */
void display_vsint8_k(vsint8 x, char *format, char *name, int k)
/* ---------------------------------------------------------- */
{
    int i;
    vsint8 T[1];
    sint8 *p = (sint8*) T;
    
    //vec_st(x, 0, T);
    _mm_store_si128(T, x);
    
    if(name != NULL) printf("%s%d", name, k);
    
    for(i=0; i<16; i++)
        printf(format, p[i]);
}
/* ------------------------------------------------------------ */
void display_vuint16_k(vuint16 x, char *format, char *name, int k)
/* ------------------------------------------------------------ */
{
    int i;
    vuint16 T[1];
    uint16 *p = (uint16*) T;
    
    //vec_st(x, 0, T);
    _mm_store_si128(T, x);
    
    if(name != NULL) printf("%s%d", name, k);
    
    for(i=0; i<8; i++)
        printf(format, p[i]);
}
/* ------------------------------------------------------------ */
void display_vsint16_k(vsint16 x, char *format, char *name, int k)
/* ------------------------------------------------------------ */
{
    int i;
    vsint16 T[1];
    sint16 *p = (sint16*) T;
    
    //vec_st(x, 0, T);
    _mm_store_si128(T, x);
    
    if(name != NULL) printf("%s%d", name, k);
    
    for(i=0; i<8; i++)
        printf(format, p[i]);
    
    printf("\n");
}
/* ------------------------------------------------------------ */
void display_vuint32_k(vuint32 x, char *format, char *name, int k)
/* ------------------------------------------------------------ */
{
    int i;
    vuint32 T[1];
    uint32 *p = (uint32*) T;
    
    //vec_st(x, 0, T);
    _mm_store_si128(T, x);
    
    if(name != NULL) printf("%s%d", name, k);
    
    for(i=0; i<4; i++)
        printf(format, p[i]);
}
/* ------------------------------------------------------------ */
void display_vsint32_k(vsint32 x, char *format, char *name, int k)
/* ------------------------------------------------------------ */
{
    int i;
    vsint32 T[1];
    sint32 *p = (sint32*) T;
    
    //vec_st(x, 0, T);
    _mm_store_si128(T, x);
    
    if(name != NULL) printf("%s%d", name, k);
    
    for(i=0; i<4; i++)
        printf(format, p[i]);
}
/* -------------------------------------------------------------- */
void display_vfloat32_k(vfloat32 x, char *format, char *name, int k)
/* -------------------------------------------------------------- */
{
    int i;
    vfloat32 T[1];
    float32 *p = (float32*) T;
    
    
    //vec_st(x, 0, T);
    _mm_store_ps((float*)T, x);
    
    if(name != NULL) printf("%s%d", name, k);
    
    for(i=0; i<4; i++)
        printf(format, p[i]);
}

/* ------------------ */
/* --- 1D Altivec --- */
/* ------------------ */

/* ----------------------------------------------------------------------- */
void display_vui8vector(vuint8 *vX, int j0, int j1, char *format, char *name)
/* ----------------------------------------------------------------------- */
{
    int j;
    vuint8 x;
    
    if(name != NULL)
        printf("%s\n", name);
    
    
    for(j=j0; j<=j1; j++) {
        //x = vec_ld(0, &vX[j]);
        x = _mm_load_si128(&vX[j]);
        display_vuint8(x, format, NULL);
    }
    printf("\n");
}
/* ----------------------------------------------------------------------- */
void display_vsi8vector(vsint8 *vX, int j0, int j1, char *format, char *name)
/* ----------------------------------------------------------------------- */
{
    int j;
    vsint8 x;
    
    if(name != NULL)
        printf("%s\n", name);
    
    
    for(j=j0; j<=j1; j++) {
        //x = vec_ld(0, &vX[j]);
        x = _mm_load_si128(&vX[j]);
        display_vsint8(x, format, NULL);
    }
    printf("\n");
}
/* ------------------------------------------------------------------------- */
void display_vui16vector(vuint16 *vX, int j0, int j1, char *format, char *name)
/* ------------------------------------------------------------------------- */
{
    int j;
    vuint16 x;
    
    if(name != NULL)
        printf("%s\n", name);
    
    
    for(j=j0; j<=j1; j++) {
        //x = vec_ld(0, &vX[j]);
        x = _mm_load_si128(&vX[j]);
        display_vuint16(x, format, NULL);
    }
    printf("\n");
}
/* ------------------------------------------------------------------------- */
void display_vsi16vector(vsint16 *vX, int j0, int j1, char *format, char *name)
/* ------------------------------------------------------------------------- */
{
    int j;
    vsint16 x;
    
    if(name != NULL)
        printf("%s\n", name);
    
    
    for(j=j0; j<=j1; j++) {
        //x = vec_ld(0, &vX[j]);
        x = _mm_load_si128(&vX[j]);
        display_vsint16(x, format, NULL);
    }
    printf("\n");
}
/* ------------------------------------------------------------------------- */
void display_vui32vector(vuint32 *vX, int j0, int j1, char *format, char *name)
/* ------------------------------------------------------------------------- */
{
    int j;
    vuint32 x;
    
    if(name != NULL)
        printf("%s\n", name);
    
    
    for(j=j0; j<=j1; j++) {
        //x = vec_ld(0, &vX[j]);
        x = _mm_load_si128(&vX[j]);
        display_vuint32(x, format, NULL);
    }
    printf("\n");
}
/* ------------------------------------------------------------------------- */
void display_vsi32vector(vsint32 *vX, int j0, int j1, char *format, char *name)
/* ------------------------------------------------------------------------- */
{
    int j;
    vsint32 x;
    
    if(name != NULL)
        printf("%s\n", name);
    
    
    for(j=j0; j<=j1; j++) {
        //x = vec_ld(0, &vX[j]);
        x = _mm_load_si128(&vX[j]);
        display_vsint32(x, format, NULL);
    }
    printf("\n");
}
/* ------------------------------------------------------------------------- */
void display_vf32vector(vfloat32 *vX, int j0, int j1, char *format, char *name)
/* ------------------------------------------------------------------------- */
{
    int j;
    vfloat32 x;
    
    if(name != NULL)
        printf("%s\n", name);
    
    
    for(j=j0; j<=j1; j++) {
        //x = vec_ld(0, &vX[j]);
        x = _mm_load_ps((float*)&vX[j]);
        display_vfloat32(x, format, "\n");
    }
    printf("\n");
}
/* ------------------ */
/* --- 2D Altivec --- */
/* ------------------ */
/* --------------------------------------------------------------------------------------- */
void display_vui8matrix(vuint8 **vX, int i0, int i1, int j0, int j1, char *format, char *name)
/* --------------------------------------------------------------------------------------- */
{
    int i, j;
    vuint8 x;
    //vuint8 *vXi;
    //uint8 *sXi;
    
    if(name != NULL)
        printf("%s\n", name);
    
    /*for(i=i0; i<=i1; i++) {
        vXi = vX[i];
        sXi = (uint8*) &vX[i][j0]; sXi -= j0;  // Beware of pointer arithmetic
        for(j=j0; j<=j1; j++) {
            printf(format, sXi[j]);
        }
        printf("\n");
    }
    printf("\n");*/
    
    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            //x = vec_ld(0, &vX[i][j]);
            x = _mm_load_si128(&vX[i][j]);
            display_vuint8(x, format, NULL);
        }
        printf("\n");
    }
    printf("\n");
}
/* ---------------------------------------------------------------------------------------- */
void display_vsi8matrix(vsint8 **vX, int i0, int i1, int j0, int j1, char *format, char *name)
/* ---------------------------------------------------------------------------------------- */
{
    int i, j;
    vsint8 x;
    
    if(name != NULL)
        printf("%s\n", name);
    
    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            //x = vec_ld(0, &vX[i][j]);
            x = _mm_load_si128(&vX[i][j]);
            display_vsint8(x, format, NULL);
        }
        printf("\n");
    }
    printf("\n");
}
/* ------------------------------------------------------------------------------------------ */
void display_vui16matrix(vuint16 **vX, int i0, int i1, int j0, int j1, char *format, char *name)
/* ------------------------------------------------------------------------------------------ */
{
    int i, j;
    vuint16 x;
    
    if(name != NULL)
        printf("%s\n", name);
    
    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            //x = vec_ld(0, &vX[i][j]);
            x = _mm_load_si128(&vX[i][j]);
            display_vuint16(x, format, NULL);
        }
        printf("\n");
    }
    printf("\n");
}
/* ------------------------------------------------------------------------------------------ */
void display_vsi16matrix(vsint16 **vX, int i0, int i1, int j0, int j1, char *format, char *name)
/* ------------------------------------------------------------------------------------------ */
{
    int i, j;
    vsint16 x;
    
    if(name != NULL)
        printf("%s\n", name);
    
    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            //x = vec_ld(0, &vX[i][j]);
            x = _mm_load_si128(&vX[i][j]);
            display_vsint16(x, format, NULL);
        }
        printf("\n");
    }
    printf("\n");
}
/* ------------------------------------------------------------------------------------------ */
void display_vui32matrix(vuint32 **vX, int i0, int i1, int j0, int j1, char *format, char *name)
/* ------------------------------------------------------------------------------------------ */
{
    int i, j;
    vuint32 x;
    //vuint32 *vXi;
    //uint32 *sXi;
    
    if(name != NULL)
        printf("%s\n", name);
    
    /*for(i=i0; i<=i1; i++) {
        vXi = vX[i];
    sXi = (uint32*) &vX[i][j0]; sXi -= j0;  // Beware of pointer arithmetic
    for(j=j0; j<=j1; j++) {
        printf(format, sXi[j]);
    }
    printf("\n");
    }
printf("\n");*/
    
    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            //x = vec_ld(0, &vX[i][j]);
            x = _mm_load_si128(&vX[i][j]);
            display_vuint32(x, format, NULL);
        }
        printf("\n");
    }
    printf("\n");
}
/* ------------------------------------------------------------------------------------------ */
void display_vsi32matrix(vsint32 **vX, int i0, int i1, int j0, int j1, char *format, char *name)
/* ------------------------------------------------------------------------------------------ */
{
    int i, j;
    vsint32 x;
    
    if(name != NULL)
        printf("%s\n", name);
    
    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            //x = vec_ld(0, &vX[i][j]);
            x = _mm_load_si128(&vX[i][j]);
            display_vsint32(x, format, NULL);
        }
        printf("\n");
    }
    printf("\n");
}
/* ------------------------------------------------------------------------------------------ */
void display_vf32matrix(vfloat32 **vX, int i0, int i1, int j0, int j1, char *format, char *name)
/* ------------------------------------------------------------------------------------------ */
{
    int i, j;
    vfloat32 x;
    //vfloat32 *vXi;
    //float32 *sXi;
    
    if(name != NULL)
        printf("%s\n", name);
    
    /*for(i=i0; i<=i1; i++) {
        vXi = vX[i];
        sXi = (uint8*) &vX[i][j0]; sXi -= j0;  // Beware of pointer arithmetic
        for(j=j0; j<=j1; j++) {
            printf(format, sXi[j]);
        }
        printf("\n");
    }
    printf("\n");*/
    
    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            //x = vec_ld(0, &vX[i][j]);
            x = _mm_load_ps((float*)&vX[i][j]);
            display_vfloat32(x, format, NULL);
        }
        printf("\n");
    }
    printf("\n");
}
