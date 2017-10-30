/* --------------- */
/* --- nrdef.h --- */
/* --------------- */


#ifndef __NR_DEF_H__
#define __NR_DEF_H__

typedef unsigned char byte;

typedef          char  int8;
typedef unsigned char uint8;
typedef   signed char sint8;

typedef          short  int16;
typedef unsigned short uint16;
typedef   signed short sint16;

typedef          int  int32;
typedef unsigned int uint32;
typedef   signed int sint32;

typedef float         float32;
typedef double        float64;
typedef struct { byte r; byte g; byte b;} rgb8;

#endif // __NR_DEF_H__
