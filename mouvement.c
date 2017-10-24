#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#include "nrdef.h"
#include "nrutil.h"

#define N 2
#define VMIN 2
#define VMAX 255 //V est entre 2 et 2^m-1 avec m le nombre de bits des donnees ici 8 => https://hal.inria.fr/hal-01130889/document

__attribute__ ((always_inline))int inline min(int a, int b)
{
    return a<b ? a:b;
}



__attribute__ ((always_inline))int inline max(int a, int b)
{
    return a>b ? a:b;
}


uint8 ** routine_FrameDifference(uint8 **I0, uint8 **I1, uint8 **E0, long nrl,long nrh,long ncl,long nch, int seuil)
{
    //m[nrl..nrh][ncl..nch]

    uint8 **O0 = ui8matrix(nrl, nrh, ncl, nch);
    for(int i = nrl; i < nrh; i++ )
    {
        for(int j = ncl; j < nch; j++)
        {
            O0[i][j] = abs(I0[i][j] - I1[i][j]);

        }
    }
    for(int i = nrl; i < nrh; i++ )
    {
        for(int j = ncl; j < nch; j++)
        {
            if(O0[i][j] < seuil)
                E0[i][j] = 0;
            else
                E0[i][j] = 255;

        }
    }
    free_ui8matrix(O0, nrl, nrh, ncl, nch);
    return E0;
}


void routine_SigmaDelta_step0(uint8** I, uint8 **M, uint8 **V, long nrl, long nrh, long ncl, long nch)
{
    for(int i = nrl; i < nrh; i++ )
    {
        for(int j = ncl; j < nch; j++)
        {
            M[i][j] = I[i][j];
            V[i][j] = 10; //Au depart a VMIN mais il y avait beaucoup de mouvement des le debut, a 10 ça marche mieux
        }
    }
}

void routine_SigmaDelta_1step(uint8 **I0, uint8 **I1, uint8**V0, uint8 **V1, uint8**M0, uint8 **M1, uint8 **E0,  long nrl, long nrh, long ncl, long nch )
{
    uint8 **O0 = ui8matrix(nrl, nrh, ncl, nch);

    for(int i = nrl; i < nrh; i++ ) //Step1 Mt Estimation
    {
        for(int j = ncl; j < nch; j++)
        {
            if(M1[i][j] < I0[i][j])
                M0[i][j]  = M1[i][j] + 1;

            else if(M1[i][j] < I0[i][j])
                M0[i][j] = M1[i][j] - 1;

            else
                M0[i][j] = M1[i][j];


        }
    }


    for(int i = nrl; i < nrh; i++)//Step 2 difference Computation
    {
        for(int j = ncl; j < nch; j++)
        {
            O0[i][j] = abs(M0[i][j] - I0[i][j]);
        }
    }


    for(int i = nrl; i < nrh; i++)//Step 3 Update and clamping
    {
        for(int j = ncl; j < nch; j++)
        {
            if(V1[i][j] < N * O0[i][j])
                V0[i][j]  = V1[i][j] + 1;

            else if(V1[i][j] < N * O0[i][j])
                V0[i][j] = V1[i][j] - 1;

            else
                V0[i][j] = V1[i][j];

            V0[i][j] = max( min(V0[i][j], VMAX), VMIN);


        }
    }

    for(int i = nrl; i < nrh; i++)//Step 4 Et estimation
    {
        for(int j = ncl; j < nch; j++)
        {
            if(O0[i][j] < V0[i][j])
                E0[i][j] = 0;
            else
                E0[i][j] = 255;
        }
    }
    free_ui8matrix(O0, nrl, nrh, ncl, nch);

}
