#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#include "nrdef.h"
#include "nrutil.h"
#define BORD 2

void erosion3x3(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch)
{
    int i,j,k,l;
    uint8 result;
    for(i=nrl; i<=nrh; i++)
    {
        for(j=ncl; j<=nch; j++)
        {
            result = 255;
            for(k = i-1 ; k < i+2 ; k++)
            {
                for(l = j-1; l<j+2 ; l++)
                {
                    result &= X[k][l];
                }
            }
            Y[i][j] = result;
        }
    }
}

void erosion3x3_bin(ulong32** X,ulong32** Y, long bi0,long bi1,long bj0,long bj1)
{
    int i,j,k,l;
    ulong32 result, left, right;
    for(i=bi0; i<=bi1; i++)
    {
        for(j=bj0; j<=bj1; j++)
        {
            result = 0;
            result = ~result;//Pour avoir tous les bits a 1
            result&=X[i-1][j];
            result&=X[i][j];
            result&=X[i+1][j];

            left = (X[i-1][j-1] & X[i][j-1] & X[i+1][j-1]);
            left = ( (result >> 1) & ~(1<<(TLONG-1)) ) | (left & 1) << (TLONG-1);

            right = (X[i-1][j+1] & X[i][j+1] & X[i+1][j+1]);
            right = (result << 1) | (right>>(TLONG-1) & 1) ;

            result &= right & left;


            Y[i][j] = result;
        }
    }
}

void dilatation3x3(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch)
{
    int i,j,k,l;
    uint8 result;
    for(i=nrl; i<=nrh; i++)
    {
        for(j=ncl; j<=nch; j++)
        {
            result = 0;
            for(k = i-1 ; k < i+2 ; k++)
            {
                for(l = j-1; l<j+2 ; l++)
                {
                    result |= X[k][l];
                }
            }
            Y[i][j] = result;
        }
    }

}
void dilatation3x3_bin(ulong32** X,ulong32** Y, long bi0,long bi1,long bj0,long bj1)
{
    int i,j,k,l;
    ulong32 result, left, right;
    for(i=bi0; i<=bi1; i++)
    {
        for(j=bj0; j<=bj1; j++)
        {
            result = 0;
            result|=X[i-1][j];
            result|=X[i][j];
            result|=X[i+1][j];

            left = (X[i-1][j-1] | X[i][j-1] | X[i+1][j-1]);
            left = ( (result >> 1) & ~(1<<(TLONG-1) )) | (left & 1) << (TLONG-1);

            right = (X[i-1][j+1] | X[i][j+1] | X[i+1][j+1]);
            right = (result << 1) | (right>>(TLONG-1) & 1) ;
            result |= right | left;



            Y[i][j] = result;
        }
    }

}

void fermeture3x3(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch)
{
    uint8 **O0 = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    dilatation3x3(X, O0, nrl, nrh, ncl, nch);
    erosion3x3(O0, Y,nrl,nrh,ncl,nch);
    free_ui8matrix(O0, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
}

void ouverture3x3(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch)
{
    uint8 **O0 = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    erosion3x3(X, O0,nrl,nrh,ncl,nch);
    dilatation3x3(O0, Y, nrl, nrh, ncl, nch);
    free_ui8matrix(O0, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
}

void ouverture3x3_bin(ulong32** X,ulong32** Y, ulong32** O0, long bi0,long bi1,long bj0,long bj1)
{
    erosion3x3_bin(X, O0,bi0,bi1,bj0,bj1);
    dilatation3x3_bin(O0, Y, bi0,bi1,bj0,bj1);
}

void fermeture3x3_bin(ulong32** X,ulong32** Y, ulong32** O0, long bi0,long bi1,long bj0,long bj1)
{
    dilatation3x3_bin(X, O0,bi0,bi1,bj0,bj1);
    erosion3x3_bin(O0, Y,bi0,bi1,bj0,bj1);
}
void erosion5x5(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch)
{
    int i,j,k,l;
    uint8 result;
    for(i=nrl; i<=nrh; i++)
    {
        for(j=ncl; j<=nch; j++)
        {
            result = 255;
            for(k = i-2 ; k < i+3 ; k++)
            {
                for(l = j-2; l<j+3 ; l++)
                {
                    result &= X[k][l];
                }
            }
            Y[i][j] = result;
        }
    }
}

void dilatation5x5(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch)
{
    int i,j,k,l;
    uint8 result;
    for(i=nrl; i<=nrh; i++)
    {
        for(j=ncl; j<=nch; j++)
        {
            result = 0;
            for(k = i-2 ; k < i+3 ; k++)
            {
                for(l = j-2; l<j+3 ; l++)
                {
                    result |= X[k][l];
                }
            }
            Y[i][j] = result;
        }
    }
}

void fermeture5x5(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch)
{
    uint8 **O0 = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    dilatation5x5(X, O0, nrl, nrh, ncl, nch);
    erosion5x5(O0, Y,nrl,nrh,ncl,nch);
    free_ui8matrix(O0, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
}


//Loop unwinding + rotation de registre + scalarisation
void erosion3x3_line(uint8** X,uint8** Y, int i,long ncl,long nch)
{
    uint8 result;
    int j = ncl;
    uint8 Xim1jm1, Xim1j, Xim1jp1;

    uint8 Xijm1, Xij, Xijp1;

    uint8 Xip1jm1, Xip1j, Xip1jp1;

    Xim1jm1 = X[i-1][j-1];
    Xijm1=X[i][j-1];
    Xip1jm1=X[i+1][j-1];

    Xim1j=X[i-1][j];
    Xij=X[i][j];
    Xip1j=X[i+1][j];


    for(j=ncl; j<=nch; j++)
    {
        result = 255;

        Xim1jp1=X[i-1][j+1];
        Xijp1=X[i][j+1];
        Xip1jp1=X[i+1][j+1];

        result &= Xim1jm1 & Xim1j & Xim1jp1 & Xijm1 & Xij & Xijp1 &
                  Xip1jm1 & Xip1j & Xip1jp1 ;//Calculs
        Y[i][j] = result;

        Xim1jm1 = Xim1j;
        Xijm1 = Xij;
        Xip1jm1 = Xip1j;

        Xim1j = Xim1jp1;
        Xij = Xijp1;
        Xip1j = Xip1jp1;
    }

}

//Loop unwinding + rotation de registre + scalarisation
void dilatation3x3_line(uint8** X,uint8** Y, int i,long ncl,long nch)
{
    uint8 result;
    int j = ncl;
    uint8 Xim1jm1, Xim1j, Xim1jp1;

    uint8 Xijm1, Xij, Xijp1;

    uint8 Xip1jm1, Xip1j, Xip1jp1;

    Xim1jm1 = X[i-1][j-1];
    Xijm1=X[i][j-1];
    Xip1jm1=X[i+1][j-1];

    Xim1j=X[i-1][j];
    Xij=X[i][j];
    Xip1j=X[i+1][j];


    for(j=ncl; j<=nch; j++)
    {

        result = 0;

        Xim1jp1=X[i-1][j+1];
        Xijp1=X[i][j+1];
        Xip1jp1=X[i+1][j+1];

        result |= Xim1jm1 | Xim1j | Xim1jp1 | Xijm1 | Xij | Xijp1 |
                  Xip1jm1 | Xip1j | Xip1jp1 ;//Calculs
        Y[i][j] = result;

        //j dans j-1
        Xim1jm1 = Xim1j;
        Xijm1 = Xij;
        Xip1jm1 = Xip1j;

        //j+1 dans j
        Xim1j = Xim1jp1;
        Xij = Xijp1;
        Xip1j = Xip1jp1;


    }

}

void ouverture3x3_pipe(uint8** X,uint8** Y, uint8** O0, long nrl,long nrh,long ncl,long nch)
{
    int i = nrl, j = 0;
    //Creation de la premiere ligne
    erosion3x3_line(X, O0, i, ncl, nch);

    for(i = nrl; i <= nrh-1; i++)
    {
        erosion3x3_line(X, O0, i+1, ncl, nch);
        dilatation3x3_line(O0, Y, i, ncl, nch);
    }
    dilatation3x3_line(O0, Y, nrh, ncl, nch);
}

void fermeture3x3_pipe(uint8** X,uint8** Y,uint8** O0, long nrl,long nrh,long ncl,long nch)
{
    int i = nrl, j;
    //Creation de la premiere ligne
    dilatation3x3_line(X, O0, i, ncl, nch);

    for(i = nrl; i <= nrh-1; i++)
    {
        dilatation3x3_line(X, O0, i+1, ncl, nch);
        erosion3x3_line(O0, Y, i, ncl, nch);

    }
    erosion3x3_line(O0, Y, nrh, ncl, nch);
}

void ouverture5x5(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch)
{
    uint8 **O0 = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    erosion5x5(X, O0,nrl,nrh,ncl,nch);
    dilatation5x5(O0, Y, nrl, nrh, ncl, nch);
    free_ui8matrix(O0, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
}
