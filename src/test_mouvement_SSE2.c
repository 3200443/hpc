#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vnrutil.h"
#include "nrutil.h"
#include "mouvement_SSE2.h"
#include "morpho.h"

#define NBIMAGES 199
//I0 = It et I1 = It-1 : pareil pour tout

void init_tab(vuint8 **vX1, uint8 **Itm1, int nrl, int nrh, int ncl, int nch)
{

    vuint8 T[1];
    uint8 *p = (uint8*) T;
    int cpt = 0;

    for(int i = nrl; i <=nrh; i++)
    {
        for(int j = ncl; j <= nch; j++)
        {
            for(int k = 0; k < 16; k++){
                p[k] = Itm1[i][j*16+k];
            }
            vX1[i][j] = T[0];
        }
    }
}


void test_routine_FrameDifference_SSE2(int seuil)
{
    char nomImageLoad[50];// = "car3/car_3";
    char nomImageSave[50];// = "car3Sigma/car_3"
    long nrl, nrh, ncl, nch;
    uint8 **Itm1 =  LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **It = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et = ui8matrix(nrl, nrh, ncl, nch);
    int nrow=nrh-nrl+1,ncol=nch-ncl+1;
    int n = nrow*ncol;
    //int si0, si1, sj0, sj1;
    int vi0, vi1, vj0, vj1;
    s2v(nrl, nrh, ncl, nch, 16, &vi0, &vi1, &vj0, &vj1);
    vuint8 ** vX1 = vui8matrix_s(nrl, nrh, ncl, nch);
    init_tab(vX1, Itm1, vi0, vi1, vj0, vj1);    //On fait la copie de l'image dans une matrice SIMD

    /*display_vui8matrix(vX1, vi0, vi1, vj0, vj1, "%d ", "test2");
    printf("\n\n\n\n\n\n\n\n==============================================IMAGE REELLE %d=================================\n\n\n\n\n\n", nch);
    for(int i = nrl; i <= nrh; i++ )
    {
        for(int j = ncl; j <= nch; j++)
        {
            printf("%d ",Itm1[i][j]);
        }
        printf("\n");
    }*/



    /*for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);

        routine_FrameDifference(It, Itm1, Et, nrl,nrh,ncl,nch, seuil);
        sprintf(nomImageSave, "car3Frame/car_3%03d.pgm", i);
        SavePGM_ui8matrix(Et, nrl, nrh, ncl, nch, nomImageSave);
        memcpy(Itm1[nrl], It[nrl], sizeof(uint8)*(nrow*ncol));
    }*/
    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl, nrh, ncl, nch );

}
