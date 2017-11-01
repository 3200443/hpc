#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vnrutil.h"
#include "nrutil.h"
#include "mouvement_SSE2.h"
#include "morpho.h"

#define NBIMAGES 199
//I0 = It et I1 = It-1 : pareil pour tout

void MatScal2MatSIMD(vuint8 **vX1, uint8 **Itm1, int vi0, int vi1, int vj0, int vj1)
{
    vuint8 x;
    vuint8 T[1];
    uint8 *p = (uint8*) T;
    int cpt = 0;

    for(int i = vi0; i <=vi1; i++)
    {
        for(int j = vj0; j <= vj1; j++)
        {
            for(int k = 0; k < 16; k++)
            {
                p[k] = Itm1[i][j*16+k];
            }
            vX1[i][j] = T[0];
        }
    }
}

void MatSIMD2MatScal(vuint8 **vX1, uint8 **Itm1, int vi0, int vi1, int vj0, int vj1)
{
    vuint8 T[1];
    vuint8 x;

    uint8 *p = (uint8*) T;
    int cpt = 0;

    for(int i = vi0; i <=vi1; i++)
    {
        for(int j = vj0; j <= vj1; j++)
        {
            x = _mm_load_si128(&vX1[i][j]);
            _mm_store_si128(T, x);
            for(int k = 0; k < 16; k++)
            {
                Itm1[i][j*16+k] = p[k];
            }
        }
    }
}

void test_routine_FrameDifference_SSE2(int seuil)
{
    char nomImageLoad[50];// = "car3/car_3";
    char nomImageSave[50];// = "car3Sigma/car_3"
    long nrl, nrh, ncl, nch;

    //Partie scalaire
    uint8 **Itm1 =  LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **It = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et = ui8matrix(nrl, nrh, ncl, nch);

    // Partie vecteur
    int vi0, vi1, vj0, vj1; //Indices SIMD
    s2v(nrl, nrh, ncl, nch, 16, &vi0, &vi1, &vj0, &vj1); //Recuperation des seuils SIMD
    int nrow=vi1-vi0+1,ncol=vj1-vj0+1;

    vuint8 ** vXtm1 = vui8matrix_s(nrl, nrh, ncl, nch); //Creation d'une matrice SIMD avec les indices scalaires
    vuint8 ** vXt = vui8matrix_s(nrl, nrh, ncl, nch);
    vuint8 ** vXEt = vui8matrix_s(nrl, nrh, ncl, nch);
    vuint8 seuilSIMD = init_vuint8(seuil); //Copie du seuil dans un vecteur SIMD

    MatScal2MatSIMD(vXtm1, Itm1, vi0, vi1, vj0, vj1);    //On fait la copie de l'image dans une matrice SIMD


    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);
        MatScal2MatSIMD(vXt, It,  vi0, vi1, vj0, vj1);

        routine_FrameDifference_SSE2(vXt, vXtm1, vXEt, vi0, vi1, vj0, vj1, seuilSIMD);
        MatSIMD2MatScal(vXEt, Et, vi0, vi1, vj0, vj1);    //On fait la copie d'une matrice SIMD dans une image normale
        sprintf(nomImageSave, "car3FrameSIMD/car_3%03d.pgm", i);
        SavePGM_ui8matrix(Et, nrl, nrh, ncl, nch, nomImageSave);
        memcpy(vXtm1[vi0], vXt[vi0], sizeof(vuint8)*(nrow*ncol));
    }
    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl, nrh, ncl, nch );
    free_vui8matrix(vXtm1, vi0, vi1, vj0, vj1);
    free_vui8matrix(vXt, vi0, vi1, vj0, vj1);
    free_vui8matrix(vXEt, vi0, vi1, vj0, vj1);



}
