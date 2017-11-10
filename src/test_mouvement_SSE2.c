#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vnrutil.h"
#include "nrutil.h"
#include "mouvement_SSE2.h"
#include "morpho.h"
#include "morpho_simd.h"

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

void test_unitaire_SD_SSE2()
{
    /*Test du step1*/
    printf("Test Step 1 avec(Mtm1-It) \n127-128, 255-255, 0-0, 0-1, 1-0, 255-254, 254-255,50-60,70-45, 50-50, + Cas non critiques\n");
    vuint8 Mtm1 = init_vuint8_all(127, 255, 0, 0, 1, 255, 254, 50, 70, 50, 104, 17, 195, 6, 90, 156);
    vuint8 It = init_vuint8_all(128, 255, 0, 1, 0, 254, 255, 60, 45, 50, 133, 67, 149, 198, 191, 68);
    vuint8 tmpMt;
    vuint8 un = init_vuint8(1);

    vuint8 maxSChar = init_vuint8(128);

    display_vuint8(It," %d ","ItStep1");
    printf("\n");
    display_vuint8(Mtm1," %d ","Mtm1Step1");
    printf("\n");


    vuint8 Mtm1Plus1 = _mm_add_epi8(Mtm1, un);
    vuint8 Mtm1Moins1 = _mm_sub_epi8(Mtm1, un);

    vuint8 res = _mm_cmplt_epi8(_mm_sub_epi8(Mtm1, maxSChar), _mm_sub_epi8(It, maxSChar));
    tmpMt = _mm_or_si128(_mm_and_si128(res, Mtm1Plus1), _mm_andnot_si128(res, Mtm1)); //Mtm1< It

    res = _mm_cmpgt_epi8(_mm_sub_epi8(Mtm1, maxSChar), _mm_sub_epi8(It, maxSChar));
    tmpMt = _mm_or_si128(_mm_and_si128(res, Mtm1Moins1), _mm_andnot_si128(res, tmpMt)); // //Mtm1 > It

    display_vuint8(tmpMt, " %d ","Mt");
    printf("\n");

    printf("Resultat attendu de Mt\n: 128, 255, 0, 2, 254, 255, 51, 69, 50, 105, 18, 194, 7, 91, 155\n");







    /*vuint8 pixelNoir = init_vuint8(0);
    vuint8 pixelBlanc = init_vuint8(255);
    vuint8 tmpEt;
    vuint8 un = init_vuint8(1);
    vuint8 VMAXSIMD = init_vuint8(VMAX);
    vuint8 VMINSIMD = init_vuint8(VMIN);*/



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

void test_routine_FrameDifference_SSE2M(int seuil)
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
    //
    vuint8 ** vXEt1 = vui8matrix_s(nrl, nrh, ncl, nch);
    vuint8 ** vXEt2 = vui8matrix_s(nrl, nrh, ncl, nch);
    //
    vuint8 seuilSIMD = init_vuint8(seuil); //Copie du seuil dans un vecteur SIMD

    MatScal2MatSIMD(vXtm1, Itm1, vi0, vi1, vj0, vj1);    //On fait la copie de l'image dans une matrice SIMD


    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);
        MatScal2MatSIMD(vXt, It,  vi0, vi1, vj0, vj1);

        routine_FrameDifference_SSE2(vXt, vXtm1, vXEt, vi0, vi1, vj0, vj1, seuilSIMD);
        //
        dilatation3x3_SIMD(vXEt,vXEt1,vi0,vi1,vj0,vj1);
        erosion3x3_SIMD(vXEt1,vXEt2,vi0,vi1,vj0,vj1);
        //
        MatSIMD2MatScal(vXEt2, Et, vi0, vi1, vj0, vj1);    //On fait la copie d'une matrice SIMD dans une image normale
        sprintf(nomImageSave, "car3FrameSIMD_M/car_3%03d.pgm", i);
        SavePGM_ui8matrix(Et, nrl, nrh, ncl, nch, nomImageSave);
        memcpy(vXtm1[vi0], vXt[vi0], sizeof(vuint8)*(nrow*ncol));
    }
    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl, nrh, ncl, nch );
    free_vui8matrix(vXtm1, vi0, vi1, vj0, vj1);
    free_vui8matrix(vXt, vi0, vi1, vj0, vj1);
    free_vui8matrix(vXEt, vi0, vi1, vj0, vj1);
    //
    free_vui8matrix(vXEt1, vi0, vi1, vj0, vj1);
    free_vui8matrix(vXEt2, vi0, vi1, vj0, vj1);
    //


}


void test_routine_sigmaDelta_SSE2()
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

    vuint8 ** vXMt = vui8matrix_s(nrl, nrh, ncl, nch);
    vuint8 ** vXMtm1 = vui8matrix_s(nrl, nrh, ncl, nch);

    vuint8 ** vXVt = vui8matrix_s(nrl, nrh, ncl, nch);
    vuint8 ** vXVtm1 = vui8matrix_s(nrl, nrh, ncl, nch);


    MatScal2MatSIMD(vXtm1, Itm1, vi0, vi1, vj0, vj1);    //On fait la copie de l'image dans une matrice SIMD
    routine_SigmaDelta_step0SSE2(vXtm1, vXMtm1,vXVtm1, vi0, vi1, vj0, vj1);

    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);
        MatScal2MatSIMD(vXt, It,  vi0, vi1, vj0, vj1);

        routine_SigmaDelta_1stepSSE2(vXt, vXtm1, vXVt, vXVtm1, vXMt, vXMtm1, vXEt, vi0, vi1, vj0, vj1);
        MatSIMD2MatScal(vXEt, Et, vi0, vi1, vj0, vj1);    //On fait la copie d'une matrice SIMD dans une image normale
        sprintf(nomImageSave, "car3SigmaSIMD/car_3%03d.pgm", i);
        SavePGM_ui8matrix(Et, nrl, nrh, ncl, nch, nomImageSave);
        memcpy(vXtm1[vi0], vXt[vi0], sizeof(vuint8)*(nrow*ncol));
        memcpy(vXVtm1[vi0], vXVt[vi0], sizeof(vuint8)*(nrow*ncol));
        memcpy(vXMtm1[vi0], vXMt[vi0], sizeof(vuint8)*(nrow*ncol));

    }
    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl, nrh, ncl, nch );
    free_vui8matrix(vXtm1, vi0, vi1, vj0, vj1);
    free_vui8matrix(vXt, vi0, vi1, vj0, vj1);
    free_vui8matrix(vXVtm1, vi0, vi1, vj0, vj1);
    free_vui8matrix(vXVt, vi0, vi1, vj0, vj1);
    free_vui8matrix(vXMtm1, vi0, vi1, vj0, vj1);
    free_vui8matrix(vXMt, vi0, vi1, vj0, vj1);
    free_vui8matrix(vXEt, vi0, vi1, vj0, vj1);



}
