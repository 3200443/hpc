#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vnrutil.h"
#include "nrutil.h"
#include "mouvement_SSE2.h"
#include "morpho.h"
#include "morpho_simd.h"
#include "mymacro.h"

#define N 2
#define NBIMAGES 199
#define VMIN 20
#define VMAX 240
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



/*
void test_unitaire_FD_SSE2()
{
    vuint8 tmpOt;
    vuint8 pixelNoir = init_vuint8(0);
    vuint8 pixelBlanc = init_vuint8(255);
    vuint8 maxInt8 = init_vuint8(127);
    //test seuil=20
    vuint8 testMin1 = init_vuint8_all(130,       140,       130,         140,         130,         140,         255, 244, 24, 10, 163,  3,  18, 211, 15, 254);
    vuint8 testMin2 = init_vuint8_all(110,       160,       109,         159,         111,         161,   0, 244,  2, 45, 146, 89, 145, 211, 48,   2);
                                    // +seuil, -seuil, +seuil-1, -seuil-1, +seuil+1, -seuil+1, extremes, 0, lambda...
   vuint8 testOt = init_vuint8_all(  20,        20,        21,          19,          19,          21,         255,   0, 22, 35,  17, 86, 127,   0, 33, 252);
   //testcmplt pour seuil = 20;
   vuint8 testcmplt = init_vuint8_all(0,          0,         0,        255,         255,           0,           0, 255,  0,  0, 255,  0,   0, 255,  0,   0);
   vuint8 testSeuil = init_vuint8(20);

   tmpOt = _mm_min_epu8(_mm_sub_epi8(testMin1,testMin2), _mm_sub_epi8(testMin2, testMin1) ); //min(a-b,b-a) donne la valeur absolue car on peut pas avoir de valeurs negatives
   display_vuint8(tmpOt, " %d ", "Min");
   printf("\n");
   display_vuint8(testOt, " %d ", "Ce qu'on doit avoir");
   printf("\n");

   vuint8 res = _mm_cmplt_epi8(_mm_sub_epi8(testOt, maxInt8), _mm_sub_epi8(testSeuil, maxInt8)); //Met 255 si inferieur au seuil et 0 sinon
   display_vuint8(testSeuil, " %d ", "Seuil");
   printf("\n");
   display_vuint8(res," %d ", "Res");
   printf("\n");
   display_vuint8(testcmplt, " %d ", "Ce qu'on doit avoir");
   printf("\n");

   res = _mm_andnot_si128(res, pixelBlanc);
   display_vuint8(res," %d ", "ResNot");
   printf("\n");

}
*/

void test_unitaire_FD_SSE2()
{
    vuint8 tmpOt;
    vuint8 pixelNoir = init_vuint8(0);
    vuint8 pixelBlanc = init_vuint8(255);
    vuint8 maxInt8 = init_vuint8(128);

    /*
        //test seuil=20
        vuint8 test1 = init_vuint8_all(    130, 140, 130, 140, 130, 140, 255, 244, 224,  10, 163,   3,  18, 211,  15, 254);
        vuint8 test2 = init_vuint8_all(    110, 160, 109, 159, 111, 161,   0, 244,  96,  45, 146,  89, 145, 211,  48,   2);
                                        // +seuil, -seuil, +seuil-1, -seuil-1, +seuil+1, -seuil+1, extremes, 0, lambda...
       vuint8 testOt = init_vuint8_all(    20,  20,  21,  19,  19,  21, 255,   0, 128,  35,  17,  86, 127,   0,  33, 252);
       //testcmplt pour seuil = 20;
       vuint8 testcmplt = init_vuint8_all(255, 255, 255,   0,   0, 255, 255,   0, 255, 255,   0, 255, 255,   0, 255, 255);
       vuint8 testSeuil = init_vuint8(20);
    */

    //test seuil=10
    vuint8 test1 = init_vuint8_all(    130, 140, 130, 140, 130, 140, 255, 244, 224,  10, 163,   3,  18, 211,  15, 254);
    vuint8 test2 = init_vuint8_all(    120, 150, 119, 149, 121, 151,   0, 244,  96,  45, 146,  89, 145, 211,  48,   2);
    vuint8 testOt = init_vuint8_all(    10,  10,  11,   9,   9,  11, 255,   0, 128,  35,  17,  86, 127,   0,  33, 252);
    vuint8 testcmplt = init_vuint8_all(  0,   0,   0, 255, 255,   0,   0, 255,   0,   0,   0,   0,   0, 255,   0,   0);
    vuint8 testSeuil = init_vuint8(10);


    vuint8 minI = _mm_min_epu8(test1, test2);
    vuint8 maxI = _mm_max_epu8(test1, test2);

    tmpOt = _mm_sub_epi8(maxI, minI);
    display_vuint8(tmpOt, " %d ", "Ot    ");
    printf("\n");
    display_vuint8(testOt, " %d ", "Ot ref");
    printf("\n");

    vuint8 res = _mm_cmplt_epi8(_mm_sub_epi8(testOt, maxInt8), _mm_sub_epi8(testSeuil, maxInt8)); //Met 255 si inferieur au seuil et 0 sinon
    display_vuint8(res," %d ", "Res    ");
    printf("\n");
    display_vuint8(testcmplt, " %d ", "Res ref");
    printf("\n");

    res = _mm_andnot_si128(res, pixelBlanc);
    display_vuint8(res," %d ", "ResNot ");
    printf("\n");
}

void test_unitaire_SD_SSE2()
{
    printf("\n\n\n");
    /*Test du step1*/
    //printf("Test Step 1 avec(Mtm1-It) \n127-128, 255-255, 0-0, 0-1, 1-0, 255-254, 254-255,50-60,70-45, 50-50, + Cas non critiques\n");
    vuint8 Mtm1 = init_vuint8_all(127, 255, 0, 0, 1, 255, 254, 50, 70, 50, 104, 17, 195,   6,  90, 156);
    vuint8 It = init_vuint8_all(  128, 255, 0, 1, 0, 254, 255, 60, 45, 50, 133, 67, 149, 198, 191, 68);
    vuint8 tmpMt;
    vuint8 un = init_vuint8(1);

    vuint8 maxSChar = init_vuint8(128);

    printf("=========== Step 1 =========\n");
    display_vuint8(It," %d ","ItStep1\n");
    printf("\n");
    display_vuint8(Mtm1," %d ","Mtm1Step1\n");
    printf("\n");


    vuint8 Mtm1Plus1 = _mm_add_epi8(Mtm1, un);
    vuint8 Mtm1Moins1 = _mm_sub_epi8(Mtm1, un);

    vuint8 res = _mm_cmplt_epi8(_mm_sub_epi8(Mtm1, maxSChar), _mm_sub_epi8(It, maxSChar));
    tmpMt = _mm_or_si128(_mm_and_si128(res, Mtm1Plus1), _mm_andnot_si128(res, Mtm1)); //Mtm1< It

    res = _mm_cmpgt_epi8(_mm_sub_epi8(Mtm1, maxSChar), _mm_sub_epi8(It, maxSChar));
    tmpMt = _mm_or_si128(_mm_and_si128(res, Mtm1Moins1), _mm_andnot_si128(res, tmpMt)); // //Mtm1 > It

    display_vuint8(tmpMt, " %d ","Mt\n");
    printf("\n");

    printf("Resultat attendu de Mt:\n 128, 255, 0, 1, 0, 254, 255, 51, 69, 50, 105, 18, 194, 7, 91, 155\n");

    printf("\n\n ========== Step 2 =========\n\n");
    It = init_vuint8_all(    130, 140, 130, 140, 130, 140, 255, 244, 224,  10, 163,   3,  18, 211,  15, 254);
    tmpMt = init_vuint8_all(    120, 150, 119, 149, 121, 151,   0, 244,  96,  45, 146,  89, 145, 211,  48,   2);
    vuint8 max = _mm_max_epu8(It,tmpMt);
    vuint8 min = _mm_min_epu8(It, tmpMt);
    vuint8 tmpOt = _mm_sub_epi8(max, min); //Le max - min donne la valeur absolue
    display_vuint8(It," %d ","ItStep2\n");
    printf("\n");
    display_vuint8(tmpMt," %d ","MtStep2\n");
    printf("\n");
    display_vuint8(tmpOt, " %d ","OtStep2\n");
    printf("\n");
    printf( "Resultat attendu de Ot:\n 10, 10, 11, 9, 9, 11, 255, 0, 128, 35, 17, 86, 127, 0, 33, 252\n");


    printf("\n\n ========== Step 3 =========\n\n");
    tmpOt = init_vuint8_all(0, 1, 100, 127, 128, 255, 2, 244,  96,  45, 146,  89, 145, 211,  48,   9);
    vuint8 NfoisOt = init_vuint8(0);
    for(int k = 0; k < N; k++)
    {
        NfoisOt = _mm_adds_epu8(NfoisOt, tmpOt);
    }
    display_vuint8(tmpOt, " %d ","OtStep3\n");
    printf("\n");
    display_vuint8(NfoisOt, " %d ","2*OtStep3\n");
    printf("\n");

    printf("Resultat attendu de NfoisOt:\n 0, 2, 200, 254, 255, 255, 4, 255, 192, 90, 255, 178, 255, 255, 96, 18\n");

    printf("\n\n ========== Step 3 Clamping =========\n\n");

    vuint8 VMAXSIMD = init_vuint8(VMAX);
    vuint8 VMINSIMD = init_vuint8(VMIN);
    vuint8 tmpVt = init_vuint8_all(0, 1, 5, 19, 20, 21, 22, 239,  240,  241, 255,  127, 128, 129,  130,   126);

    vuint8 tmpResVt = _mm_max_epu8(_mm_min_epu8(tmpVt, VMAXSIMD), VMINSIMD);

    display_vuint8(tmpVt, " %d ","VtStep3Avant\n");
    printf("\n");

    display_vuint8(tmpResVt, " %d ","VtStep3Apres\n");
    printf("\n");
    //VMIN = 20 et VMAX = 240 A 240 pour verifier que le max est correctement fait

    printf("Resultat attendu de Vt:\n 20, 20, 20, 20, 20, 21, 22, 239, 240, 240, 240, 127, 128, 129, 130, 126\n");


    printf("\n\n ========== Step 4 =========\n\n");

    vuint8 pixelBlanc = init_vuint8(255);
    tmpOt = init_vuint8_all(0, 127, 128, 128, 127, 10, 49, 255,    0,  230, 25,  147, 145, 255,  254,   126);
    tmpVt = init_vuint8_all(0, 128, 127, 128, 127, 50, 50,   0,  255,  241, 135,  127, 144, 254,  255,   125);
    res = _mm_cmplt_epi8(_mm_sub_epi8(tmpOt, maxSChar), _mm_sub_epi8(tmpVt, maxSChar)); //Met 255 si inferieur au seuil et 0 sinon
    res = _mm_andnot_si128(res, pixelBlanc);

    display_vuint8(tmpOt," %d ", " OtStep4\n");
    printf("\n");
    display_vuint8(tmpVt, " %d ", "VtStep4\n");
    printf("\n");

    display_vuint8(res," %d ", "Et\n");
    printf("\n");

    printf("Resultat attendu de Et:\n 255, 0, 255, 255, 255, 0, 0, 255, 0, 0, 0, 255, 255, 255, 0, 255\n");






    /*vuint8 pixelNoir = init_vuint8(0);
    vuint8 pixelBlanc = init_vuint8(255);
    vuint8 tmpEt;
    vuint8 un = init_vuint8(1);
    vuint8 VMAXSIMD = init_vuint8(VMAX);
    vuint8 VMINSIMD = init_vuint8(VMIN);*/



}

void test_routine_FrameDifference_SSE2(int seuil)
{

    /////////////// Pour le cycle par point////////////
    double cycles;

    char *format = "%6.2f \n";
    double cycleTotal = 0;
    int iter, niter = 2;
    int run, nrun = 5;
    double t0, t1, dt, tmin, t;
    ///////////////////////////////////////////////


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

    vuint8 ** vItm1 = vui8matrix_s(nrl, nrh, ncl, nch); //Creation d'une matrice SIMD avec les indices scalaires
    vuint8 ** vIt = vui8matrix_s(nrl, nrh, ncl, nch); //vIt=(uint8**)It devrait fonctionner
    vuint8 ** vEt = vui8matrix_s(nrl, nrh, ncl, nch);
    vuint8 seuilSIMD = init_vuint8(seuil); //Copie du seuil dans un vecteur SIMD

    MatScal2MatSIMD(vItm1, Itm1, vi0, vi1, vj0, vj1);    //On fait la copie de l'image dans une matrice SIMD


    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);
        MatScal2MatSIMD(vIt, It,  vi0, vi1, vj0, vj1);

        //routine_FrameDifference_SSE2(vIt, vItm1, vEt, vi0, vi1, vj0, vj1, seuilSIMD);
        CHRONO(routine_FrameDifference_SSE2(vIt, vItm1, vEt, vi0, vi1, vj0, vj1, seuilSIMD), cycles);
        cycleTotal+=cycles;

        MatSIMD2MatScal(vEt, Et, vi0, vi1, vj0, vj1);    //On fait la copie d'une matrice SIMD dans une image normale
        sprintf(nomImageSave, "car3FrameSIMD/car_3%03d.pgm", i);
        SavePGM_ui8matrix(Et, nrl, nrh, ncl, nch, nomImageSave);
        dup_vui8matrix(vIt, vi0, vi1, vj0, vj1, vItm1);
    }


    cycleTotal/=NBIMAGES;
    cycleTotal/=((nch+1)*(nrh+1));
    BENCH(printf("Cycles FD_SSE2 = "));
    BENCH(printf(format, cycleTotal));


    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl, nrh, ncl, nch );
    free_vui8matrix(vItm1, vi0, vi1, vj0, vj1);
    free_vui8matrix(vIt, vi0, vi1, vj0, vj1);
    free_vui8matrix(vEt, vi0, vi1, vj0, vj1);



}


void test_routine_FrameDifference_SSE2M(int seuil)
{
    /////////////// Pour le cycle par point////////////
    double cycles;

    char *format = "%6.2f \n";
    double cycleTotal = 0;
    int iter, niter = 2;
    int run, nrun = 5;
    double t0, t1, dt, tmin, t;
    ///////////////////////////////////////////////


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

    vuint8 ** vItm1 = vui8matrix_s(nrl, nrh, ncl, nch); //Creation d'une matrice SIMD avec les indices scalaires
    vuint8 ** vIt = vui8matrix_s(nrl, nrh, ncl, nch);
    vuint8 ** vEt = vui8matrix_s(nrl-1, nrh+1, ncl-1, nch+1);
    //
    vuint8 ** vEt1 = vui8matrix_s(nrl-1, nrh+1, ncl-1, nch+1);
    //vuint8 ** vEt1 = vui8matrix_s(nrl, nrh, ncl, nch);
    vuint8 ** vEt2 = vui8matrix_s(nrl-1, nrh+1, ncl-1, nch+1);
    //
    vuint8 seuilSIMD = init_vuint8(seuil); //Copie du seuil dans un vecteur SIMD

    MatScal2MatSIMD(vItm1, Itm1, vi0, vi1, vj0, vj1);    //On fait la copie de l'image dans une matrice SIMD


    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);
        MatScal2MatSIMD(vIt, It,  vi0, vi1, vj0, vj1);

        routine_FrameDifference_SSE2(vIt, vItm1, vEt, vi0, vi1, vj0, vj1, seuilSIMD);
        //
        //dilatation3x3_SIMD_B(vEt,vEt1,vi0,vi1,vj0,vj1);
        //erosion3x3_SIMD_B(vEt1,vEt2,vi0,vi1,vj0,vj1);
        //fermeture3x3_SIMD_B(vEt,vEt1,vi0,vi1,vj0,vj1);
        CHRONO(fermeture3x3_SIMD(vEt,vEt1,vi0,vi1,vj0,vj1), cycles);
        cycleTotal+=cycles;
        //
        MatSIMD2MatScal(vEt1, Et, vi0, vi1, vj0, vj1);    //On fait la copie d'une matrice SIMD dans une image normale
        sprintf(nomImageSave, "car3FrameSIMD_M/car_3%03d.pgm", i);
        SavePGM_ui8matrix(Et, nrl, nrh, ncl, nch, nomImageSave);
        memcpy(vItm1[vi0], vIt[vi0], sizeof(vuint8)*(nrow*ncol));
    }

    cycleTotal/=NBIMAGES;
    cycleTotal/=((nch+1)*(nrh+1));
    BENCH(printf("Cycles FermetureSSE2 = "));
    BENCH(printf(format, cycleTotal));


    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl, nrh, ncl, nch );
    free_vui8matrix(vItm1, vi0, vi1, vj0, vj1);
    free_vui8matrix(vIt, vi0, vi1, vj0, vj1);
    free_vui8matrix(vEt, vi0-1, vi1+1, vj0-1, vj1+1);
    //
    //
    free_vui8matrix(vEt1, vi0-1, vi1+1, vj0-1, vj1+1);
    //free_vui8matrix(vEt1, vi0, vi1, vj0, vj1);
    free_vui8matrix(vEt2, vi0-1, vi1+1, vj0-1, vj1+1);
    //


}


void test_routine_sigmaDelta_SSE2()
{


    /////////////// Pour le cycle par point////////////
    double cycles;

    char *format = "%6.2f \n";
    double cycleTotal = 0;
    int iter, niter = 2;
    int run, nrun = 5;
    double t0, t1, dt, tmin, t;
    ///////////////////////////////////////////////



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

    vuint8 ** vItm1 = vui8matrix_s(nrl, nrh, ncl, nch); //Creation d'une matrice SIMD avec les indices scalaires
    vuint8 ** vIt = vui8matrix_s(nrl, nrh, ncl, nch);
    vuint8 ** vEt = vui8matrix_s(nrl, nrh, ncl, nch);

    vuint8 ** vMt = vui8matrix_s(nrl, nrh, ncl, nch);
    vuint8 ** vMtm1 = vui8matrix_s(nrl, nrh, ncl, nch);

    vuint8 ** vVt = vui8matrix_s(nrl, nrh, ncl, nch);
    vuint8 ** vVtm1 = vui8matrix_s(nrl, nrh, ncl, nch);


    MatScal2MatSIMD(vItm1, Itm1, vi0, vi1, vj0, vj1);    //On fait la copie de l'image dans une matrice SIMD
    routine_SigmaDelta_step0SSE2(vItm1, vMtm1,vVtm1, vi0, vi1, vj0, vj1);

    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);
        MatScal2MatSIMD(vIt, It,  vi0, vi1, vj0, vj1);

        //routine_SigmaDelta_1stepSSE2(vIt, vItm1, vVt, vVtm1, vMt, vMtm1, vEt, vi0, vi1, vj0, vj1);
        CHRONO(routine_SigmaDelta_1stepSSE2(vIt, vItm1, vVt, vVtm1, vMt, vMtm1, vEt, vi0, vi1, vj0, vj1), cycles);
        cycleTotal+=cycles;


        MatSIMD2MatScal(vEt, Et, vi0, vi1, vj0, vj1);    //On fait la copie d'une matrice SIMD dans une image normale
        sprintf(nomImageSave, "car3SigmaSIMD/car_3%03d.pgm", i);
        SavePGM_ui8matrix(Et, nrl, nrh, ncl, nch, nomImageSave);
        dup_vui8matrix(vIt, vi0, vi1, vj0, vj1, vItm1);
        dup_vui8matrix(vVt, vi0, vi1, vj0, vj1, vVtm1);
        dup_vui8matrix(vMt, vi0, vi1, vj0, vj1, vMtm1);



    }

    cycleTotal/=NBIMAGES;
    cycleTotal/=((nch+1)*(nrh+1));
    BENCH(printf("Cycles SD_SSE2 = "));
    BENCH(printf(format, cycleTotal));



    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl, nrh, ncl, nch );
    free_vui8matrix(vItm1, vi0, vi1, vj0, vj1);
    free_vui8matrix(vIt, vi0, vi1, vj0, vj1);
    free_vui8matrix(vVtm1, vi0, vi1, vj0, vj1);
    free_vui8matrix(vVt, vi0, vi1, vj0, vj1);
    free_vui8matrix(vMtm1, vi0, vi1, vj0, vj1);
    free_vui8matrix(vMt, vi0, vi1, vj0, vj1);
    free_vui8matrix(vEt, vi0, vi1, vj0, vj1);



}
