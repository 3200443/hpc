#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vnrutil.h"
#include "nrutil.h"
#include "morpho_simd.h"
#include "test_mouvement_SSE2.h"
#include "test_morpho.h"
#include "mouvement.h"

#include "mymacro.h"
#include "morpho.h"
#define NBIMAGES 199
#define BORD 2


void test_morpho3x3simd()
{
    ///////////////////////////////////////////////


    char nomImageLoad[50];// = "car3/car_3";
    char nomImageSave[50];// = "car3Sigma/car_3"
    long nrl, nrh, ncl, nch;

    //Partie scalaire
    uint8 **Itm1 =  LoadPGM_ui8matrix("testmorpho/image_depart.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **Et = ui8matrix(nrl, nrh, ncl, nch);
    // Partie vecteur
    vuint8 ** vEt = vui8matrix_s(nrl, nrh, ncl, nch);
    int vi0, vi1, vj0, vj1; //Indices SIMD
    s2v(nrl, nrh, ncl, nch, 16, &vi0, &vi1, &vj0, &vj1); //Recuperation des seuils SIMD
    int nrow=vi1-vi0+1,ncol=vj1-vj0+1;

    vuint8 ** vItm1 = vui8matrix_s(nrl, nrh, ncl, nch); //Creation d'une matrice SIMD avec les indices scalaires
    MatScal2MatSIMD(vItm1, Itm1, vi0, vi1, vj0, vj1);

    //routine_FrameDifference_SSE2(vIt, vItm1, vEt, vi0, vi1, vj0, vj1, seuilSIMD);

    fermeture3x3_SIMD(vItm1, vEt, vi0, vi1, vj0, vj1);
    MatSIMD2MatScal(vEt, Et, vi0, vi1, vj0, vj1);    //On fait la copie d'une matrice SIMD dans une image normale
    SavePGM_ui8matrix(Et, nrl, nrh, ncl, nch, "testmorpho/fermeture_produite.pgm");

    ouverture3x3_SIMD(vItm1, vEt, vi0, vi1, vj0, vj1);
    MatSIMD2MatScal(vEt, Et, vi0, vi1, vj0, vj1);    //On fait la copie d'une matrice SIMD dans une image normale
    SavePGM_ui8matrix(Et, nrl, nrh, ncl, nch, "testmorpho/ouverture_produite.pgm");



    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl, nrh, ncl, nch );
    free_vui8matrix(vItm1, vi0, vi1, vj0, vj1);
    free_vui8matrix(vEt, vi0, vi1, vj0, vj1);



}

void test_routine_FrameDifferenceMorpho3x3ouverture(int seuil)
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
    uint8 **Itm1 =  LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **It = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8 **Et1 = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);

    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);

        routine_FrameDifference(It, Itm1, Et, nrl,nrh,ncl,nch, seuil);
        sprintf(nomImageSave, "car3Frame3x3O/car_3%03d.pgm", i);
        //ouverture3x3(Et,Et1, nrl,nrh,ncl,nch);
        CHRONO(ouverture3x3(Et,Et1, nrl,nrh,ncl,nch), cycles);
        cycleTotal+=cycles;

        SavePGM_ui8matrix(Et1, nrl, nrh, ncl, nch, nomImageSave);
        copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
    }


    cycleTotal/=NBIMAGES;
    cycleTotal/=((nch+1)*(nrh+1));
    BENCH(printf("Cycles Ouverture = "));
    BENCH(printf(format, cycleTotal));

    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );
    free_ui8matrix(Et1, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );

}

void convFormatToBin(uint8 ** EtAvant, ulong32 **EtApres, long nrl, long nrh, long ncl, long nch)
{
    int i, j, k ;
    for (i = nrl; i <= nrh; i++)
    {
        k = 0;
        EtApres[i][k] = 0;//Comme ca on n'a qu'a mettre les pixels a 1 quand il faut

        for(j = ncl; j <= nch; j++)
        {
            if(j%TLONG == 0)
            {
                k++;
                EtApres[i][k] = 0;//Comme ca on n'a qu'a mettre les pixels a 1 quand il faut

            }
            if(EtAvant[i][j] == 255)
                EtApres[i][k] |= 1 << ((TLONG-1)-j%TLONG);//Commence a ecrire par la gauche pour avoir la meme representation que les pixels
        }
    }
}

void convFormatToChar(ulong32 ** EtAvant, uint8 **EtApres, long nrl, long nrh, long ncl, long nch)
{
    int i, j, k ;
    int val;
    for (i = nrl; i <= nrh; i++)
    {
        k = 0;
        for(j = ncl; j <= nch; j++)
        {
            EtApres[i][j] = 0;
            val = 0;
            if(j%TLONG == 0)
            {
                k++;
            }
            val = (EtAvant[i][k] >> ((TLONG-1)-j%TLONG) ) & 1;
            EtApres[i][j] = (val==1 ? 255:0);

        }
    }

}
void test_routine_FrameDifferenceMorpho3x3ouverture_bin(int seuil)
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
    uint8 **Itm1 =  LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **It = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8 **Et1 = ui8matrix(nrl, nrh, ncl, nch); //Pas besoin de bords car on va juste recopier dedans


    long bi0 = nrl, bi1 = nrh, bj0 = ncl/TLONG, bj1 = nch/TLONG;
    ulong32 **EtBin = long64matrix(bi0 - BORD, bi1 + BORD, bj0 - 1, bj1 + 1);//Ajout de bords et allocation de matrice etiquette binaire
    ulong32 **Et1Bin = long64matrix(bi0 - BORD, bi1 + BORD, bj0 - 1, bj1 + 1);
    ulong32 **O0 =  long64matrix(bi0 - BORD, bi1 + BORD, bj0 - 1, bj1 + 1);

    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);

        routine_FrameDifference(It, Itm1, Et, nrl,nrh,ncl,nch, seuil);
        sprintf(nomImageSave, "car3Frame3x3O_bin/car_3%03d.pgm", i);
        convFormatToBin(Et, EtBin, nrl, nrh, ncl, nch);
        //ouverture3x3(Et,Et1, nrl,nrh,ncl,nch);
        CHRONO(ouverture3x3_bin(EtBin,Et1Bin, O0, bi0,bi1,bj0,bj1), cycles);
        cycleTotal+=cycles;

        convFormatToChar(Et1Bin, Et1, nrl, nrh, ncl, nch);
        SavePGM_ui8matrix(Et1, nrl, nrh, ncl, nch, nomImageSave);
        copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
    }


    cycleTotal/=NBIMAGES;
    cycleTotal/=((nch+1)*(nrh+1));
    BENCH(printf("Cycles Ouverture_bin = "));
    BENCH(printf(format, cycleTotal));

    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );
    free_ui8matrix(Et1, nrl, nrh, ncl, nch );

    free_long64matrix(EtBin, bi0 - BORD, bi1 + BORD, bj0 - 1, bj1 + 1 );
    free_long64matrix(Et1Bin, bi0 - BORD, bi1 + BORD, bj0 - 1, bj1 + 1);
    free_long64matrix(O0, bi0 - BORD, bi1 + BORD, bj0 - 1, bj1 + 1);


}

void test_routine_FrameDifferenceMorpho3x3fermeture_bin(int seuil)
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
    uint8 **Itm1 =  LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **It = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8 **Et1 = ui8matrix(nrl, nrh, ncl, nch); //Pas besoin de bords car on va juste recopier dedans


    long bi0 = nrl, bi1 = nrh, bj0 = ncl/TLONG, bj1 = nch/TLONG;
    ulong32 **EtBin = long64matrix(bi0 - BORD, bi1 + BORD, bj0 - 1, bj1 + 1);//Ajout de bords et allocation de matrice etiquette binaire
    ulong32 **Et1Bin = long64matrix(bi0 - BORD, bi1 + BORD, bj0 - 1, bj1 + 1);
    ulong32 **O0 =  long64matrix(bi0 - BORD, bi1 + BORD, bj0 - 1, bj1 + 1);

    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);

        routine_FrameDifference(It, Itm1, Et, nrl,nrh,ncl,nch, seuil);
        sprintf(nomImageSave, "car3Frame3x3F_bin/car_3%03d.pgm", i);
        convFormatToBin(Et, EtBin, nrl, nrh, ncl, nch);
        //ouverture3x3(Et,Et1, nrl,nrh,ncl,nch);
        CHRONO(ouverture3x3_bin(EtBin,Et1Bin, O0, bi0,bi1,bj0,bj1), cycles);
        cycleTotal+=cycles;

        convFormatToChar(Et1Bin, Et1, nrl, nrh, ncl, nch);
        SavePGM_ui8matrix(Et1, nrl, nrh, ncl, nch, nomImageSave);
        copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
    }


    cycleTotal/=NBIMAGES;
    cycleTotal/=((nch+1)*(nrh+1));
    BENCH(printf("Cycles fermeture_bin = "));
    BENCH(printf(format, cycleTotal));

    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );
    free_ui8matrix(Et1, nrl, nrh, ncl, nch );

    free_long64matrix(EtBin, bi0 - BORD, bi1 + BORD, bj0 - 1, bj1 + 1 );
    free_long64matrix(Et1Bin, bi0 - BORD, bi1 + BORD, bj0 - 1, bj1 + 1);
    free_long64matrix(O0, bi0 - BORD, bi1 + BORD, bj0 - 1, bj1 + 1);


}

void test_routine_FrameDifferenceMorpho3x3fermeture(int seuil)
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
    uint8 **Itm1 =  LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **It = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8 **Et1 = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);

    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);

        routine_FrameDifference(It, Itm1, Et, nrl,nrh,ncl,nch, seuil);
        sprintf(nomImageSave, "car3Frame3x3F/car_3%03d.pgm", i);
        //fermeture3x3(Et,Et1, nrl,nrh,ncl,nch);
        CHRONO(fermeture3x3(Et,Et1, nrl,nrh,ncl,nch), cycles);
        cycleTotal+=cycles;

        SavePGM_ui8matrix(Et1, nrl, nrh, ncl, nch, nomImageSave);
        copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
    }

    cycleTotal/=NBIMAGES;
    cycleTotal/=((nch+1)*(nrh+1));
    BENCH(printf("Cycles Fermeture = "));
    BENCH(printf(format, cycleTotal));

    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );
    free_ui8matrix(Et1, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );

}
void test_routine_FrameDifferenceMorpho3x3ouvertureFermeture(int seuil)
{
    /*pour commit */
    char nomImageLoad[50];// = "car3/car_3";
    char nomImageSave[50];// = "car3Sigma/car_3"
    long nrl, nrh, ncl, nch;
    uint8 **Itm1 =  LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **It = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8 **Et1 = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);

    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);

        routine_FrameDifference(It, Itm1, Et, nrl,nrh,ncl,nch, seuil);
        sprintf(nomImageSave, "car3Frame3x3OF/car_3%03d.pgm", i);
        ouverture3x3(Et,Et1, nrl,nrh,ncl,nch);
        fermeture3x3(Et1,Et, nrl,nrh,ncl,nch);
        SavePGM_ui8matrix(Et, nrl, nrh, ncl, nch, nomImageSave);
        copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
    }

    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );
    free_ui8matrix(Et1, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );

}

void test_routine_FrameDifferenceMorpho3x3fermetureOuverture(int seuil)
{
    /*pour commit */
    char nomImageLoad[50];// = "car3/car_3";
    char nomImageSave[50];// = "car3Sigma/car_3"
    long nrl, nrh, ncl, nch;
    uint8 **Itm1 =  LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **It = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8 **Et1 = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);

    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);

        routine_FrameDifference(It, Itm1, Et, nrl,nrh,ncl,nch, seuil);
        sprintf(nomImageSave, "car3Frame3x3FO/car_3%03d.pgm", i);
        fermeture3x3(Et,Et1, nrl,nrh,ncl,nch);
        ouverture3x3(Et1,Et, nrl,nrh,ncl,nch);
        SavePGM_ui8matrix(Et, nrl, nrh, ncl, nch, nomImageSave);
        copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
    }

    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );
    free_ui8matrix(Et1, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );

}

void test_routine_FrameDifferenceMorpho3x3fermeturefermeture(int seuil)
{
    /*pour commit */
    char nomImageLoad[50];// = "car3/car_3";
    char nomImageSave[50];// = "car3Sigma/car_3"
    long nrl, nrh, ncl, nch;
    uint8 **Itm1 =  LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **It = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8 **Et1 = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);

    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);

        routine_FrameDifference(It, Itm1, Et, nrl,nrh,ncl,nch, seuil);
        sprintf(nomImageSave, "car3Frame3x3FF/car_3%03d.pgm", i);

        fermeture3x3(Et,Et1, nrl,nrh,ncl,nch);
        fermeture3x3(Et1,Et, nrl,nrh,ncl,nch);
        SavePGM_ui8matrix(Et, nrl, nrh, ncl, nch, nomImageSave);
        copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
    }

    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );
    free_ui8matrix(Et1, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );

}


void test_routine_FrameDifferenceMorpho3x3fermeture_pipe(int seuil)
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
    uint8 **Itm1 =  LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **It = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8 **Et1 = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8 **O0 = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);//Matrice intermediaire pour la morpho

    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);

        routine_FrameDifference(It, Itm1, Et, nrl,nrh,ncl,nch, seuil);
        sprintf(nomImageSave, "car3Frame3x3F_pipe/car_3%03d.pgm", i);
        //fermeture3x3(Et,Et1, nrl,nrh,ncl,nch);
        CHRONO(fermeture3x3_pipe(Et,Et1, O0, nrl,nrh,ncl,nch), cycles);
        cycleTotal+=cycles;

        SavePGM_ui8matrix(Et1, nrl, nrh, ncl, nch, nomImageSave);
        copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
    }

    cycleTotal/=NBIMAGES;
    cycleTotal/=((nch+1)*(nrh+1));
    BENCH(printf("Cycles FermeturePipe = "));
    BENCH(printf(format, cycleTotal));

    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );
    free_ui8matrix(Et1, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );
    free_ui8matrix(O0, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );


}

void test_routine_FrameDifferenceMorpho3x3ouverture_pipe(int seuil)
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
    uint8 **Itm1 =  LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **It = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8 **Et1 = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8 **O0 = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);

    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);

        routine_FrameDifference(It, Itm1, Et, nrl,nrh,ncl,nch, seuil);
        sprintf(nomImageSave, "car3Frame3x3O_pipe/car_3%03d.pgm", i);
        //ouverture3x3(Et,Et1, nrl,nrh,ncl,nch);
        CHRONO(ouverture3x3_pipe(Et,Et1, O0, nrl,nrh,ncl,nch), cycles);
        cycleTotal+=cycles;

        SavePGM_ui8matrix(Et1, nrl, nrh, ncl, nch, nomImageSave);
        copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
    }


    cycleTotal/=NBIMAGES;
    cycleTotal/=((nch+1)*(nrh+1));
    BENCH(printf("Cycles OuverturePipe = "));
    BENCH(printf(format, cycleTotal));

    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );
    free_ui8matrix(Et1, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );
    free_ui8matrix(O0, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );

}
