#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vnrutil.h"
#include "nrutil.h"
#include "morpho.h"
#include "morpho_simd.h"
#include "test_mouvement_SSE2.h"
#include "test_morpho.h"


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