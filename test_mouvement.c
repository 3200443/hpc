#include <stdio.h>
#include <stdlib.h>
#include "nrdef.h"
#include "nrutil.h"
#include <string.h>
#include "mouvement.h"
#define NBIMAGES 199
//I0 = It et I1 = It-1 : pareil pour tout

void test_routine_sigmaDelta()
{
    //m[nrl..nrh][ncl..nch]
    char nomImageLoad[50];// = "car3/car_3";
    char nomImageSave[50];// = "car3Sigma/car_3"
    long nrl, nrh, ncl, nch;


    sprintf(nomImageLoad,"car3/car_3000.pgm");//Image a t-1
    uint8 **I1 =  LoadPGM_ui8matrix(nomImageLoad, &nrl, &nrh, &ncl, &nch);
    uint8 **I0 = ui8matrix(nrl, nrh, ncl, nch);

    int nrow=nrh-nrl+1,ncol=nch-ncl+1;

    uint8 **E0 = ui8matrix(nrl, nrh, ncl, nch);


    uint8 **V0 = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **V1 =  ui8matrix(nrl, nrh, ncl, nch);

    uint8 **M0 = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **M1 =  ui8matrix(nrl, nrh, ncl, nch);

    routine_SigmaDelta_step0(I1, M1, V1, nrl, nrh, ncl, nch);

    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad,"car3/car_3%03d.pgm",i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, I0);
        routine_SigmaDelta_1step(I0, I1, V0, V1, M0, M1, E0, nrl, nrh, ncl, nch);
        sprintf(nomImageSave,"car3Sigma/car_3%03d.pgm",i);
        SavePGM_ui8matrix(E0, nrl, nrh, ncl, nch, nomImageSave);//Copie de t a t-1
        memcpy(M1[nrl], M0[nrl], sizeof(uint8)*(nrow*ncol));
        memcpy(V1[nrl], V0[nrl], sizeof(uint8)*(nrow*ncol));
        memcpy(I1[nrl], I0[nrl], sizeof(uint8)*(nrow*ncol));

    }
    free_ui8matrix(E0, nrl, nrh, ncl, nch);
    free_ui8matrix(I0, nrl, nrh, ncl, nch);
    free_ui8matrix(I1, nrl, nrh, ncl, nch);
    free_ui8matrix(V0, nrl, nrh, ncl, nch);
    free_ui8matrix(V1, nrl, nrh, ncl, nch);
    free_ui8matrix(M0, nrl, nrh, ncl, nch);
    free_ui8matrix(M1, nrl, nrh, ncl, nch);





}

void test_routine_FrameDifference(int seuil)
{
    char nomImageLoad[50];// = "car3/car_3";
    char nomImageSave[50];// = "car3Sigma/car_3"
    long nrl, nrh, ncl, nch;
    uint8 **I1 =  LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **I0 = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **E0 = ui8matrix(nrl, nrh, ncl, nch);

    int nrow=nrh-nrl+1,ncol=nch-ncl+1;

    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, I0);

        routine_FrameDifference(I0, I1, E0, nrl,nrh,ncl,nch, seuil);
        sprintf(nomImageSave, "car3Frame/car_3%03d.pgm", i);
        SavePGM_ui8matrix(E0, nrl, nrh, ncl, nch, nomImageSave);
        memcpy(I1[nrl], I0[nrl], sizeof(uint8)*(nrow*ncol));
    }
    free_ui8matrix(I0, nrl, nrh, ncl, nch );
    free_ui8matrix(I1, nrl, nrh, ncl, nch );
    free_ui8matrix(E0, nrl, nrh, ncl, nch );
}


int main(int argc, char* argv[])
{

    test_routine_FrameDifference(atoi(argv[1]));
    test_routine_sigmaDelta();





    return 0;
}
