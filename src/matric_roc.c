#include <stdio.h>
#include <stdlib.h>

#include "nrdef.h"
#include "nrutil.h"
#include "matric_roc.h"
#define IMAGEMIN 121
#define IMAGEMAX 149

void creation_matrices_ROC(char nomDossier[])
{
	long nrl, nrh, ncl, nch;
	char nomImageLoad[100];
	uint8 **ImageVerite =  LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
	uint8 **ImageTest = ui8matrix(nrl, nrh, ncl, nch);
	int matRoc[2][2] = {0};
	for(int i = IMAGEMIN; i <= IMAGEMAX; i++)
	{
        sprintf(nomImageLoad,"verite/car_3%03d.pgm",i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, ImageVerite);
        sprintf(nomImageLoad,"%s/car_3%03d.pgm",nomDossier, i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, ImageTest);


        for(int i = nrl; i <= nrh; i++)
        {
        	for(int j = ncl; j <= nch; j++)
        	{
        		if(ImageTest[i][j] == 255 && ImageVerite[i][j] == 255)
                    matRoc[0][0]+=1; //VP
                else if(ImageTest[i][j] == 0 && ImageVerite[i][j] == 255)
                    matRoc[0][1]+=1; //FN
                else if(ImageTest[i][j] == 255 && ImageVerite[i][j] == 0)
                    matRoc[1][0]+=1; //FP
                else if(ImageTest[i][j] == 0 && ImageVerite[i][j] == 0)
                    matRoc[1][1]+=1; //VN
            }
        }
    }


    double rapport = 1.0*(matRoc[0][0]+matRoc[1][1])/(matRoc[0][1]+matRoc[1][0]);
    printf("Matrice ROC = \n%d %d\n%d %d\n",matRoc[0][0], matRoc[0][1], matRoc[1][0], matRoc[1][1]);
    printf("Rapport = %d/%d = %f\n",matRoc[1][1]+matRoc[0][0], matRoc[0][1]+matRoc[1][0],rapport);
    free_ui8matrix(ImageVerite, nrl, nrh, ncl, nch);
    free_ui8matrix(ImageTest, nrl, nrh, ncl, nch);
}
