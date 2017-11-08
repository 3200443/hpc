#include <stdio.h>
#include <stdlib.h>

#include "nrdef.h"
#include "nrutil.h"
#include "matric_roc.h"


void creation_matrices_ROC(char fichierVerite[], char fichierTest[])
{
    long nrl, nrh, ncl, nch;
    uint8 **ImageVerite =  LoadPGM_ui8matrix(fichierVerite, &nrl, &nrh, &ncl, &nch);
    uint8 **ImageTest =  LoadPGM_ui8matrix(fichierTest, &nrl, &nrh, &ncl, &nch);
    int matRoc[2][2] = {0};

    for(int i = nrl; i <= nrh; i++)
    {
        for(int j = ncl; j <= nch; j++)
        {
            if(ImageTest[i][j] == 255 && ImageVerite[i][j] == 255)
                matRoc[0][0]+=1; //VP
            else if(ImageTest[i][j] == 255 && ImageVerite[i][j] == 0)
                matRoc[0][1]+=1; //FN
            else if(ImageTest[i][j] == 0 && ImageVerite[i][j] == 255)
                matRoc[1][0]+=1; //FP
            else if(ImageTest[i][j] == 0 && ImageVerite[i][j] == 0)
                matRoc[1][1]+=1;
        }
    }
    printf("Matrice ROC = \n%d %d\n%d %d\n",matRoc[0][0], matRoc[0][1], matRoc[1][0], matRoc[1][1]);
}
