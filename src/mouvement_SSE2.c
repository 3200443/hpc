#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#include "nrdef.h"
#include "vnrutil.h"




void routine_FrameDifference_SSE2(vuint8 **It, vuint8 **Itm1, vuint8 **Et, long vi0,long vi1,long vj0,long vj1, vuint8 seuil)
{
    vuint8 ** vXOt = vui8matrix(vi0, vi1, vj0, vj1);
    vuint8 tmpIt;
    vuint8 tmpItm1;
    vuint8 tmpOt;
    vuint8 pixelNoir = init_vuint8(0);
    vuint8 pixelBlanc = init_vuint8(255);
    vuint8 tmpEt;
    int test = 0;
    for(int i = vi0; i <= vi1; i++ )
    {
        for(int j = vj0; j <= vj1; j++)
        {
            tmpIt = _mm_load_si128(&It[i][j]);
            tmpItm1 = _mm_load_si128(&Itm1[i][j]);
            tmpOt = _mm_min_epu8(_mm_sub_epi8(tmpIt,tmpItm1), _mm_sub_epi8(tmpItm1, tmpIt) ); //min(a-b,b-a) donne la valeur absolue car on peut pas avoir de valeurs negatives

            _mm_store_si128(&vXOt[i][j], tmpOt); //Sauvegarde de l'image de difference

        }
    }
    test = 0;
    for(int i = vi0; i <= vi1; i++ )
    {
        for(int j = vj0; j <= vj1; j++)
        {
            tmpOt = _mm_load_si128(&vXOt[i][j]);
            vuint8 res = _mm_cmplt_epi8(tmpOt, seuil); //Met 1 si inferieur au seuil et 0 si superieur

            vuint8 dest = _mm_or_si128(_mm_and_si128(res, pixelNoir), _mm_andnot_si128(res, pixelBlanc)); //Pixel noir si res a 1 et pixel blanc si res a 0
            _mm_store_si128(&Et[i][j], dest);

        }
    }
    free_vui8matrix(vXOt, vi0, vi1, vj0, vj1);





}
