#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#include "nrdef.h"
#include "vnrutil.h"

#define N 2
#define VMIN 20
#define VMAX 255 //V est entre 2 et 2^m-1 avec m le nombre de bits des donnees ici 8 => https://hal.inria.fr/hal-01130889/document
#define VINI 35


void routine_FrameDifference_SSE2(vuint8 **It, vuint8 **Itm1, vuint8 **Et, long vi0,long vi1,long vj0,long vj1, vuint8 seuil)
{
    vuint8 tmpIt;
    vuint8 tmpItm1;
    vuint8 tmpOt;
    vuint8 pixelBlanc = init_vuint8(255);
    vuint8 tmpEt;
    vuint8 maxSChar = init_vuint8(128);
    for(int i = vi0; i <= vi1; i++ )
    {
        for(int j = vj0; j <= vj1; j++)
        {
            //Calcul de Ot, image de difference
            tmpIt = _mm_load_si128(&It[i][j]);
            tmpItm1 = _mm_load_si128(&Itm1[i][j]);
            vuint8 max = _mm_max_epu8(tmpIt,tmpItm1);
            vuint8 min = _mm_min_epu8(tmpIt, tmpItm1);
            tmpOt = _mm_sub_epi8(max, min);

            //Si Ot < 0, on a 255, sinon 0 donc on inverse pour avoir 255 sur Et quand Ot>=0 et 0 pour Ot < 0
            vuint8 res = _mm_cmplt_epi8(_mm_sub_epi8(tmpOt, maxSChar), _mm_sub_epi8(seuil, maxSChar)); //Met 1 si inferieur au seuil et 0 si superieur

            res = _mm_andnot_si128(res, pixelBlanc);
            _mm_store_si128(&Et[i][j], res);

        }
    }

}


void routine_SigmaDelta_step0SSE2(vuint8** I, vuint8 **M, vuint8 **V, long vi0, long vi1, long vj0, long vj1)
{
    vuint8 tmpM;
    vuint8 tmpV;
    vuint8 tmpI;
    vuint8 ecartTypeIni = init_vuint8(VINI);
    for(int i = vi0; i <= vi1; i++ )
    {
        for(int j = vj0; j <= vj1; j++)
        {

            tmpI = _mm_load_si128(&I[i][j]); //M[i][j] = I[i][j];
            _mm_store_si128(&M[i][j], tmpI);
            _mm_store_si128(&V[i][j], ecartTypeIni);    //V[i][j] = VINI; //Au depart a VMIN mais il y avait beaucoup de mouvement des le debut, a 10 ça marche mieux
        }
    }
}



void routine_SigmaDelta_1stepSSE2(vuint8 **It, vuint8 **Itm1, vuint8**Vt, vuint8 **Vtm1, vuint8**Mt, vuint8 **Mtm1, vuint8 **Et,  long vi0, long vi1, long vj0, long vj1 )
{
    vuint8 tmpIt, tmpMt, tmpVt;
    vuint8 tmpItm1, tmpMtm1, tmpVtm1;
    vuint8 tmpOt;
    vuint8 pixelBlanc = init_vuint8(255);
    vuint8 tmpEt;
    vuint8 un = init_vuint8(1);
    vuint8 VMAXSIMD = init_vuint8(VMAX);
    vuint8 VMINSIMD = init_vuint8(VMIN);
    vuint8 maxSChar = init_vuint8(128);
    //Les comparaisons se font en signe donc il faut sub 128 pour que ça fasse une comparaison correcte
    //255 devient 128, 128 => 0 et 0 => -128
    //Si on ne fait pas ca, on a 128 < 128 qui est faux

    for(int i = vi0; i <= vi1; i++ )
    {
        for(int j = vj0; j <= vj1; j++)
        {

            //Step1
            tmpMtm1 = _mm_load_si128(&Mtm1[i][j]);
            tmpIt = _mm_load_si128(&It[i][j]);
            tmpVtm1 = _mm_load_si128(&Vtm1[i][j]);

            vuint8 Mtm1Plus1 = _mm_add_epi8(tmpMtm1, un);
            vuint8 Mtm1Moins1 = _mm_sub_epi8(tmpMtm1, un);
            vuint8 NfoisOt = init_vuint8(0);

            vuint8 res = _mm_cmplt_epi8(_mm_sub_epi8(tmpMtm1, maxSChar), _mm_sub_epi8(tmpIt, maxSChar));
            tmpMt = _mm_or_si128(_mm_and_si128(res, Mtm1Plus1), _mm_andnot_si128(res, tmpMtm1)); //Mtm1< It

            res = _mm_cmpgt_epi8(_mm_sub_epi8(tmpMtm1, maxSChar), _mm_sub_epi8(tmpIt, maxSChar));
            tmpMt = _mm_or_si128(_mm_and_si128(res, Mtm1Moins1), _mm_andnot_si128(res, tmpMt)); // //Mtm1 > It


            //Step 2 Calcul matrice difference |Mt-It|
            vuint8 max = _mm_max_epu8(tmpIt,tmpMt);
            vuint8 min = _mm_min_epu8(tmpIt, tmpMt);
            tmpOt = _mm_sub_epi8(max, min); //Le max - min donne la valeur absolue
            //Step 3 Vt Update and clamping
            for(int k = 0; k < N; k++)
            {
                NfoisOt = _mm_adds_epi8(NfoisOt, tmpOt);
            }

            vuint8 Vtm1Plus1 = _mm_add_epi8(tmpVtm1, un);
            vuint8 Vtm1Moins1 = _mm_sub_epi8(tmpVtm1, un);

            res = _mm_cmplt_epi8(_mm_sub_epi8(tmpVtm1, maxSChar), _mm_sub_epi8(NfoisOt, maxSChar));//On soustrait 128 car la comparaison est signee
            tmpVt = _mm_or_si128(_mm_and_si128(res, Vtm1Plus1), _mm_andnot_si128(res, tmpVtm1)); //Vtm1< N*Ot

            res = _mm_cmpgt_epi8(_mm_sub_epi8(tmpVtm1, maxSChar), _mm_sub_epi8(NfoisOt, maxSChar));
            tmpVt = _mm_or_si128(_mm_and_si128(res, Vtm1Moins1), _mm_andnot_si128(res, tmpVt)); // //Vtm1 > N* Ot

            tmpVt = _mm_max_epu8(_mm_min_epu8(tmpVt, VMAXSIMD), VMINSIMD);

            //Step 4: Et estimation
            res = _mm_cmplt_epi8(_mm_sub_epi8(tmpOt,maxSChar), _mm_sub_epi8(tmpVt,maxSChar)); //Met 255 si inferieur a Vt et 0 si superieur

            vuint8 dest = _mm_andnot_si128(res, pixelBlanc);//Inverse les 255 et 0 pour avoir la bonne couleur de pixel

            _mm_store_si128(&Et[i][j], dest);
            _mm_store_si128(&Vt[i][j], tmpVt);
            _mm_store_si128(&Mt[i][j], tmpMt);
        }
    }

}
