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
	vuint8 ** vXOt = vui8matrix(vi0, vi1, vj0, vj1);
	vuint8 tmpIt;
	vuint8 tmpItm1;
	vuint8 tmpOt;
	vuint8 pixelNoir = init_vuint8(0);
	vuint8 pixelBlanc = init_vuint8(255);
	vuint8 tmpEt;
	vuint8 maxInt8 = init_vuint8(127);
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
    for(int i = vi0; i <= vi1; i++ )
    {
    	for(int j = vj0; j <= vj1; j++)
    	{
    		tmpOt = _mm_load_si128(&vXOt[i][j]);
            vuint8 res = _mm_cmplt_epi8(_mm_sub_epi8(tmpOt, maxInt8), _mm_sub_epi8(seuil, maxInt8)); //Met 1 si inferieur au seuil et 0 si superieur
            //display_vuint8(res," %d ", "Res");
            //vuint8 dest = _mm_or_si128(_mm_and_si128(res, pixelNoir), _mm_andnot_si128(res, pixelBlanc)); //Pixel noir si res a 1 et pixel blanc si res a 0
            res = _mm_andnot_si128(res, pixelBlanc);
            _mm_store_si128(&Et[i][j], res);

        }
    }
    free_vui8matrix(vXOt, vi0, vi1, vj0, vj1);
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
	vuint8 ** vXOt = vui8matrix(vi0, vi1, vj0, vj1);
	vuint8 tmpIt, tmpMt, tmpVt;
	vuint8 tmpItm1, tmpMtm1, tmpVtm1;
	vuint8 tmpOt;
	vuint8 pixelNoir = init_vuint8(0);
	vuint8 pixelBlanc = init_vuint8(255);
	vuint8 tmpEt;
	vuint8 un = init_vuint8(1);
	vuint8 VMAXSIMD = init_vuint8(VMAX);
	vuint8 VMINSIMD = init_vuint8(VMIN);
    vuint8 maxInt8 = init_vuint8(127);
    //Les comparaisons se font en signe donc il faut sub 127 pour que ça fasse une comparaison correcte
    //255 devient 127, 127 => 0 et 0 => -127
    //Si on ne fait pas ca, on a 127 < 128 qui est faux

    for(int i = vi0; i <= vi1; i++ )//Step1
    {
    	for(int j = vj0; j <= vj1; j++)
    	{


    		tmpMtm1 = _mm_load_si128(&Mtm1[i][j]);
    		tmpIt = _mm_load_si128(&It[i][j]);
    		vuint8 Mtm1Plus1 = _mm_add_epi8(tmpMtm1, un);
    		vuint8 Mtm1Moins1 = _mm_sub_epi8(tmpMtm1, un);

    		vuint8 res = _mm_cmplt_epi8(_mm_sub_epi8(tmpMtm1, maxInt8), _mm_sub_epi8(tmpIt, maxInt8));
            tmpMt = _mm_or_si128(_mm_and_si128(res, Mtm1Plus1), _mm_andnot_si128(res, tmpMtm1)); //Mtm1< It

            res = _mm_cmpgt_epi8(_mm_sub_epi8(tmpMtm1, maxInt8), _mm_sub_epi8(tmpIt, maxInt8));
            tmpMt = _mm_or_si128(_mm_and_si128(res, Mtm1Moins1), _mm_andnot_si128(res, tmpMt)); // //Mtm1 > It
            _mm_store_si128(&Mt[i][j], tmpMt);
        }
    }

    for(int i = vi0; i <= vi1; i++ )//Step1
    {
    	for(int j = vj0; j <= vj1; j++)
    	{
    		tmpIt = _mm_load_si128(&It[i][j]);
    		tmpMt = _mm_load_si128(&Mt[i][j]);
            tmpOt = _mm_min_epu8(_mm_sub_epi8(tmpMt,tmpIt), _mm_sub_epi8(tmpIt, tmpMt) ); //min(a-b,b-a) donne la valeur absolue car on peut pas avoir de valeurs negatives
            _mm_store_si128(&vXOt[i][j], tmpOt); //Sauvegarde de l'image de difference
        }
    }

    for(int i = vi0; i <= vi1; i++ )//Step1
    {
    	for(int j = vj0; j <= vj1; j++)
    	{
    		tmpOt = _mm_load_si128(&vXOt[i][j]);
    		tmpVtm1 = _mm_load_si128(&Vtm1[i][j]);

    		vuint8 NfoisOt = init_vuint8(0);
    		for(int k = 0; k < N; k++)
    		{
    			NfoisOt = _mm_adds_epi8(NfoisOt, tmpOt);
    		}

    		vuint8 Vtm1Plus1 = _mm_add_epi8(tmpVtm1, un);
    		vuint8 Vtm1Moins1 = _mm_sub_epi8(tmpVtm1, un);
    		vuint8 res = _mm_cmplt_epi8(_mm_sub_epi8(tmpVtm1, maxInt8), _mm_sub_epi8(NfoisOt, maxInt8));
            tmpVt = _mm_or_si128(_mm_and_si128(res, Vtm1Plus1), _mm_andnot_si128(res, tmpVtm1)); //Vtm1< N*Ot

            res = _mm_cmpgt_epi8(_mm_sub_epi8(tmpVtm1, maxInt8), _mm_sub_epi8(NfoisOt, maxInt8));
            tmpVt = _mm_or_si128(_mm_and_si128(res, Vtm1Moins1), _mm_andnot_si128(res, tmpVt)); // //Vtm1 > N* Ot


            tmpVt = _mm_max_epu8(_mm_min_epu8(tmpVt, VMAXSIMD), VMINSIMD);

            _mm_store_si128(&Vt[i][j], tmpVt);

        }
    }

    for(int i = vi0; i <= vi1; i++ )
    {
    	for(int j = vj0; j <= vj1; j++)
    	{
    		tmpOt = _mm_load_si128(&vXOt[i][j]);
    		tmpVt = _mm_load_si128(&Vt[i][j]);
            vuint8 res = _mm_cmplt_epi8(_mm_sub_epi8(tmpOt,maxInt8), _mm_sub_epi8(tmpVt,maxInt8)); //Met 1 si inferieur a Vt et 0 si superieur

            vuint8 dest = _mm_or_si128(_mm_and_si128(res, pixelNoir), _mm_andnot_si128(res, pixelBlanc)); //Pixel noir si res a 1 et pixel blanc si res a 0
            _mm_store_si128(&Et[i][j], dest);

        }
    }
    free_vui8matrix(vXOt, vi0, vi1, vj0, vj1);


}
