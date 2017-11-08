#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nrdef.h"
#include "nrutil.h"
#include "morpho_simd.h"


void erosion3x3_SIMD(vuint8 **It0,vuint8 **It1,long vi0,long vi1,long vj0,long vj1)
{
	vuint8 l1;
	vuint8 l2;
	vuint8 l3;

	vuint8 temp,result1,result2,result3;

	vuint8 left, right;

	vuint8 or_droit, or_gauche;

	temp = init_vuint8(255); //AND DONC 255 INVARIANT
	or_gauche = _mm_srli_epi16(temp,15);
	or_droit = _mm_slli_epi16(temp,15);

	//premiere ligne : prologue
	int j = vj0;
	int i = vi0;

	l1 = _mm_load_si128(&It[i+0][j]);
	l2 = _mm_load_si128(&It[i+1][j]);

	result1 =  _mm_and_si128(l1,l2);

	left = _mm_slli_epi16(result1,1);
	right = _mm_srli_epi16(result1,1);


	right = _mm_or_si128(or_droit,right);

	//ici right contient result1 decallé vers la droite avec 255 comme valeur de gauche
	for(j;j<vj1;j++)
	{
		l1 = _mm_load_si128(&It[i+0][j]);
		l2 = _mm_load_si128(&It[i+1][j]);

		result2 = _mm_and_si128(l1,l2);

		left = _mm_or_si128(left,_mm_srli_epi16(result2,15)); // On complete le vecteur gauche en ajoutant a sa doite la valeur de gauvhe du vecteur suivant

		result3 = _mm_and_si128(left,_mm_and_si128(right,result1)); // Valeur finale de IT1[i][j-1]

		right = _mm_or_si128(result2,_mm_slli_epi16(result1,15));

		result1 = result2;

		_mm_store_si128(&It[i][j-1], result3);
	}

	left = _mm_or_si128(result2,or_gauche);


	_mm_store_si128(&result3);

	j= vj0;

	l1 = _mm_load_si128(&It[i+0][j]);
	l2 = _mm_load_si128(&It[i+1][j]);
	result1 =  _mm_and_si128(l1,l2);
	i++;
	// corps de boucle
	for(i;i<vj1-1;i++)
	{
		l3 = _mm_load_si128(&It[i+1][j]);
		result2 = _mm_and_si128(result2,l3);
		left = _mm_or_si128(left,_mm_srli_epi16(result2,15)); // On complete le vecteur gauche en ajoutant a sa doite la valeur de gauvhe du vecteur suivant

	}








	int j = vj0;
	int i = vi0+1;

	//vecteur haut gauche :
	//
	l1 = _mm_load_si128(&It[i-1][j]);
	l2 = _mm_load_si128(&It[i+0][j]);
	l3 = _mm_load_si128(&It[i+1][j]); //On charge les 3 premier vecteurs


	result1 = _mm_and_si128(l1,l2);
	result2 = _mm_and_si128(result1,l3); // On fait le and verticalement

	left = _mm_slli_epi16(result2,1);
	right = _mm_srli_epi16(result2,1);

	result1 = init_vuint8(255); //AND DONC 255 INVARIANT
	result1 = _mm_slli_epi16(result1,15);
	right = _mm_or_si128(result1,right);





	for(i; i < vi1; i++ )
	{




		for(j; j< vj1; j++)
		{






			l1 = _mm_load_si128(&It[i-1][j]);
			l2 = _mm_load_si128(&It[i+0][j]);
			l3 = _mm_load_si128(&It[i+1][j]);

			result1 = _mm_and_si128(l1,l2);
			result3 = _mm_and_si128(result1,l2); //result2 contient tous les "and" verticaux

			//ici result 2 contient le vecteur precedent, result3 le vecteur courant il faut completer left
			result1 = _mm_srli_epi16(result3,15);
			right = _mm_or_si128(result1,right);

			left = _mm_slli_epi16(result2,1);
			right = _mm_srli_epi16(result2,1);


		}
	}
	




}