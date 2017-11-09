#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#include "nrdef.h"
#include "vnrutil.h"

#include "morpho_simd.h"


void erosion3x3_SIMD(vuint8 **It,vuint8 **It1,long vi0,long vi1,long vj0,long vj1)
{
	vuint8 l1;
	vuint8 l2;
	vuint8 l3;

	vuint8 result1,result2,result3;

	vuint8 temp = init_vuint8(255);
	vuint8 left,right;

	//vuint8 or_droit, or_gauche;

	// temp = init_vuint8(255); //AND DONC 255 INVARIANT
	//#CONNERIES
	// or_gauche = _mm_srli_si128(temp,15);
	// or_droit = _mm_slli_si128(temp,15);

	//premiere ligne : prologue
	int j = vj0;
	int i = vi0;

	l1 = _mm_load_si128(&It[i+0][j]);
	l2 = _mm_load_si128(&It[i+1][j]);

	result1 =  _mm_and_si128(l1,l2);

	right = _mm_srli_si128(result1,1);
	//right = _mm_or_si128(or_droit,right);

	//ici right contient result1 decallé vers la droite avec 255 comme valeur de gauche
	//
	j++;
	for(j;j<=vj1;j++)
	{
		l1 = _mm_load_si128(&It[i+0][j]);
		l2 = _mm_load_si128(&It[i+1][j]);

		result2 = _mm_and_si128(l1,l2);

		left = _mm_slli_si128(result1,1);
		left = _mm_or_si128(left,_mm_srli_si128(result2,15)); // On complete le vecteur gauche en ajoutant a sa doite la valeur de gauvhe du vecteur suivant

		result3 = _mm_and_si128(left,_mm_and_si128(right,result1)); // Valeur finale de IT1[i][j-1]

		right = _mm_srli_si128(result2,1);
		right = _mm_or_si128(right,_mm_slli_si128(result1,15));

		result1 = result2;

		_mm_store_si128(&It1[i][j-1], result3);
	}

	left = _mm_slli_si128(result1,1);
	//left = _mm_or_si128(left,or_gauche);

	result3 = _mm_and_si128(left,_mm_and_si128(right,result1)); // Valeur finale de IT1[i][j-1]

	_mm_store_si128(&It1[i][j-1],result3);

	j= vj0;
	i++;
	// corps de boucle
	for(i;i<vi1;i++)
	{
		l1 = _mm_load_si128(&It[i-1][j]);
		l2 = _mm_load_si128(&It[i+0][j]);
		l3 = _mm_load_si128(&It[i+1][j]);

		temp =  _mm_and_si128(l1,l2);
		result1 = _mm_and_si128(temp,l3);

		right = _mm_srli_si128(result1,1);

		//j vaut vj0 donc on insere un zero a gauche de right
		//right = _mm_or_si128(right,or_droit);
		
		j++;
		for(j;j<=vj1;j++)
		{
			l1 = _mm_load_si128(&It[i-1][j]);
			l2 = _mm_load_si128(&It[i+0][j]);
			l3 = _mm_load_si128(&It[i+1][j]);

			temp =  _mm_and_si128(l1,l2);
			result2 = _mm_and_si128(temp,l3);

			left = _mm_slli_si128(result1,1);
			left = _mm_or_si128(left,_mm_srli_si128(result2,15)); // On complete le vecteur gauche en ajoutant a sa doite la valeur de gauvhe du vecteur suivant

			result3 = _mm_and_si128(left,_mm_and_si128(right,result1));

			right = _mm_srli_si128(result2,1);
			right = _mm_or_si128(right,_mm_slli_si128(result1,15));

			result1 = result2;

			_mm_store_si128(&It1[i][j-1], result3);
		}

		left = _mm_slli_si128(result1,1);
		//left = _mm_or_si128(left,or_gauche);

		result3 = _mm_and_si128(left,_mm_and_si128(right,result1)); // Valeur finale de IT1[i][j-1]

		_mm_store_si128(&It1[i][j-1],result3);

		j= vj0;
	}

	l1 = _mm_load_si128(&It[i-1][j]);
	l2 = _mm_load_si128(&It[i+0][j]);

	result1 = _mm_and_si128(l1,l2);

	right = _mm_srli_si128(result1,1);

	//j vaut vj0 donc on insere un zero a gauche de right
	//right = _mm_or_si128(right,or_droit);
	
	j++;
	for(j;j<=vj1;j++)
	{
		l1 = _mm_load_si128(&It[i-1][j]);
		l2 = _mm_load_si128(&It[i+0][j]);

		result2 = _mm_and_si128(l1,l2);

		left = _mm_slli_si128(result1,1);
		left = _mm_or_si128(left,_mm_srli_si128(result2,15)); // On complete le vecteur gauche en ajoutant a sa doite la valeur de gauvhe du vecteur suivant

		result3 = _mm_and_si128(left,_mm_and_si128(right,result1)); // Valeur finale de IT1[i][j-1]

		right = _mm_srli_si128(result2,1);
		right = _mm_or_si128(right,_mm_slli_si128(result1,15));

		result1 = result2;

		_mm_store_si128(&It1[i][j-1], result3);
	}
	//derniere ligne


	left = _mm_slli_si128(result1,1);
	//left = _mm_or_si128(left,or_gauche);

	result3 = _mm_and_si128(left,_mm_and_si128(right,result1)); // Valeur finale de IT1[i][j-1]

	_mm_store_si128(&It1[i][j-1],result3);
}

void dilatation3x3_SIMD(vuint8 **It,vuint8 **It1,long vi0,long vi1,long vj0,long vj1)
{
	vuint8 l1;
	vuint8 l2;
	vuint8 l3;

	vuint8 result1,result2,result3;

	vuint8 temp = init_vuint8(255);
	vuint8 left,right;

	//vuint8 or_droit, or_gauche;

	temp = init_vuint8(255); //AND DONC 255 INVARIANT

	// or_gauche = _mm_andnot_si128( _mm_srli_si128(temp,15),init_vuint8(255));
	// or_droit = _mm_andnot_si128(  _mm_slli_si128(temp,15), init_vuint8(255));

	//premiere ligne : prologue
	int j = vj0;
	int i = vi0;

	l1 = _mm_load_si128(&It[i+0][j]);
	l2 = _mm_load_si128(&It[i+1][j]);

	result1 =  _mm_or_si128(l1,l2);

	right = _mm_srli_si128(result1,1);
	//right = _mm_or_si128(or_droit,right);

	//ici right contient result1 decallé vers la droite avec 255 comme valeur de gauche
	//
	j++;
	for(j;j<=vj1;j++)
	{
		l1 = _mm_load_si128(&It[i+0][j]);
		l2 = _mm_load_si128(&It[i+1][j]);

		result2 = _mm_or_si128(l1,l2);

		left = _mm_slli_si128(result1,1);
		left = _mm_or_si128(left,_mm_srli_si128(result2,15)); // On complete le vecteur gauche en ajoutant a sa doite la valeur de gauvhe du vecteur suivant

		result3 = _mm_or_si128(left,_mm_or_si128(right,result1)); // Valeur finale de IT1[i][j-1]

		right = _mm_srli_si128(result2,1);
		right = _mm_or_si128(right,_mm_slli_si128(result1,15));

		result1 = result2;

		_mm_store_si128(&It1[i][j-1], result3);
	}

	left = _mm_slli_si128(result1,1);
	//left = _mm_or_si128(left,or_gauche);

	result3 = _mm_or_si128(left,_mm_or_si128(right,result1)); // Valeur finale de IT1[i][j-1]

	_mm_store_si128(&It1[i][j-1],result3);

	j= vj0;
	i++;
	// corps de boucle
	for(i;i<vi1;i++)
	{
		l1 = _mm_load_si128(&It[i-1][j]);
		l2 = _mm_load_si128(&It[i+0][j]);
		l3 = _mm_load_si128(&It[i+1][j]);

		temp =  _mm_or_si128(l1,l2);
		result1 = _mm_or_si128(temp,l3);

		right = _mm_srli_si128(result1,1);

		//j vaut vj0 donc on insere un zero a gauche de right
		//right = _mm_or_si128(right,or_droit);
		
		j++;
		for(j;j<=vj1;j++)
		{
			l1 = _mm_load_si128(&It[i-1][j]);
			l2 = _mm_load_si128(&It[i+0][j]);
			l3 = _mm_load_si128(&It[i+1][j]);

			temp =  _mm_or_si128(l1,l2);
			result2 = _mm_or_si128(temp,l3);

			left = _mm_slli_si128(result1,1);
			left = _mm_or_si128(left,_mm_srli_si128(result2,15)); // On complete le vecteur gauche en ajoutant a sa doite la valeur de gauvhe du vecteur suivant

			result3 = _mm_or_si128(left,_mm_or_si128(right,result1));

			right = _mm_srli_si128(result2,1);
			right = _mm_or_si128(right,_mm_slli_si128(result1,15));

			result1 = result2;

			_mm_store_si128(&It1[i][j-1], result3);
		}

		left = _mm_slli_si128(result1,1);
		//left = _mm_or_si128(left,or_gauche);

		result3 = _mm_or_si128(left,_mm_or_si128(right,result1)); // Valeur finale de IT1[i][j-1]

		_mm_store_si128(&It1[i][j-1],result3);

		j= vj0;
	}

	l1 = _mm_load_si128(&It[i-1][j]);
	l2 = _mm_load_si128(&It[i+0][j]);

	result1 = _mm_or_si128(l1,l2);

	right = _mm_srli_si128(result1,1);

	//j vaut vj0 donc on insere un zero a gauche de right
	//right = _mm_or_si128(right,or_droit);
	
	j++;
	for(j;j<=vj1;j++)
	{
		l1 = _mm_load_si128(&It[i-1][j]);
		l2 = _mm_load_si128(&It[i+0][j]);

		result2 = _mm_or_si128(l1,l2);

		left = _mm_slli_si128(result1,1);
		left = _mm_or_si128(left,_mm_srli_si128(result2,15)); // On complete le vecteur gauche en ajoutant a sa doite la valeur de gauvhe du vecteur suivant

		result3 = _mm_or_si128(left,_mm_or_si128(right,result1)); // Valeur finale de IT1[i][j-1]

		right = _mm_srli_si128(result2,1);
		right = _mm_or_si128(right,_mm_slli_si128(result1,15));

		result1 = result2;

		_mm_store_si128(&It1[i][j-1], result3);
	}
	//derniere ligne


	left = _mm_slli_si128(result1,1);
	//left = _mm_or_si128(left,or_gauche);

	result3 = _mm_or_si128(left,_mm_or_si128(right,result1)); // Valeur finale de IT1[i][j-1]

	_mm_store_si128(&It1[i][j-1],result3);
}