#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#include "nrdef.h"
#include "vnrutil.h"

#include "morpho_simd.h"

//ui8matrix();

void erosion3x3_SIMD(vuint8 **It,vuint8 **It1,long vi0,long vi1,long vj0,long vj1)
{
	vuint8 l1;
	vuint8 l2;
	vuint8 l3;

	vuint8 result1,result2,result3;

	vuint8 temp;
	vuint8 left,right;

	//premiere ligne : prologue
	int j = vj0;
	int i = vi0;

	l1 = _mm_load_si128(&It[i+0][j]);
	l2 = _mm_load_si128(&It[i+1][j]);

	result1 = _mm_and_si128(l1,l2);

	left = _mm_slli_si128(result1,1);

	//ici right contient result1 decallé vers la droite avec 255 comme valeur de gauche
	//
	j++;
	for(j;j<=vj1;j++)
	{
		l1 = _mm_load_si128(&It[i+0][j]);
		l2 = _mm_load_si128(&It[i+1][j]);

		result2 = _mm_and_si128(l1,l2);

		right = _mm_srli_si128(result1,1);
		right = _mm_or_si128(right,_mm_slli_si128(result2,15)); // On complete le vecteur gauche en ajoutant a sa doite la valeur de gauvhe du vecteur suivant

		result3 = _mm_and_si128(left,_mm_and_si128(right,result1)); // Valeur finale de IT1[i][j-1]

		left = _mm_slli_si128(result2,1);
		left = _mm_or_si128(left,_mm_srli_si128(result1,15));

		result1 = result2;

		_mm_store_si128(&It1[i][j-1], result3);
	}

	right = _mm_srli_si128(result1,1);

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

		left = _mm_slli_si128(result1,1);

		//j vaut vj0 donc on insere un zero a gauche de right
		
		j++;
		for(j;j<=vj1;j++)
		{
			l1 = _mm_load_si128(&It[i-1][j]);
			l2 = _mm_load_si128(&It[i+0][j]);
			l3 = _mm_load_si128(&It[i+1][j]);

			temp =  _mm_and_si128(l1,l2);
			result2 = _mm_and_si128(temp,l3);

			right = _mm_srli_si128(result1,1);
			right = _mm_or_si128(right,_mm_slli_si128(result2,15)); // On complete le vecteur gauche en ajoutant a sa doite la valeur de gauvhe du vecteur suivant

			result3 = _mm_and_si128(left,_mm_and_si128(right,result1));

			left = _mm_slli_si128(result2,1);
			left = _mm_or_si128(left,_mm_srli_si128(result1,15));

			result1 = result2;

			_mm_store_si128(&It1[i][j-1], result3);
		}

		right = _mm_srli_si128(result1,1);

		result3 = _mm_and_si128(left,_mm_and_si128(right,result1)); // Valeur finale de IT1[i][j-1]

		_mm_store_si128(&It1[i][j-1],result3);

		j= vj0;
	}

	l1 = _mm_load_si128(&It[i-1][j]);
	l2 = _mm_load_si128(&It[i+0][j]);

	result1 = _mm_and_si128(l1,l2);

	left = _mm_slli_si128(result1,1);

	//j vaut vj0 donc on insere un zero a gauche de right
	
	j++;
	for(j;j<=vj1;j++)
	{
		l1 = _mm_load_si128(&It[i-1][j]);
		l2 = _mm_load_si128(&It[i+0][j]);

		result2 = _mm_and_si128(l1,l2);

		right = _mm_srli_si128(result1,1);
		right = _mm_or_si128(right,_mm_slli_si128(result2,15)); // On complete le vecteur gauche en ajoutant a sa doite la valeur de gauvhe du vecteur suivant

		result3 = _mm_and_si128(left,_mm_and_si128(right,result1)); // Valeur finale de IT1[i][j-1]

		left = _mm_slli_si128(result2,1);
		left = _mm_or_si128(left,_mm_srli_si128(result1,15));

		result1 = result2;

		_mm_store_si128(&It1[i][j-1], result3);
	}

	right = _mm_srli_si128(result1,1);

	result3 = _mm_and_si128(left,_mm_and_si128(right,result1)); // Valeur finale de IT1[i][j-1]

	_mm_store_si128(&It1[i][j-1],result3);
}

void erosion3x3_SIMD_B(vuint8 **It,vuint8 **It1,long vi0,long vi1,long vj0,long vj1)
{
	vuint8 l1;
	vuint8 l2;
	vuint8 l3;

	vuint8 result1,result2,result3;

	vuint8 temp;
	vuint8 left,right;

	//premiere ligne : prologue
	int j = vj0;
	int i = vi0;

	// corps de boucle
	for(i;i<=vi1;i++)
	{
		l1 = _mm_load_si128(&It[i-1][j]);
		l2 = _mm_load_si128(&It[i+0][j]);
		l3 = _mm_load_si128(&It[i+1][j]);

		temp =  _mm_and_si128(l1,l2);
		result1 = _mm_and_si128(temp,l3);

		left = _mm_slli_si128(result1,1);

		//j vaut vj0 donc on insere un zero a gauche de right
		
		j++;
		for(j;j<=vj1;j++) //densité arithmétique de 2.5 : (10)/4
		{
			l1 = _mm_load_si128(&It[i-1][j]);
			l2 = _mm_load_si128(&It[i+0][j]);
			l3 = _mm_load_si128(&It[i+1][j]);

			temp =  _mm_and_si128(l1,l2);
			result2 = _mm_and_si128(temp,l3);

			right = _mm_srli_si128(result1,1);
			right = _mm_or_si128(right,_mm_slli_si128(result2,15)); // On complete le vecteur gauche en ajoutant a sa doite la valeur de gauvhe du vecteur suivant

			result3 = _mm_and_si128(left,_mm_and_si128(right,result1));

			left = _mm_slli_si128(result2,1);
			left = _mm_or_si128(left,_mm_srli_si128(result1,15));

			result1 = result2;

			_mm_store_si128(&It1[i][j-1], result3);
		}

		right = _mm_srli_si128(result1,1);

		result3 = _mm_and_si128(left,_mm_and_si128(right,result1)); // Valeur finale de IT1[i][j-1]

		_mm_store_si128(&It1[i][j-1],result3);

		j= vj0;
	}
}

void dilatation3x3_SIMD(vuint8 **It,vuint8 **It1,long vi0,long vi1,long vj0,long vj1)
{
	vuint8 l1;
	vuint8 l2;
	vuint8 l3;

	vuint8 result1,result2,result3;

	vuint8 temp;
	vuint8 left,right;


	//premiere ligne : prologue
	int j = vj0;
	int i = vi0;

	l1 = _mm_load_si128(&It[i+0][j]);
	l2 = _mm_load_si128(&It[i+1][j]);

	result1 =  _mm_or_si128(l1,l2);

	left = _mm_slli_si128(result1,1);

	//ici right contient result1 decallé vers la droite avec 255 comme valeur de gauche
	//
	j++;
	for(j;j<=vj1;j++)
	{
		l1 = _mm_load_si128(&It[i+0][j]);
		l2 = _mm_load_si128(&It[i+1][j]);

		result2 = _mm_or_si128(l1,l2);

		right = _mm_srli_si128(result1,1);
		right = _mm_or_si128(right,_mm_slli_si128(result2,15)); // On complete le vecteur gauche en ajoutant a sa doite la valeur de gauvhe du vecteur suivant

		result3 = _mm_or_si128(left,_mm_or_si128(right,result1)); // Valeur finale de IT1[i][j-1]

		left = _mm_slli_si128(result2,1);
		left = _mm_or_si128(left,_mm_srli_si128(result1,15));

		result1 = result2;

		_mm_store_si128(&It1[i][j-1], result3);
	}

	right = _mm_srli_si128(result1,1);

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

		left = _mm_slli_si128(result1,1);

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

			right = _mm_srli_si128(result1,1);
			right = _mm_or_si128(right,_mm_slli_si128(result2,15)); // On complete le vecteur gauche en ajoutant a sa doite la valeur de gauvhe du vecteur suivant

			result3 = _mm_or_si128(left,_mm_or_si128(right,result1));

			left = _mm_slli_si128(result2,1);
			left = _mm_or_si128(left,_mm_srli_si128(result1,15));

			result1 = result2;

			_mm_store_si128(&It1[i][j-1], result3);
		}

		right = _mm_srli_si128(result1,1);

		result3 = _mm_or_si128(left,_mm_or_si128(right,result1)); // Valeur finale de IT1[i][j-1]

		_mm_store_si128(&It1[i][j-1],result3);

		j= vj0;
	}

	l1 = _mm_load_si128(&It[i-1][j]);
	l2 = _mm_load_si128(&It[i+0][j]);

	result1 = _mm_or_si128(l1,l2);

	left = _mm_slli_si128(result1,1);

	//j vaut vj0 donc on insere un zero a gauche de right
	
	j++;
	for(j;j<=vj1;j++)
	{
		l1 = _mm_load_si128(&It[i-1][j]);
		l2 = _mm_load_si128(&It[i+0][j]);

		result2 = _mm_or_si128(l1,l2);

		right = _mm_srli_si128(result1,1);
		right = _mm_or_si128(right,_mm_slli_si128(result2,15)); // On complete le vecteur gauche en ajoutant a sa doite la valeur de gauvhe du vecteur suivant

		result3 = _mm_or_si128(left,_mm_or_si128(right,result1)); // Valeur finale de IT1[i][j-1]

		left = _mm_slli_si128(result2,1);
		left = _mm_or_si128(left,_mm_srli_si128(result1,15));

		result1 = result2;

		_mm_store_si128(&It1[i][j-1], result3);
	}
	//derniere ligne


	right = _mm_srli_si128(result1,1);

	result3 = _mm_or_si128(left,_mm_or_si128(right,result1)); // Valeur finale de IT1[i][j-1]

	_mm_store_si128(&It1[i][j-1],result3);
}

void dilatation3x3_SIMD_B(vuint8 **It,vuint8 **It1,long vi0,long vi1,long vj0,long vj1)
{
	vuint8 l1;
	vuint8 l2;
	vuint8 l3;

	vuint8 result1,result2,result3;

	vuint8 temp;
	vuint8 left,right;


	//premiere ligne : prologue
	int j = vj0;
	int i = vi0;

	
	// corps de boucle
	for(i;i<=vi1;i++)
	{
		l1 = _mm_load_si128(&It[i-1][j]);
		l2 = _mm_load_si128(&It[i+0][j]);
		l3 = _mm_load_si128(&It[i+1][j]);

		temp =  _mm_or_si128(l1,l2);
		result1 = _mm_or_si128(temp,l3);

		left = _mm_slli_si128(result1,1);

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

			right = _mm_srli_si128(result1,1);
			right = _mm_or_si128(right,_mm_slli_si128(result2,15)); // On complete le vecteur gauche en ajoutant a sa doite la valeur de gauvhe du vecteur suivant

			result3 = _mm_or_si128(left,_mm_or_si128(right,result1));

			left = _mm_slli_si128(result2,1);
			left = _mm_or_si128(left,_mm_srli_si128(result1,15));

			result1 = result2;

			_mm_store_si128(&It1[i][j-1], result3);
		}

		right = _mm_srli_si128(result1,1);

		result3 = _mm_or_si128(left,_mm_or_si128(right,result1)); // Valeur finale de IT1[i][j-1]

		_mm_store_si128(&It1[i][j-1],result3);

		j= vj0;
	}
}

void fermeture3x3_SIMD(vuint8 **It,vuint8 **It1,long vi0,long vi1,long vj0,long vj1)
{
	vuint8 ** vXVt = vui8matrix(vi0, vi1, vj0, vj1);
	dilatation3x3_SIMD(It,vXVt,vi0,vi1,vj0,vj1);
	erosion3x3_SIMD(vXVt,It1,vi0,vi1,vj0,vj1);
	free_vui8matrix(vXVt,vi0, vi1, vj0, vj1);
}

void fermeture3x3_SIMD_B(vuint8 **It,vuint8 **It1,long vi0,long vi1,long vj0,long vj1)
{
	vuint8 ** vXVt = vui8matrix(vi0-1, vi1+1, vj0-1, vj1+1);
	dilatation3x3_SIMD_B(It,vXVt,vi0,vi1,vj0,vj1);
	erosion3x3_SIMD_B(vXVt,It1,vi0,vi1,vj0,vj1);
	free_vui8matrix(vXVt,vi0-1, vi1+1, vj0-1, vj1+1);
}

void ouverture3x3_SIMD(vuint8 **It,vuint8 **It1,long vi0,long vi1,long vj0,long vj1)
{
	vuint8 ** vXVt = vui8matrix(vi0, vi1, vj0, vj1);
	erosion3x3_SIMD(It,vXVt,vi0,vi1,vj0,vj1);
	dilatation3x3_SIMD(vXVt,It1,vi0,vi1,vj0,vj1);
	free_vui8matrix(vXVt,vi0, vi1, vj0, vj1);
}

void ouverture3x3_SIMD_B(vuint8 **It,vuint8 **It1,long vi0,long vi1,long vj0,long vj1)
{
	vuint8 ** vXVt = vui8matrix(vi0-1, vi1+1, vj0-1, vj1+1);
	erosion3x3_SIMD_B(It,vXVt,vi0,vi1,vj0,vj1);
	dilatation3x3_SIMD_B(vXVt,It1,vi0,vi1,vj0,vj1);
	free_vui8matrix(vXVt,vi0-1, vi1+1, vj0-1, vj1+1);
}

void erosion5x5_SIMD_B(vuint8 **It,vuint8 **It1,long vi0,long vi1,long vj0,long vj1)
{
	vuint8 ** vXVt = vui8matrix(vi0-1, vi1+1, vj0-1, vj1+1);
	erosion3x3_SIMD_B(It,vXVt,vi0,vi1,vj0,vj1);
	erosion3x3_SIMD_B(vXVt,It1,vi0,vi1,vj0,vj1);
	free_vui8matrix(vXVt,vi0-1, vi1+1, vj0-1, vj1+1);
}

void dilatation5x5_SIMD_B(vuint8 **It,vuint8 **It1,long vi0,long vi1,long vj0,long vj1)
{
	vuint8 ** vXVt = vui8matrix(vi0-1, vi1+1, vj0-1, vj1+1);
	dilatation3x3_SIMD_B(It,vXVt,vi0,vi1,vj0,vj1);
	dilatation3x3_SIMD_B(vXVt,It1,vi0,vi1,vj0,vj1);
	free_vui8matrix(vXVt,vi0-1, vi1+1, vj0-1, vj1+1);
}

void fermeture5x5_SIMD_B(vuint8 **It,vuint8 **It1,long vi0,long vi1,long vj0,long vj1)
{
	vuint8 ** vXVt = vui8matrix(vi0-1, vi1+1, vj0-1, vj1+1);
	dilatation5x5_SIMD_B(It,vXVt,vi0,vi1,vj0,vj1);
	erosion5x5_SIMD_B(vXVt,It1,vi0,vi1,vj0,vj1);
	free_vui8matrix(vXVt,vi0-1, vi1+1, vj0-1, vj1+1);
}

void ouverture5x5_SIMD_B(vuint8 **It,vuint8 **It1,long vi0,long vi1,long vj0,long vj1)
{
	vuint8 ** vXVt = vui8matrix(vi0-1, vi1+1, vj0-1, vj1+1);
	erosion5x5_SIMD_B(It,vXVt,vi0,vi1,vj0,vj1);
	dilatation5x5_SIMD_B(vXVt,It1,vi0,vi1,vj0,vj1);
	free_vui8matrix(vXVt,vi0-1, vi1+1, vj0-1, vj1+1);
}