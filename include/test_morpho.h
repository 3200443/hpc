#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vnrutil.h"
#include "nrutil.h"
#include "morpho.h"
#include "morpho_simd.h"


void test_morpho3x3simd();
void test_routine_FrameDifferenceMorpho3x3ouverture(int seuil);
void test_routine_FrameDifferenceMorpho3x3fermeture(int seuil);
void test_routine_FrameDifferenceMorpho3x3ouvertureFermeture(int seuil);
void test_routine_FrameDifferenceMorpho3x3fermetureOuverture(int seuil);
void test_routine_FrameDifferenceMorpho3x3fermeturefermeture(int seuil);
void test_routine_FrameDifferenceMorpho3x3fermeture_pipe(int seuil);
void test_routine_FrameDifferenceMorpho3x3ouverture_pipe(int seuil);
void test_routine_FrameDifferenceMorpho3x3ouverture_bin(int seuil);
void test_routine_FrameDifferenceMorpho3x3fermeture_bin(int seuil);

