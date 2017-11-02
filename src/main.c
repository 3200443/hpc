#include <stdio.h>
#include <stdlib.h>
#include "test_mouvement.h"
#include "test_mouvement_SSE2.h"

#define OPTI 3 //1 pour optimisation 2 sans optimisation 3 pour tout
int main(int argc, char* argv[])
{
#if OPTI & 0x1
    test_routine_FrameDifference_SSE2(10);
#endif
#if OPTI & 0x2
    test_routine_FrameDifference(10);
    test_routine_sigmaDelta();
#endif
    return 0;
}
