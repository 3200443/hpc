#include <stdio.h>
#include <stdlib.h>
#include "test_mouvement.h"
#include "test_mouvement_SSE2.h"

#define OPTI 1
int main(int argc, char* argv[])
{
#if OPTI == 1
    test_routine_FrameDifference_SSE2(10);
#else
    test_routine_FrameDifference(10);
    test_routine_sigmaDelta();
#endif
    return 0;
}
