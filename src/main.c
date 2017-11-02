#include <stdio.h>
#include <stdlib.h>
#include "test_mouvement.h"
#include "test_mouvement_SSE2.h"

#define OPTI 0

int main(int argc, char* argv[])
{
#if !OPTI
    test_routine_FrameDifference(10);
    test_routine_sigmaDelta();
#endif
#if OPTI
    test_routine_FrameDifference_SSE2(10);

#endif
    return 0;
}
