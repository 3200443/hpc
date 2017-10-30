#include <stdio.h>
#include <stdlib.h>
#include "test_mouvement.h"

//#include "test_mouvement_SSE2.h"


int main(int argc, char* argv[])
{

    test_routine_FrameDifference(atoi(argv[1]));
    test_routine_sigmaDelta();

    return 0;
}
