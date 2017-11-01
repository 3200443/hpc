#include <stdio.h>
#include <stdlib.h>
#include "test_mouvement.h"
#include "test_mouvement_SSE2.h"


int main(int argc, char* argv[])
{

    /*test_routine_FrameDifference(10);
    test_routine_sigmaDelta();
*/
    test_routine_FrameDifference_SSE2(10);
    return 0;
}
