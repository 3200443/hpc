rm -r car3F*
rm -r car3S*
rm -r obj

mkdir obj

mkdir car3Frame
mkdir car3Frame3x3F
mkdir car3Frame3x3FF
mkdir car3Frame3x3FO
mkdir car3Frame3x3O
mkdir car3Frame3x3OF
mkdir car3FrameSIMD
mkdir car3FrameSIMD_M
mkdir car3Sigma
mkdir car3SigmaSIMD
mkdir car3Frame3x3O_pipe
mkdir car3Frame3x3F_pipe
mkdir car3Frame3x3F_bin
mkdir car3Frame3x3O_bin
make clean && make
