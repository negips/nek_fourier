#!/bin/sh

rm logfile1 logfile2

./makenek couette_egv

sed -i 's/numberOfPerturbations=1/numberOfPerturbations=2/g' couette.par
sed -i 's/iff3d=no/iff3d=yes/g' couette.par

nekmpi couette 2 > logfile2

mv tmpcouette0.f00001 tmpcouette0.f00002
mv ptrcouette0.f00001 ptrcouette0.f00002
mv pticouette0.f00001 pticouette0.f00002

sed -i 's/numberOfPerturbations=2/numberOfPerturbations=1/g' couette.par
sed -i 's/iff3d=yes/iff3d=no/g' couette.par

nekmpi couette 2 > logfile1

