#!/bin/sh

rm logfile1

./makenek test_arpack

nekmpi couette 2 > logfile1

visnek couette
visnek vbacouette
