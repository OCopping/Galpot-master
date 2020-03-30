gcc -O3 -I/home/oco/gsl/include -L/home/oco/gsl/lib -c -o orbit_McMillan.o orbit_McMillan.c -lm -lgsl -lgslcblas
g++ -O3 -ffast-math -Isrc/ -c -o foo.o foo.cc  -Lobj -lPot -lOther -lm
g++ -O3 -ffast-math -Isrc/ -I/home/oco/gsl/include -L/home/oco/gsl/lib -o a.out orbit_McMillan.o foo.o -lm -lgsl -lgslcblas -Lobj -lPot -lOther
