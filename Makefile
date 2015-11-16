all:		OQCompress
OQCompress:	OQCompress.cc BGZF.h BGZF.cc
	g++ -std=c++11 -fno-strict-aliasing -Wextra -Wall -Wsign-promo -Woverloaded-virtual -Wendif-labels -march=native -g -o OQCompress OQCompress.cc BGZF.cc -lz
