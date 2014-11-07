PROGS = analyzer# nls_test

all: $(PROGS)

analyzer: FE_stress_analysis.cpp
	clear && clear && g++ -g -O2 -o analyzer FE_stress_analysis.cpp -I eigen/ && cp -f analyzer ../

#nls_test: nls_test.cpp
#	clear && clear && g++ -g -o nls_test nls_test.cpp -I /home/mullinix/Dropbox/School/AE510/eigen/

clean:
	clear && clear && rm -f $(PROGS)

