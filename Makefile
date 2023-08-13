distance_oracle.dylib: distance_oracle.cc distance_oracle.h Makefile
	clang++ -DNDEBUG -dynamiclib --std=c++17 -Ofast -Weverything -Wno-c++98-compat -o $@ $<
