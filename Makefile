distance_oracle.dylib: distance_oracle.c distance_oracle.h Makefile
	clang -DNDEBUG -dynamiclib --std=c99 -Ofast -Weverything -Wno-c99-extensions -Wno-declaration-after-statement -o $@ $<
