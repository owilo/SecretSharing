.PHONY: single

single:
	g++ src/ss/secret.cpp src/ss/image.cpp src/$(file).cpp -o $(file) -O3 -Iinclude -Ilib/include -std=c++20
