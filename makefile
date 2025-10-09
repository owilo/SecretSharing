all:
	g++ secret.cpp image.cpp sis_schemes.cpp -o sis_schemes -O3 -Iinclude