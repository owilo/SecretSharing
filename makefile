sis_schemes: secret.cpp image.cpp sis_schemes.cpp
	g++ secret.cpp image.cpp sis_schemes.cpp -o sis_schemes -O3 -Iinclude

shamir_jpeg: secret.cpp image.cpp shamir_jpeg.cpp
	g++ secret.cpp image.cpp shamir_jpeg.cpp -o shamir_jpeg -O3 -Iinclude