sis_schemes: secret.cpp image.cpp sis_schemes.cpp
	g++ secret.cpp image.cpp sis_schemes.cpp -o sis_schemes -O3 -Iinclude -std=c++20

shamir_jpeg: secret.cpp image.cpp shamir_jpeg.cpp
	g++ secret.cpp image.cpp shamir_jpeg.cpp -o shamir_jpeg -O3 -Iinclude -std=c++20

shamir_jpeg_symmetric: secret.cpp image.cpp shamir_jpeg_symmetric.cpp
	g++ secret.cpp image.cpp shamir_jpeg_symmetric.cpp -o shamir_jpeg_symmetric -O3 -Iinclude -std=c++20

shamir_jpeg_x: secret.cpp image.cpp shamir_jpeg_x.cpp
	g++ secret.cpp image.cpp shamir_jpeg_x.cpp -o shamir_jpeg_x -O3 -Iinclude -std=c++20