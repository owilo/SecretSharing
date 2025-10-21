sis_schemes: src/secret.cpp src/image.cpp src/sis_schemes.cpp
	g++ src/secret.cpp src/image.cpp src/sis_schemes.cpp -o sis_schemes -O3 -Iinclude -Ilib/include -std=c++20

shamir_jpeg: src/secret.cpp src/image.cpp src/shamir_jpeg.cpp
	g++ src/secret.cpp src/image.cpp src/shamir_jpeg.cpp -o shamir_jpeg -O3 -Iinclude -Ilib/include -std=c++20

shamir_jpeg_symmetric: src/secret.cpp src/image.cpp src/shamir_jpeg_symmetric.cpp
	g++ src/secret.cpp src/image.cpp src/shamir_jpeg_symmetric.cpp -o shamir_jpeg_symmetric -O3 -Iinclude -Ilib/include -std=c++20

shamir_jpeg_x: src/secret.cpp src/image.cpp src/shamir_jpeg_x.cpp
	g++ src/secret.cpp src/image.cpp src/shamir_jpeg_x.cpp -o shamir_jpeg_x -O3 -Iinclude -Ilib/include -std=c++20

shamir_jpeg_median: src/secret.cpp src/image.cpp src/shamir_jpeg_median.cpp
	g++ src/secret.cpp src/image.cpp src/shamir_jpeg_median.cpp -o shamir_jpeg_median -O3 -Iinclude -Ilib/include -std=c++20

jpeg_quality_factor: src/secret.cpp src/image.cpp src/jpeg_quality_factor.cpp
	g++ src/secret.cpp src/image.cpp src/jpeg_quality_factor.cpp -o jpeg_quality_factor -O3 -Iinclude -Ilib/include -std=c++20