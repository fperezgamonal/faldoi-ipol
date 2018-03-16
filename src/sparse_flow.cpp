#include <cstdlib>
#include <cassert>
#include <string>
#include <fstream>
#include <limits>
#include <cmath>
#include <iostream>
#include <cstdio>
extern "C" {
#include "iio.h"
}

static int sparse_optical_flow(char *input, int nx, int ny, float *out) {

    float x1, x2, y1, y2;
    std::string filename_sift_matches(input);
    std::ifstream file(input);
    std::string str;

    //Initialize all the the optical flow to NAN
    for (int j = 0; j < ny; j++){
        for (int i = 0; i < nx; i++){
            out[j*nx + i] = NAN;
            out[nx*ny + j*nx + i] = NAN;
        }
    }
    if (file){
        //Insert the sparse flow obtained from matches
        while (getline(file, str)){
            //Colum I0, Row I0, Colum I1, Row I1
            sscanf(str.c_str(), "%f %f %f %f\n",
                   &x1, &y1, &x2, &y2);
            float u = x2 - x1;
            float v = y2 - y1;
            int i = std::floor(x1);
            int j = std::floor(y1);
            //fprintf(stderr, "Colum x row: %d x %d u x v: (%f,%f)\n", i,j,u,v);
            out[j*nx + i] = u;
            out[nx*ny + j*nx + i] = v;
        }
        return 1;
    }else{
        std::cout << "File does not exist\n";
        std::cout << input << "\n";
        return 0;
    }
}


int main(int argc, char *argv[]) {
    // process input arguments
    if (argc != 5) {
        fprintf(stderr, "usage:\n\t%s sift_matches.txt colum row out.flo\n", *argv);
        //                          0 				1      2   3    4
        fprintf(stderr, "usage:\n\t Nargs:%d\n", argc);
        return 1;
    }

    //Read the input parameter
    char *filename_in = argv[1];
    char *filename_out = argv[4];
    int   nx = atoi(argv[2]);
    int   ny = atoi(argv[3]);
    float *out = new float[2*nx*ny];
    //Call the function that creates the sparse optical flow
    sparse_optical_flow(filename_in, nx, ny, out);
    //Save the optical flow
    iio_save_image_float_split(filename_out, out, nx, ny, 2);

    delete [] out;
    return 0;
}

