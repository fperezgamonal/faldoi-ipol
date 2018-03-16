// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2017, Onofre Martorell <onofremartorelln@gmail.com>
// All rights reserved.

#ifndef FALDOI
#define FALDOI


#include <stdlib.h>
#include "parameters.h"
#include "utils_preprocess.h"
#include <vector>
#include <iostream>
#include <sstream>
//#include <direct.h>

using namespace std;
using namespace std::chrono;

string get_working_path(){
    char temp[MAXPATHLEN];
    return ( getcwd(temp, MAXPATHLEN) ? string( temp ) : string("") );
}

class splitstring : public string {
    vector<string> flds;
public:
    splitstring(char *s) : string(s) { }
    vector<string>& split(char delim, int rep=0);
};

vector<string>& splitstring::split(char delim, int rep) {
    if (!flds.empty()) flds.clear();
    string work = data();
    string buf = "";
    int i = 0;
    while (i < work.length()) {
        if (work[i] != delim)
            buf += work[i];
        else if (rep == 1) {
            flds.push_back(buf);
            buf = "";
        } else if (buf.length() > 0) {
            flds.push_back(buf);
            buf = "";
        }
        i++;
    }
    if (!buf.empty())
        flds.push_back(buf);
    return flds;
}



int main(int argc, char *argv[]) {


    int method = M_TVL1_OCC;

    //M_TVL1       0
    //M_TVL1_W     1
    //M_NLTVL1     2
    //M_NLTVL1_W   3
    //M_TVCSAD     4
    //M_TVCSAD_W   5
    //M_NLTVCSAD   6
    //M_NLTVCSAD_W 7
    //M_TVL1_OCC   8

    bool sparse_fl = True;
    bool local_of = True;
    bool global_of = True;


    system_clock::time_point today = system_clock::now();
    time_t tt;

    tt = system_clock::to_time_t ( today );
    cerr << "today is: " << ctime(&tt);

    // process input
    vector<string> args(argv, argv + argc);


    auto windows_radio  = pick_option(args, "wr", "5"); // Windows radius for local minimization
    auto warps_val      = pick_option(args, "warps", to_string(PAR_DEFAULT_NWARPS_GLOBAL)); // Warpings global
    auto method_ev      = pick_option(args, "vm", to_string(method));

    auto filename_gt        = pick_option(args, "ev",  ""); // Methods
    auto filename_params    = pick_option(args, "p",  ""); // File of parameters
    auto num_file_params    = pick_option(args, "nump",  "");

//    if (args.size() != 1) {
//        fprintf(stderr, "Usage: %lu  ims.txt in_flow.flo  out.flo [-m] val [-w] val \n", args.size());

//        return EXIT_FAILURE;
//    }

    string method_extension = "";

    if (method == M_TVL1){
        method_extension = "tvl1";
    }else{
        if (method == M_TVL1_OCC){
            method_extension = "tvl1_occ";
        }
    }
    string folder_params =  "../Parameters_optimization_files/params_";

    const string& filename_images = args[1];
    ifstream infile(filename_images);
    int num_files = 0;
    string line;
    while(getline(infile, line)){

        ++num_files;

        if (num_files == 1){
            filename_i0  = line;
        }else{
            if (num_files == 2){
                filename_i1  = line;
            }
        }
    }

    infile.close();


    //    splitstring s("Humpty Dumpty sat on a wall.   Humpty Dumpty had a great fall");
    //    cout << s << endl;

    //    vector<string> flds = s.split(' ');
    //    for (int k = 0; k < flds.size(); k++)
    //        cout << k << " => " << flds[k] << endl;

    //    cout << endl << "with repeated delimiters:" << endl;
    //    vector<string> flds2 = s.split(' ', 1);
    //    for (int k = 0; k < flds2.size(); k++)
    //        cout << k << " => " << flds2[k] << endl;

    splitstring s(filename_i0);
    vector<string> output = s.split('.');
    s2(output[2]);
    vector<string> output2 = s.split('/');
    string sequence = output2[output2.size() - 1];
    string core_name1 = output2[output2.size()];


    splitstring s(filename_i1);
    output = s.split('.');
    s2(output[2]);
    output2 = s.split('/');
    string core_name2 = output2[output2.size()];

    //sequence = filename_i0.split('.')[2].split('/')[-2]
    //core_name1 = filename_i0.split('.')[2].split('/')[-1]
    //core_name2 = filename_i1.split('.')[2].split('/')[-1]


    //C++ program names
    string sparse_flow = "../build/sparse_flow";
    string match_propagation = "../build/local_faldoi";
    string of_var = "../build/global_faldoi";
    //std::string evaluation = "../build/evaluation";


    //Set the main directory that contains all the stuff
    root_path = get_working_path();

    string f_path = "../Results/" + sequence + "/";



    system("mkdir " + f_path);

    //if(!(filename_params == "")){
    //        iteration_params = "_" + str(filename_params.split('/')[-1].split('.')[0].split('_')[1]).zfill(2)
    //}else{
    //        iteration_params = ""
    //}


    //#Set the images input.
    int width_im;
    int height_im;
    //im_name0 = os.path.abspath(data[0])
    //im_name1 = os.path.abspath(data[1])
    //#Get the image size
    //with open(im_name1, 'r') as f:
    //    image = Image.open(f)
    //    width_im = image.size[0]
    //    height_im = image.size[1]

    //os.chdir(binary_path)

    string match_name_1 = f_path + core_name1 + "_exp_mt_1.txt";
    string sparse_name_1 = f_path + core_name1 + "_exp_mt_1.flo";
    string sparse_in_1 = f_path + core_name1 + "_exp_mt_1_saliency_out_cut.flo";


    string match_name_2 = f_path + core_name2 + "_exp_mt_2.txt";
    string sparse_name_2 = f_path + core_name2 + "_exp_mt_2.flo";
    string sparse_in_2 = f_path + core_name2 + "_exp_mt_2_saliency_out_cut.flo";

    string region_growing = f_path + core_name1 + "_rg_" + method_extension + iteration_params + ".flo";
    string sim_value = f_path + core_name1 + "_exp_sim_" + method_extension + iteration_params + ".tiff";
    string var_flow = f_path + core_name1 + "_exp_var_" + method_extension + iteration_params + ".flo";

    string occlusions_rg = f_path + core_name1 + "_rg_occ_" + method_extension + iteration_params + ".png";
    string occlusions_var = f_path + core_name1 + "_var_occ_" + method_extension + iteration_params + ".png" ;


    //#Create a sparse flow from the deepmatching matches.
    //print('Creating sparse from matches')
    if (sparse_fl){
        string param = sparse_in_1 + width_im + height_im + sparse_name_1;
        string command_line = sparse_flow + param + "\n";
        system(command_line);
    }



    if (sparse_fl){
        string param = sparse_in_2 + width_im + height_im + sparse_name_2;
        string command_line = sparse_flow + param + "\n";
        system(command_line);
    }

    //#Create a dense flow from a sparse set of initial seeds


    //#print(command_line)
    if (local_of){
        string options = "";
        cout << "Computing local faldoi";
        if (!(filename_params == "")){
            options = "-m " + method_ev + " -wr " + windows_radio + " -p " + filename_params;
        }else{
            options = "-m " + method_ev + " -wr " + windows_radio;
        }
        string param = args.file_images + sparse_name_1 + sparse_name_2 + region_growing + sim_value + occlusions_rg + options + "\n";
        string command_line = match_propagation + param + "\n";
        system(command_line);
    }

    //#Put the dense flow as input for a variational method

    //# Tv-l2 coupled 0 Du 1

    //#print(command_line)
    if( global_of){
        string options = "";
        cout << "Computing global faldoi";
        if (!(filename_params == "")){
             options = "-m " + method_ev + "-w " + warps_val + " -p " + filename_params;
        }else{
            options = "-m " + method_ev + "-w " + warps_val;
        }
        string param = args.file_images + region_growing + var_flow + occlusions_rg + occlusions_var + options;
        string command_line = of_var + param + "\n";
        system(command_line);
    }
    //#Evaluate results of method
    //    if  not filename_gt == '':
    //        command_line = evaluation + ' ' + var_flow + ' ' + filename_gt
    //    #print(command_line)
    //        os.system(command_line)
    //      today = system_clock::now();

    tt = system_clock::to_time_t ( today );
    cerr << "today is: " << ctime(&tt);
    return EXIT_SUCCESS;
}

#endif
