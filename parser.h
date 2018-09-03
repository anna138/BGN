//
// main.cpp
// Developer: Steve Jankly <spjankly@cpp.edu>
// University: California State Polytechnic University, Pomona
// Project: GenomeCrypt
//

/*
#include "gpu_enc.cuh"
#include "gpu_utilities.cuh"
#include "PaillierCrypt.h"
*/
/*
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
*/
/*
#include <SDKDDKVer.h>
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
*/

double encrypt_genomic_data_sam(std::string filename, int numberSequences, PublicKey * pk);
/*
int main()
{
    //initializeGPU(0);
    //printDeviceProperties();

    std::cout << "Press enter to continue..." << std::endl;
    std::cin.get();

    std::string filename = "HG00239.sam";

    // CPU
    double t_avg_cpu = encrypt_genomic_data_sam(filename, 20, false);
*/  
    /*
    std::cout << "Press enter to continue..." << std::endl;
    std::cin.get();

     GPU
    double t_avg_gpu = encrypt_genomic_data_sam(filename, 20, true);

    std::cout << std::endl;
    double speedup = t_avg_cpu / t_avg_gpu;
    std::cout << "speedup = " << speedup << std::endl;

    std::cout << "Press enter to continue..." << std::endl;
    std::cin.get();
    */
/*
    return 0;
}

double encrypt_genomic_data_sam(std::string filename, int numberSequences, bool useGPU)
{
    
    //LARGE_INTEGER freq, t0, t1; // ticks/s, ticks
    //QueryPerformanceFrequency(&freq);
    //double elapsedMul = 1.0e3 / freq.QuadPart; // 1ms base

    const char nucleobase[] = { 'A', 'C', 'G', 'T' };
    int nbCode[] = { 0, 1, 2, 3 };
    std::ifstream genomeFile;
    genomeFile.open(filename, std::ifstream::in);
    std::string line;
    //PaillierCrypt pc;
    //pc.keygen(PaillierKeyLength_32bit);

    const unsigned int maxSequenceLength = 100;
    //BigUnsigned* data_enc = new BigUnsigned[maxSequenceLength];

    //std::ofstream samencFile;
    //samencFile.open(filename + std::string(".enc"));

    std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(12);
   
    // process tab delimited lines
    std::string item;
    int sequence = 0;
    //BigUnsigned t_total; // in ns
    unsigned int n_total = 0;
    while (std::getline(genomeFile, line) && (sequence < numberSequences)) {
        // ignore header lines
        if (line[0] == '@') {
            // skipping output of header to file
            continue;
        }
        std::stringstream ls(line);
        // go through first 9 items
        for (unsigned int i = 0; i < 9; i++) {
            std::getline(ls, item, '\t');
          //  samencFile << item << '\t';
        }
        // 10th item is data to encrypt
        std::getline(ls, item, '\t');
        unsigned int n = item.length();
        for (unsigned int i = 0; i < n; i++) { // print original
            std::cout << item[i];
        }
        std::cout << std::endl;
        std::cout << ":";
        //samencFile << ":";
        */
        /*
        if (useGPU) { // GPU
            QueryPerformanceCounter(&t0);
            gpu_enc(data_enc, &pc, item.c_str(), n);
            QueryPerformanceCounter(&t1);
        }*/
        /*
        if(true) { // CPU
         //   QueryPerformanceCounter(&t0);
            for (unsigned int i = 0; i < n; i++) {
                char c = item.c_str()[i];
                unsigned int x =
                    c == nucleobase[0] ? nbCode[0] :
                    c == nucleobase[1] ? nbCode[1] :
                    c == nucleobase[2] ? nbCode[2] :
                    c == nucleobase[3] ? nbCode[3] : '?';
                //data_enc[i] = pc.encrypt(x);
            }
            //QueryPerformanceCounter(&t1);
        }
        */
        /*
        for (unsigned int i = 0; i < n; i++) {
            std::cout << data_enc[i] << ":";
            samencFile << data_enc[i] << ":";
        }*/
        /*samencFile << '\t';
        n_total += n;
        t_total += static_cast<unsigned int>((t1.QuadPart - t0.QuadPart) * elapsedMul * 1000000.0); // add ns
        */
        // go through rest of items
      //  while (std::getline(ls, item, '\t')) {
        //    samencFile << item << '\t';
       // }
        /*
        std::cout << std::endl;
    //    samencFile << std::endl;
        sequence++;
    }

    std::cout << std::endl;
    double t_avg = 0;//(t_total / n_total).toUnsignedInt() / 1e9;
    //std::cout << (useGPU ? "GPU" : "CPU") << " test: t_avg: " << t_avg << std::endl;

    genomeFile.close();
   // samencFile.close();
    //delete[] data_enc;

    return t_avg;
}
    */
