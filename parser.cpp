//
// parser.cpp
// Developer: Steve Jankly <spjankly@cpp.edu>
// University: California State Polytechnic University, Pomona
//
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include "ciphertext.h"
#include "plaintext.h"
#include "bgn.h"
#include "parser.h"
#include <cstdio>

#include <gmp.h>
#include <sodium.h>
#include <pbc/pbc.h>

//#include <SDKDDKVer.h>
//#define WIN32_LEAN_AND_MEAN
//#include <windows.h>
/*
 double encrypt_genomic_data_sam(std::string filename, int numberSequences, bool useGPU);

 int main()
 {
 //initializeGPU(0);
 //printDeviceProperties();

 std::cout << "Press enter to continue..." << std::endl;
 std::cin.get();

 std::string filename = "HG00239.sam";

 // CPU
 double t_avg_cpu = encrypt_genomic_data_sam(filename, 20, false);

 *//*
 std::cout << "Press enter to continue..." << std::endl;
 std::cin.get();

 GPU
 double t_avg_gpu = encrypt_genomic_data_sam(filename, 20, true);

 std::cout << std::endl;
 double speedup = t_avg_cpu / t_avg_gpu;
 std::cout << "speedup = " << speedup << std::endl;

 std::cout << "Press enter to continue..." << std::endl;
 std::cin.get();

 return 0;
 }*/

double encrypt_genomic_data_sam(std::string filename, int numberSequences,
                                PublicKey * pk)
{

    //LARGE_INTEGER freq, t0, t1; // ticks/s, ticks
    //QueryPerformanceFrequency(&freq);
    //double elapsedMul = 1.0e3 / freq.QuadPart; // 1ms base

    const char nucleobase[] = { 'A', 'C', 'G', 'T' };
    int nbCode[] = { 0, 1, 2, 3 };
    std::ifstream genomeFile;
    genomeFile.open(filename, std::ifstream::in);
    std::string line;

    //bgn.keygen(PaillierKeyLength_32bit);

    const unsigned int maxSequenceLength = 200;
    Ciphertext** data_enc = new Ciphertext*[maxSequenceLength];
    double countA = 0.0, countC = 0.0, countG = 0.0, countT = 0.0;
    FILE *samencFile;
    samencFile = fopen("genomic.enc", "w");
    FILE *countFile;
    countFile = fopen("genocount.enc", "w");
    Ciphertext** count_enc = new Ciphertext*[4];
    //Need to add this in instead of ofstream because the 
    //element_out_str() in pbc library will read file pointers not ofstream

    std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(12);

    // process tab delimited lines
    std::string item;
    int sequence = 0;
    //BigUnsigned t_total; // in ns
    //   unsigned int n_total = 0;
    while(std::getline(genomeFile, line) && (sequence < numberSequences)) {
        // ignore header lines
        if(line[0] == '@') {
            // skipping output of header to file
            continue;
        }
        std::stringstream ls(line);
        // go through first 9 items
        for(unsigned int i = 0; i < 9; i++) {
            std::getline(ls, item, '\t');
            for(unsigned int j = 0; j < item.length(); j++) {
                fprintf(samencFile, " %c \t", item[j]);
            }
        }
        // 10th item is data to encrypt
        std::getline(ls, item, '\t');
        unsigned int n = item.length();
        for(unsigned int i = 0; i < n; i++) {  // print original
            std::cout << item[i];
        }
        std::cout << std::endl;
        std::cout << ":";
        fprintf(samencFile, ":");

        /* if (useGPU) { // GPU
         QueryPerformanceCounter(&t0);
         gpu_enc(data_enc, &pc, item.c_str(), n);
         QueryPerformanceCounter(&t1);
         }*/
        mpf_t mf;
        mpf_init(mf);
        if(true) {  // CPU
            //   QueryPerformanceCounter(&t0);
            for(unsigned int i = 0; i < n; i++) {
                char c = item.c_str()[i];
                unsigned int x = c == nucleobase[0] ? nbCode[0] :
                                 c == nucleobase[1] ? nbCode[1] :
                                 c == nucleobase[2] ? nbCode[2] :
                                 c == nucleobase[3] ? nbCode[3] : '?';
                mpf_set_ui(mf, x);
                if(x == 0) {
                    countA++;
                } else if(x == 1) {
                    countC++;
                } else if(x == 2) {
                    countG++;
                } else if(x == 3) {
                    countT++;
                }
                //Encrypts all ACGTs
                data_enc[i] = Encrypt(NewPlaintext(pk, mf), pk);
            }
            //QueryPerformanceCounter(&t1);
        }
        //inserts all Encrypted ACGTs into samencFile
        for(unsigned int i = 0; i < n; i++) {
            for(int j = 0; j < data_enc[i]->Degree; j++) {
                //element_printf(" %B:",  data_enc[i]->Coefficients[j]);
                element_fprintf(samencFile, "%B:",
                                data_enc[i]->Coefficients[j]);

            }
        }
        fprintf(samencFile, "\t");

        //n_total += n;
        //t_total += static_cast<unsigned int>((t1.QuadPart - t0.QuadPart) * elapsedMul * 1000000.0); // add ns

        // go through rest of items
        while(std::getline(ls, item, '\t')) {
            fprintf(samencFile, " %s \t", item.c_str());
        }
        // std::cout << std::endl;
        fprintf(samencFile, "\n");
        sequence++;
    }
    //Encrypt total count of A,C,G and T
    mpf_t mf;
    mpf_init(mf);

    mpf_set_d(mf, countA);
    count_enc[0] = Encrypt(NewPlaintext(pk, mf), pk);
    mpf_set_d(mf, countC);
    count_enc[1] = Encrypt(NewPlaintext(pk, mf), pk);
    mpf_set_d(mf, countG);
    count_enc[2] = Encrypt(NewPlaintext(pk, mf), pk);
    mpf_set_d(mf, countT);
    count_enc[3] = Encrypt(NewPlaintext(pk, mf), pk);

    //QueryPerformanceCounter(&t1);

    //inserts all counts into countFile
    for(unsigned int i = 0; i < 4; i++) {
        for(int j = 0; j < count_enc[i]->Degree; j++) {
            element_fprintf(countFile, "%B\n", count_enc[i]->Coefficients[j]);
        }
        fprintf(countFile, ":\n");

    }

    std::cout << std::endl;
    double t_avg = 0;        //(t_total / n_total).toUnsignedInt() / 1e9;
    //std::cout << (useGPU ? "GPU" : "CPU") << " test: t_avg: " << t_avg << std::endl;

    genomeFile.close();
    fclose(samencFile);
    fclose(countFile);
    delete[] data_enc;
    delete[] count_enc;
    return t_avg;
}
