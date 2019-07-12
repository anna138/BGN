//
// main.cpp
//
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <time.h>
#include <sys/time.h>

#include <gmp.h>
#include <sodium.h>
#include <pbc/pbc.h>

#include "plaintext.h"
#include "bgn.h"
#include "ciphertext.h"
#include "parser.h"

//pbc_param_t params;

void fishersExact(PublicKey *pk, SecretKey * sk, Ciphertext * alleleA,
                  Ciphertext * alleleB, Ciphertext * alleleC,
                  Ciphertext * alleleD, Ciphertext *aB, Ciphertext *aC,
                  Ciphertext *bD, Ciphertext *cD, Ciphertext *total);

//Ciphertext * bayesCompute( );

//bool bayesComparison( );

int main()
{
    clock_t t = 0;
    int keyBits = 32;  // length of q1 and q2
    mpz_t messageSpace;
    int ms[10] = { 1009, 2003, 4001, 8009, 10007, 100003, 1000003, 10000019,
            100000007, 1000000007 };
    mpz_init_set_si(messageSpace, ms[0]);
    int polyBase = 2;  // base for the ciphertext polynomial
    int fpScaleBase = 2;
    double fpPrecision = 0.00000001;
    PublicKey pk;
    SecretKey sk;

    printf("Running tests...\n");
    printf("----------\n");
    printf("NewKeyGen: Test adjusting key bits\n");
    keyBits = 32;  // 32 - 2048
    while(keyBits <= 2048) {
        t = clock();
        NewKeyGen(keyBits, messageSpace, polyBase, fpScaleBase, fpPrecision,
                  false, &pk, &sk);
        t = clock() - t;
        printf("NewKeyGen(keyBits=%d): t = %f\n", keyBits,
               (float) t / CLOCKS_PER_SEC);
        keyBits *= 2;
    }
    printf("----------\n");

    printf("NewKeyGen: Test adjusting message space\n");
    keyBits = 1024;
    for(int i = 0; i < 10; i++) {
        mpz_init_set_si(messageSpace, ms[i]);
        t = clock();
        NewKeyGen(keyBits, messageSpace, polyBase, fpScaleBase, fpPrecision,
                  false, &pk, &sk);
        t = clock() - t;
        printf("NewKeyGen(ms=%d): t = %f\n", ms[i], (float) t / CLOCKS_PER_SEC);
    }
    printf("----------\n");

    printf("NewKeyGen: Test adjusting poly base\n");
    keyBits = 1024;
    mpz_init_set_si(messageSpace, ms[0]);
    polyBase = 2;
    for(int i = 0; i < 10; i++) {
        t = clock();
        NewKeyGen(keyBits, messageSpace, polyBase, fpScaleBase, fpPrecision,
                  false, &pk, &sk);
        t = clock() - t;
        printf("NewKeyGen(polyBase=%d): t = %f\n", polyBase,
               (float) t / CLOCKS_PER_SEC);
        polyBase++;
    }
    printf("----------\n");
  
    //gotta change polyBase to 2 or 3 and generate a key
    //otherwise you are left with a polybase of 11
    //and a polyBase of 11 cannot generate all numbers
    polyBase = 3;
    NewKeyGen(keyBits, messageSpace, polyBase, fpScaleBase, fpPrecision,
              false, &pk, &sk);

    /* Plaintext* pt;
     mpf_t m;
     mpf_init(m);
     mpf_set_d(m,100000000000);
     pt = NewPlaintext(&pk,m);
     Ciphertext *ct = Encrypt(pt,&pk);
     pt = Decrypt(ct,&pk,&sk);
     std::cout << String(pt) << std::endl;
     mpf_set_d(m,0.375);
     pt = NewPlaintext(&pk,m);
     Ciphertext *ct1 = Encrypt(pt,&pk);
     Ciphertext ct2,ct3;
     pt = EMult(ct,ct1,&pk,&ct2, &sk);
     std::cout << String(pt) << std::endl;
     Ciphertext *ct4 = Encrypt(pt,&pk);
     pt = EMult(ct1,ct4,&pk,&ct3, &sk);
     std::cout << String(pt) << std::endl;
     */
    std::string filename = "../Genomic_Data/HG00176.sam";

    printf("encrypt_genomic_data_sam: Running test...\n");
    t = clock();
    double t_avg_cpu = encrypt_genomic_data_sam(filename, 5, &pk);
    printf("%f\n", t_avg_cpu);
    t = clock() - t;
    printf("It took encrypt_genomic_data (%f seconds).\n",
           ((float) t) / CLOCKS_PER_SEC);
    printf("----------\n");

    std::string countallele = "genocount.enc";

    std::ifstream countF;
    countF.open(countallele, std::ifstream::in);
    std::string line;
    std::string line2;
    int count = 0;
    int degree = 0;
    element_t tea;
    Ciphertext allele0;
    Ciphertext allele1;
    Ciphertext allele2;
    Ciphertext allele3;
    element_init_same_as(tea, pk.G1);
    //Will make a function for this code
    //A reader for the Encrypted seq
    while(count < 4) {
        std::getline(countF, line);
        element_set_str(tea, line.c_str(), 10);
        if(line == ":") {
            count++;
            degree = 0;
            continue;
        }
        if(count == 0) {
            element_init_same_as(allele0.Coefficients[degree], tea);
            element_set(allele0.Coefficients[degree], tea);
            allele0.Degree = degree + 1;
        }
        if(count == 1) {
            element_init_same_as(allele1.Coefficients[degree], tea);
            element_set(allele1.Coefficients[degree], tea);
            allele1.Degree = degree + 1;
        }
        if(count == 2) {
            element_init_same_as(allele2.Coefficients[degree], tea);
            element_set(allele2.Coefficients[degree], tea);
            allele2.Degree = degree + 1;
        }
        if(count == 3) {
            element_init_same_as(allele3.Coefficients[degree], tea);
            element_set(allele3.Coefficients[degree], tea);
            allele3.Degree = degree + 1;
        }
        degree++;
    }
    countF.close();
    allele0.ScaleFactor = 0;
    allele0.L2 = false;

    allele1.ScaleFactor = 0;
    allele1.L2 = false;

    allele2.ScaleFactor = 0;
    allele2.L2 = false;

    allele3.ScaleFactor = 0;
    allele3.L2 = false;

    Ciphertext *aB = nullptr, *aC = nullptr, *bD = nullptr, *cD = nullptr,
            *total = nullptr;
    //Plaintext *pt0, *pt1, *pt2, *pt3, *ptotal;
    fishersExact(&pk, &sk, &allele0, &allele1, &allele2, &allele3, aB, aC, bD,
                 cD, total);
    //  */
    /* pt0 = Decrypt(aB, &pk, &sk);
     pt1 = Decrypt(aC, &pk, &sk);
     pt2 = Decrypt(bD, &pk, &sk);
     pt3 = Decrypt(cD, &pk, &sk);
     ptotal = Decrypt(total, &pk, &sk);
     */

    /*Tingting, I needed to put Decrypt in the multiply function in order for it to work.
     Otherwise it would segfault in main. This was my issue I mentioned earlier.
     I tried debugging and cleaning mememory leaks on valgrind.
     I even tried to allocate memory for the Ciphertext on the stack and heap.
     All did not work and I was unable to acces the Ciphertext's Coefficients from EMult to Decrypt.
     The same issue that occurs in the EMult function occurs in Fishers functions
     */

    destructPublicKey(&pk);
    destructSecretKey(&sk);

    mpz_clear(messageSpace);

    /*free(pt);
     //free(pt2);
     //free(pt3);
     */

    cleanTablrus();
    cleanParams();

    return 0;

}

void fishersExact(PublicKey * pk, SecretKey * sk, Ciphertext * alleleA,
                  Ciphertext * alleleB, Ciphertext * alleleC,
                  Ciphertext * alleleD, Ciphertext *aB, Ciphertext *aC,
                  Ciphertext *bD, Ciphertext *cD, Ciphertext *total)
{
    //write a reader for the seq

    clock_t t;

    Plaintext *pt0, *pt1, *pt2, *pt3, *ptotal;

    t = clock();
    aB = EAdd(alleleA, alleleB, pk);
    aC = EAdd(alleleA, alleleC, pk);
    bD = EAdd(alleleB, alleleD, pk);
    cD = EAdd(alleleC, alleleD, pk);
    total = EAdd(aC, bD, pk);
    t = clock() - t;
    printf("It took Fishers(%f seconds).\n", ((float) t) / CLOCKS_PER_SEC);

    t = clock();
    pt0 = Decrypt(aB, pk, sk);
    t = clock() - t;
    printf("It took Decrypt(%f seconds).\n", ((float) t) / CLOCKS_PER_SEC);
    pt1 = Decrypt(aC, pk, sk);
    pt2 = Decrypt(bD, pk, sk);
    pt3 = Decrypt(cD, pk, sk);
    ptotal = Decrypt(total, pk, sk);
    std::cout << String(pt0) << std::endl;
    std::cout << String(pt1) << std::endl;
    std::cout << String(pt2) << std::endl;
    std::cout << String(pt3) << std::endl;
    std::cout << String(ptotal) << std::endl;
    Ciphertext fs;
    Ciphertext * ci = EMult(aC, bD, pk, &fs, sk);

    pt1 = Decrypt(ci, pk, sk);

    std::cout << String(pt1) << std::endl;

}
/*
 Ciphertext * bayesCompute(PublicKey pk, Ciphertext * alleleA, Ciphertext * alleleB, Ciphertext * alleleNB){ //P(A|B)
 Ciphertext * ct;
 ct = EAdd(alleleB, alleleNB, &pk);

 i}

 bool bayesComparison(percentA, percentB){
 ;
 }
 */
