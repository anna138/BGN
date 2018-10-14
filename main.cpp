#define PBC_DEBUG
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <gmp.h>
#include <pbc.h>
#include <sodium.h>
#include <math.h>
#include <cstdlib>
#include <iostream>
//pbc_param_t params;
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include "plaintext.h"
#include "bgn.h"
#include "ciphertext.h"
#include "parser.h"

void fishersExact(PublicKey *pk,SecretKey * sk, Ciphertext * alleleA, Ciphertext * alleleB, Ciphertext * alleleC, Ciphertext * alleleD, Ciphertext *aB, Ciphertext *aC, Ciphertext *bD, Ciphertext *cD, Ciphertext *total);

//Ciphertext * bayesCompute( );

//bool bayesComparison( );

int main(){
        std::cout << "Press enter to continue..." << std::endl;
        std::cin.get();

        std::string filename = "HG00176.sam";
        
        int keyBits = 256; // length of q1 and q2
        mpz_t messageSpace;
        mpz_init_set_si(messageSpace, 1201);//1021
	
        int polyBase = 2; // base for the ciphertext polynomial
        int fpScaleBase = 2;
        double fpPrecision = 0.00000001;
        PublicKey pk;
        SecretKey sk;
        clock_t t;
        t = clock();
        NewKeyGen(keyBits, messageSpace, polyBase, fpScaleBase, fpPrecision, false, &pk, &sk);        
        t = clock() - t;
        printf ("It took NewKeyGen (%f seconds).\n",((float)t)/CLOCKS_PER_SEC);
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
        t = clock();
        double t_avg_cpu = encrypt_genomic_data_sam(filename,5, &pk);
        t = clock() - t;
        printf ("It took encrypt_genomic_data (%f seconds).\n",((float)t)/CLOCKS_PER_SEC);
        
        
        std::string countallele = "genocount.enc"; 
        
        std::ifstream countF;
        countF.open(countallele, std::ifstream::in);
        std::string line;
        std::string line2;
        int count = 0;
        int degree = 0; 
        element_t tea;
        Ciphertext allele0; Ciphertext allele1; Ciphertext allele2; Ciphertext allele3;
        element_init_same_as(tea, pk.G1);
        //Will make a function for this code  
        //A reader for the Encrypted seq
        while(count < 4){ 
        std::getline(countF, line);
        element_set_str(tea, line.c_str(), 10);
                if(line == ":"){
                        count++;
                        degree = 0; 
                        continue;
                }
                if(count == 0){
                        element_init_same_as(allele0.Coefficients[degree], tea);
                        element_set(allele0.Coefficients[degree], tea);
                        allele0.Degree = degree + 1;
                }
                if(count == 1){
                        element_init_same_as(allele1.Coefficients[degree], tea); 
                        element_set(allele1.Coefficients[degree], tea);
                        allele1.Degree = degree + 1;
                }
                if(count == 2){
                        element_init_same_as(allele2.Coefficients[degree], tea); 
                        element_set(allele2.Coefficients[degree], tea);
                        allele2.Degree = degree + 1;
                }
                if(count == 3){
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
        
        Ciphertext *aB, *aC, *bD, *cD, *total;
        Plaintext *pt0, *pt1, *pt2, *pt3, *ptotal;
        fishersExact(&pk,&sk, &allele0, &allele1, &allele2, &allele3, aB, aC, bD, cD, total);
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

void fishersExact(PublicKey * pk,SecretKey * sk ,Ciphertext * alleleA, Ciphertext * alleleB, Ciphertext * alleleC, Ciphertext * alleleD, Ciphertext *aB, Ciphertext *aC, Ciphertext *bD, Ciphertext *cD, Ciphertext *total){
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
        printf ("It took Fishers(%f seconds).\n",((float)t)/CLOCKS_PER_SEC);
        
        
        t = clock();
        pt0 = Decrypt(aB, pk, sk);
        t = clock() - t;
        printf ("It took Decrypt(%f seconds).\n",((float)t)/CLOCKS_PER_SEC);
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
        Ciphertext * ci = EMult(aC,bD, pk,&fs,sk);
	printf("here\n");	
	element_printf("[Element: %d] %B \n\n\n", 0, ci->Coefficients[1]);                
//	element_printf("[Element: %d] %B \n\n\n", 0, ci->Coefficients[1]);                
//	element_printf("[Element: %d] %B \n\n\n", 0, ci->Coefficients[1]);                
        for(int i = 0; i < ci->Degree; i++)
		element_printf("element[%d]: %B \n", i, ci->Coefficients[i]);
        
	element_t co;
//	element_init_same_as(co,ci->Coefficients[0]);
        for(int i = 0; i < ci->Degree; i++)
		element_printf("element[%d]: %B \n", i, ci->Coefficients[i]);
	//ci = Copy(&fs);
	pt1 = Decrypt(ci,pk,sk);
        
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
