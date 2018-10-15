#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <iostream>

#include "plaintext.h"
#include "bgn.h"
#include "ciphertext.h"

#include <gmp.h>
#include <sodium.h>
#include <pbc/pbc.h>

//pbc_param_t params;

void fishersExact(int64_t alleleA, int64_t allele2, int64_t allele3, int64_t allele4, double &ab, double &ac, double &bd, double &cd);

Ciphertext * bayesCompute( );

bool bayesComparison( );

int main(){
        int keyBits = 512; // length of q1 and q2
        mpz_t messageSpace;
        mpz_init_set_si(messageSpace, 1021);
        int polyBase = 2; // base for the ciphertext polynomial
        int fpScaleBase = 2;
        double fpPrecision = 0.0001;
        PublicKey pk;
        SecretKey sk;
        NewKeyGen(keyBits, messageSpace, polyBase, fpScaleBase, fpPrecision, true, &pk, &sk);        
        Ciphertext * ct;
        Ciphertext * ct2; 
        //Ciphertext * ct3 = (Ciphertext *)malloc(sizeof(Ciphertext));
        Ciphertext ct3;
        Ciphertext ct4;
        //Ciphertext * ct4 = (Ciphertext *)malloc(sizeof(Ciphertext));
        Plaintext * pt;
        Plaintext * pt2;
        Plaintext * pt3;
        mpf_t message;
        double var1 = 123.0, var2 = 123.0;
        //printf("This System uses Polybase: %d \n", polyBase); 
        mpf_init_set_d(message, var1);
        pt = NewPlaintext(&pk, message);
      
        //gmp_printf("Main: Encryption function\n");
        ct = Encrypt(pt, &pk);
        //gmp_printf("Main: Decryption function\n"); 
        //pt2 = Decrypt(ct, &pk, &sk);
        
        //printf("Testing Multiplication:\n");
        //int i;
        //for(i = 0; i < 120; i++){
        mpf_init_set_d(message, var2);
        pt = NewPlaintext(&pk, message);
        printf("Main: Testing Addition Function\n");
        ct2 = Encrypt(pt, &pk); 
        ct = EAdd(ct, ct2, &pk);
        pt3 = Decrypt(ct, &pk, &sk);

        /*
        printf("Multiply Function: \n");
        pt3 = EMult(ct, ct2, &pk, &ct3, &sk); 
        */
        /*Tingting, I needed to put Decrypt in the multiply function in order for it to work. 
        Otherwise it would segfault in main. This was my issue I mentioned earlier. 
        I tried debugging and cleaning mememory leaks on valgrind. 
        I even tried to allocate memory for the Ciphertext on the stack and heap. 
        All did not work and I was unable to acces the Ciphertext's Coefficients from EMult to Decrypt.
        */
        
        std::cout << std::endl << "Adding "<< var1 <<" and " << var2  << std::endl;
        std::cout << "The result is: "<< String(pt3) << std::endl << std::endl;
        //printf("%s\n\n",String(pt3));
        //pt3 = EMult(&ct3, ct2, &pk, &ct4, &sk);
        //printf("%s\n\n",String(pt3));
        //EMult(ct, ct2, &pk, &ct3, &sk);
        //EMult(ct, ct2, &pk, ct4);
        printf("SAME LINE Main this Degree L54: %d \n", ct->Degree);
        element_printf("Main this L56: %B \n",ct->Coefficients[0]);
        printf("SAME LINE Main this Degree L56: %d \n", ct->Degree); 
        element_printf("Main this L56: %B \n",ct->Coefficients[0]);
        std::cout << std::endl << "Adding "<< var1 <<" and " << var2  << std::endl;
        std::cout << "The result is: "<< String(pt3) << std::endl << std::endl;
        /*
        printf("SAME LINE CT3 Main this Degree L59: %d \n", ct->Degree);
        element_printf("Main this L56: %B \n",ct->Coefficients[0]);
        //printf("SAME LINE Main this Degree L56: %d \n", ct3.Degree);
        //element_printf("Main this L56: %B \n",ct3.Coefficients[0]);
        //pt3 = Decrypt(&ct3, &pk, &sk);       
        printf("Adding %f to %f is: ", var1, var2);
        printf("%s\n\n",String(pt3));
        //}*/
        /*
        element_clear(pk.G1);
        element_clear(pk.P);
        element_clear(pk.Q);
        */
        destructPublicKey(&pk);
        destructSecretKey(&sk);
        
        destructCiphertext(ct);
        destructCiphertext(ct2);
        //destructCiphertext(&ct3);
        mpz_clear(messageSpace);
        
        mpf_clear(message);
        free(pt);
        //free(pt2);
        //free(pt3);
        

        cleanTablrus();
        cleanParams();
        return 0;

}

void fishersExact(int64_t alleleA, int64_t alleleB, int64_t alleleNA, int64_t alleleNB, double &aB, double &aNB, double &bA, double &bNA){
        //write a reader for the seq
        ;
}


Ciphertext * bayesCompute( ){
        ;
}

bool bayesComparison( ){
        ;
}
