/*This is BGN.c*/
#include <inttypes.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <gmp.h>
#include <pbc.h>
#include <sodium.h>
#include <math.h>
#include <cstdlib>
#include <iostream>
pbc_param_t params;
#include "bgn.h"
#include "plaintext.h"

void constructPublicKey(PublicKey *that, pairing_t* par, element_t g1, element_t p1, element_t q1, mpz_t n1, mpz_t t1, int PB, int FPS, double FPP, bool DTRM){
        that->Pairing = par;
        element_init_same_as(that->G1, g1);
        element_set(that->G1, g1);
        element_init_same_as(that->P, p1);
        element_set(that->P, p1);
        element_init_same_as(that->Q, q1);
        element_set(that->Q, q1);

        mpz_init_set(that->N, n1);
        mpz_init_set(that->T, t1);
        that->PolyBase = PB;
        that->FPScaleBase = FPS;
        that->FPPrecision = FPP;
        that->Deterministic = DTRM; 
}
void destructPublicKey(PublicKey *that){

        element_clear(that->G1);
        element_clear(that->P);
        element_clear(that->Q);
        
        mpz_clear(that->N);
        mpz_clear(that->T);
}
void cleanParams(){
        pbc_param_clear(params);
}
void constructSecretKey(SecretKey *that, mpz_t s1, int PB){
        mpz_init_set(that->Key, s1);
        that->PolyBase = PB;
}
void destructSecretKey(SecretKey *that){
        mpz_clear(that->Key);
}
void constructPlaintext(Plaintext *that, PublicKey *pk, int64_t *coeff, int degree, int sfactor){
        that->Pk = pk;
        int i;
        for(i = 0; i < 128; i++){
                that->Coefficients[i] = 0;
        }

        for(i = 0; i < degree; i++){
                that->Coefficients[i] = coeff[i];
        //        printf("% " PRId64, that->Coefficients[i]);
        }
//        printf("\n\n");
        that->Degree = degree; 
        that->ScaleFactor = sfactor;
}
double max(double x, double y){
        return(x > y) ? x : y;
}
void NewKeyGen(int keyBits, mpz_t T, int polyBase, int fpScaleBase, double fpPrecision, bool deterministic, PublicKey * pk, SecretKey * sk){
        if(keyBits < 16){
                perror("key bits must be >= 16 bits in length");
        }
        
        mpz_t q1;  /*first random prime for product of N*/ 
        mpz_t q2;  /*second random prime (secret key) for product N*/
        mpz_t N;   /* N = r*q */
        element_t P;     /*field element*/
        element_t Q;     /*field element*/
        mpz_init(q1);
        mpz_init(q2);
        /*generate a new random prime r*/
        if(sodium_init() < 0){
                perror("Sodium Lib not initialized");
        }
      
        uint32_t randseed1, randseed2; /*Random numbers to get the q1 and q2, 
                                 sodium only has uint32*/
        randseed1 = randombytes_random();/*this provides a more 
                                           random state than time(NULL) */
        randseed2 = randombytes_random();
        
        /*This process is to get a random number and pass that number into a algorithm 
          for selecting. This algorithm in function randstate uses randomness like 
          Mersenne Twister with the state. the randstate has both the algorithm */

        gmp_randstate_t gmprand;/*this will have an algorithm selection and 
                                  state. you can use the gmp_randseed*/
        gmp_randinit_mt(gmprand);/*this will initialize state with Mersenne Twister*/
        gmp_randseed_ui(gmprand, randseed1);/*Set seed to the state for randomness*/
        mpz_urandomb(q1,gmprand,keyBits); /*gives you mpz from the selected algo)*/
        
        do{
                mpz_nextprime(q1, q1); /*Take the q1 and put it in q1*/
        }while(mpz_probab_prime_p(q1,50) < 1); 
        gmp_randseed_ui(gmprand,randseed2);
        mpz_urandomb(q2,gmprand,keyBits);
        
        do{
                mpz_nextprime(q2, q2); /*Take the q2 and put it in q2*/
        }while(mpz_probab_prime_p(q2,50) < 1);
        
        /*this is for the error if the message is greater*/
        if(mpz_cmp(q1, T) < 0 || mpz_cmp(q2, T) < 0){
                perror("The message space is greater than the group order");
                exit(1);
        }
        
        /*compute product for N*/
        mpz_init(N);
        mpz_mul(N,q1,q2);
        mpz_mul(N, T, N); //Check with line N.Mul(N,T)
        pbc_param_init_a1_gen(params, N);
        
        pairing_t pairing;
        pairing_init_pbc_param(pairing, params);
        element_t G1;
        element_init_G1(G1, pairing);
        mpz_t l;
        mpz_init(l);
        mpz_t * ptr = parseLFromPBCParams(params); 
        mpz_set(l, *ptr);
        /*obtain 1 generated from the pbc library 
        a is a small number*/
        element_random(G1);
        element_init_same_as(P,G1);
        element_set(P,G1);
        mpz_t j;
        mpz_init_set_ui(j, 4);
        mpz_mul(l,l,j);
        element_pow_mpz(P,P,l);
        mpz_t * ptr2 = newCryptoRandom(N);
        element_pow_mpz(P,P,q2);
        element_init_same_as(Q,G1);
        element_set(Q,G1);
        element_pow_mpz(P,P,*ptr2);
        element_pow_mpz(Q, Q, q2);
        element_pow_mpz(Q, Q, T);
        constructPublicKey(pk, &pairing, G1, P, Q, N, T, polyBase, fpScaleBase, fpPrecision,deterministic);
        constructSecretKey(sk, q1, polyBase);

        mpz_clear(l);
        mpz_clear(q1);
        mpz_clear(q2);
        mpz_clear(*ptr);
        free(ptr);
        mpz_clear(*ptr2);
        free(ptr2);
        element_clear(G1);
        element_clear(Q);
        mpz_clear(N);
        element_clear(P);
        mpz_clear(j);
        gmp_randclear(gmprand);
        
        // add the Public and Secret Key matching

        computeEncodingTable(pk);
}

/*Encrypt a given plaintext (integer or rational) polynomial with the public key pk)*/
Ciphertext * Encrypt(Plaintext *pt, PublicKey * pk){
          element_t * encryptedCoefficients = (element_t *)(malloc(pt->Degree*(sizeof(element_t))));
          
          int i;
          mpz_t coeff; /*struct Plaintext pt;*/
          element_t * temp;
          for(i = 0; i < pt->Degree; i++){
                if(pt->Coefficients[i] < 0){
                        mpz_init_set_d(coeff, (double)(-1 * pt->Coefficients[i]));
                        temp = ESubElements(*encryptZero(pk), *EncryptElement(coeff,pk),pk);
                        element_init_same_as(encryptedCoefficients[i], *temp);
                        element_set(encryptedCoefficients[i], *temp);
                }
                else{
                        mpz_init_set_d(coeff, (double)(pt->Coefficients[i]));
                        temp = EncryptElement(coeff,pk);
                        element_init_same_as(encryptedCoefficients[i], *temp);
                        element_set(encryptedCoefficients[i],*temp);
                }
          }
          Ciphertext * ct = (Ciphertext *)malloc(1*sizeof(Ciphertext));
          constructCiphertext(ct,encryptedCoefficients, pt->Degree, pt->ScaleFactor, false); 
          
          mpz_clear(coeff);
          for(i = 0; i < pt->Degree; i++){
                element_clear(encryptedCoefficients[i]);
          }
          return ct; 
}

Ciphertext * AInv (Ciphertext *ct, PublicKey * pk){
        if(ct->L2){
                return aInvL2(ct,pk);
        }
        int degree;
        degree = ct->Degree;
        element_t * result = (element_t*)malloc((degree)*(sizeof(element_t)));
        int i;
        for(i = degree - 1; i >=0; i--){
                element_set(result[i], *ESubElements(*encryptZero(pk), ct->Coefficients[i],pk));
        }
        Ciphertext * cipher = (Ciphertext*)malloc(1*sizeof(Ciphertext));
        constructCiphertext(cipher,result, ct->Degree, ct->ScaleFactor, ct->L2);
        element_clear(*result);
        return cipher; 
}

element_t*  AInvElement(element_t el, PublicKey * pk){
        element_t * temp = (element_t*)malloc(sizeof(element_t));
        element_init_same_as(el, *ESubL2Elements(*encryptZero(pk),el,pk));
        element_init_same_as(*temp, el);
        return temp;
}

element_t* AInvElementL2(element_t el, PublicKey * pk){ 
        element_t * temp = (element_t*)malloc(sizeof(element_t));
        element_init_same_as(el, *ESubL2Elements(*ToDeterministicL2Element(*encryptZero(pk),pk),el,pk));
        element_init_same_as(*temp, el);
        return temp;
}

Ciphertext * EAdd(Ciphertext * ciphertext1, Ciphertext * ciphertext2, PublicKey * pk){
        if(ciphertext1->L2|| ciphertext2->L2){
                if(!ciphertext1->L2){
                        return EAddL2(MakeL2(ciphertext1,pk), ciphertext2,pk); 
                }
                if(!ciphertext2->L2){
                        return EAddL2(ciphertext1, MakeL2(ciphertext2,pk),pk);
                }
                ciphertext1 = EAddL2(ciphertext1, ciphertext2, pk);
                return ciphertext1;
        }
        Ciphertext * ct1 = (Ciphertext*)malloc(sizeof(Ciphertext));
        Ciphertext * ct2 = (Ciphertext*)malloc(sizeof(Ciphertext));
    
        ct1 = Copy(ciphertext1);
        ct2 = Copy(ciphertext2);
        
        alignCiphertexts(ct1, ct2, false,pk);

        
        int degree = (int)(max((double)(ct1->Degree), (double)(ct2->Degree))); 
        element_t * result = (element_t*)malloc((degree)*(sizeof(element_t))); 
        int i;
        for(i = 0; i < degree; i++){
                if(ct2->Degree > i && ct1->Degree > i){
                        element_t * temp = EAddElements(ct1->Coefficients[i], ct2->Coefficients[i], pk); 
                        element_init_same_as(result[i], *temp);
                        element_set(result[i], *temp);
                        continue; 
                }
                if(i >= ct2->Degree){
                        
                        element_init_same_as(result[i], ct1->Coefficients[i]);
                        element_set(result[i], ct1->Coefficients[i]);
                }
                if(i >= ct1->Degree){

                        element_init_same_as(result[i], ct2->Coefficients[i]);
                        element_set(result[i], ct2->Coefficients[i]);
                }
        }
        Ciphertext * ct = (Ciphertext*)malloc(sizeof(Ciphertext));
        constructCiphertext(ct, result, degree, ct1->ScaleFactor, ct1->L2);
        return ct;
}
mpz_t* DecodeSign(mpz_t m, PublicKey * pk){
        mpz_t * temp = (mpz_t*)malloc(sizeof(mpz_t));
        mpz_t threshold;
        mpz_t j;
        mpz_init(threshold);
        mpz_init_set_ui(j, 2);
        mpz_cdiv_q(threshold, pk->T, j);
        if((mpz_cmp(m, threshold)) >= 1 ){
                mpz_sub(m, m, pk->T);
        }
        mpz_init_set(*temp, m);
        gmp_printf("temp: %Zd\n",*temp);
        return temp;
}
Plaintext * Decrypt(Ciphertext * ct, PublicKey * pk, SecretKey * sk){
        if(ct->L2){
                return decryptL2(ct, pk, sk); 
        }
        int size;
        size = ct->Degree;
        int64_t * plaintextCoeffs = (int64_t *)malloc((size) * (sizeof(int64_t)));
        int i;
        
        for(i = 0; i < ct->Degree; i++){
                plaintextCoeffs[i] = (int64_t)mpz_get_d(*DecodeSign(*DecryptElement(ct->Coefficients[i], pk, false,sk), pk));
        }

        Plaintext * pt = (Plaintext*)malloc(sizeof(Plaintext));
        constructPlaintext(pt, pk, plaintextCoeffs, size, ct->ScaleFactor);
        return pt;
}
    
Ciphertext * aInvL2(Ciphertext * ct, PublicKey * pk){
        int degree = ct->Degree;
        element_t * result = (element_t*)malloc((degree)*(sizeof(element_t)));
        int i;
        for(i = degree-1; i >=0; i--){
                element_set(result[i], *ESubL2Elements(*encryptZeroL2(pk), ct->Coefficients[i], pk));
        }
        Ciphertext * c = (Ciphertext*)malloc(sizeof(Ciphertext));   
        constructCiphertext(c, result, ct->Degree, ct->ScaleFactor, ct->L2);
        return c;
}
mpz_t* DecryptElement(element_t  el, PublicKey * pk, bool failed, SecretKey *sk){
        element_t gsk;
        element_t csk;
        element_init_same_as(gsk, pk->G1);
        element_init_same_as(csk, pk->G1);
        element_set(gsk, pk->G1);
        element_set(csk, pk->G1);
        element_pow_mpz(gsk, pk->P, sk->Key);
        element_pow_mpz(csk, el, sk->Key);
        mpz_t * pt = RecoverMessageWithDL(gsk, csk, false, pk);
        return pt; 
}
mpz_t* DecryptElementL2 (element_t el, PublicKey * pk, bool failed, SecretKey *sk){
        element_t gsk; 
        element_t csk;
        pairing_t pair;
        pairing_init_pbc_param(pair, params);
        element_init_GT(gsk, pair);
        pairing_apply(gsk, pk->P, pk->P, pair);
       
        element_pow_mpz(gsk, gsk, sk->Key);
        element_init_same_as(csk, gsk); //el
        element_set (csk, gsk); //el
        element_pow_mpz(csk, el , sk->Key);
       
        mpz_t * pt = (mpz_t*)malloc(sizeof(mpz_t));
        mpz_init(*pt);
        mpz_set(*pt,*RecoverMessageWithDL(gsk, csk, false, pk));
        return pt;
}
Plaintext * decryptL2(Ciphertext * ct, PublicKey * pk, SecretKey * sk){ 
        int size;
        size = ct->Degree;
        int64_t * plaintextCoeffs = (int64_t *)malloc((size)*(sizeof(int64_t)));
        int i;
        
        for(i = 0; i < ct->Degree; i++){
                plaintextCoeffs[i] = (int64_t)mpz_get_d(*DecodeSign(*DecryptElementL2(ct->Coefficients[i], pk, false,sk), pk));
        }
        Plaintext * pt = (Plaintext *)malloc(sizeof(Plaintext));
        constructPlaintext(pt, pk, plaintextCoeffs, ct->Degree, ct->ScaleFactor);
        return pt;
}

Ciphertext * EAddL2(Ciphertext * ciphertext1, Ciphertext * ciphertext2, PublicKey *pk){
        Ciphertext * ct1 = (Ciphertext*)malloc(sizeof(Ciphertext));
        Ciphertext * ct2 = (Ciphertext*)malloc(sizeof(Ciphertext));
        
        alignCiphertexts(ct1,ct2, true, pk);

        int degree;
        degree = (int)(max((double)(ct1->Degree), (double)(ct2->Degree)));
        element_t * result = (element_t*)malloc(degree*(sizeof(element_t)));
        int i;
        for(i = degree-1; i >= 0; i--){
                if(i >= ct2->Degree){
                        element_set(result[i], ct1->Coefficients[i]);
                        continue;
                }
                if(i >= ct1->Degree){
                        element_set(result[i], ct2->Coefficients[i]);
                        continue;
                }
                element_set(result[i], *EAddL2Elements(ct1->Coefficients[i], ct2->Coefficients[i], pk));
        }
        Ciphertext * ct = (Ciphertext*)malloc(sizeof(Ciphertext));
        constructCiphertext(ct, result, degree, ct1->ScaleFactor, ct1->L2);
        return ct;
}
Ciphertext * EMultC(Ciphertext * ct, mpf_t constant, PublicKey * pk){
        if(ct->L2){       
                return eMultCL2(ct, constant, pk);
        }
        return eMultC(ct, constant,pk);
}
Ciphertext * eMultC(Ciphertext * ct, mpf_t constant, PublicKey *pk){

        bool isNegative = mpf_cmp_d(constant,0.0) < 0;
        if (isNegative){
                mpf_neg(constant,constant);
        }
        Plaintext * poly;
        poly = NewUnbalancedPlaintext(pk, mpf_get_d(constant));
        
        int degree = ct->Degree + poly->Degree;

        element_t *result = (element_t*)malloc(degree*sizeof(element_t));
        
        element_t zero;
        element_init_same_as(zero, pk->G1);
        int i;
        for(i = 0; i < degree; i++){
                element_set(result[i], zero);
        }
        int k;
        for(i = ct->Degree-1; i >= 0; i--){
                for(k = poly->Degree-1; k >=0; k--){
                        int index = i+k;
                        element_t coeff;
                        element_init_same_as(coeff, zero);
                        mpz_t temp;
                        mpz_init(temp);
                        mpz_set_d(temp, poly->Coefficients[k]);
                        element_set(coeff,*EMultCElement(ct->Coefficients[i],temp,pk));
                        element_set(result[index], *EAddElements(result[index], coeff, pk));
                }
        }
        
        Ciphertext * product = (Ciphertext*)malloc(sizeof(Ciphertext));
        constructCiphertext(product,result, degree, ct->ScaleFactor+poly->ScaleFactor, ct->L2);

        if (isNegative){
                return AInv(product, pk);
        }
        return product;
}

Ciphertext * eMultCL2(Ciphertext * ct, mpf_t constant, PublicKey *pk){
        
        bool isNegative = mpf_cmp_d(constant,0.0) < 0;
        if (isNegative){
                mpf_neg(constant,constant);
        }
        Plaintext * poly;
        poly = NewUnbalancedPlaintext(pk, mpf_get_d(constant));
        int degree = ct->Degree + poly->Degree;

        element_t *result = (element_t*)malloc(degree*sizeof(element_t));
        
        element_t zero;
        element_init_same_as(zero, pk->G1);
        int i;
        for(i = 0; i < degree; i++){
                element_init_GT(result[i],*pk->Pairing);
                element_init_same_as(result[i],result[i]);
        }
        int k;
        for(i = ct->Degree-1; i >= 0; i--){
                for(k = poly->Degree-1; k >=0; k--){
                        int index = i+k;
                        element_t coeff;
                        element_init_same_as(coeff, zero);
                        mpz_t temp;
                        mpz_init(temp);
                        mpz_set_d(temp, poly->Coefficients[k]);
                        element_set(coeff, *EMultCElement(ct->Coefficients[i],temp,pk));
                        element_set(result[index],*EAddElements(result[index], coeff, pk));
                }
        }
        
        Ciphertext * product = (Ciphertext*)malloc(sizeof(Ciphertext));

        constructCiphertext(product,result, degree, ct->ScaleFactor+poly->ScaleFactor, ct->L2);
        if (isNegative){
                return AInv(product,pk);
        }
        return product;
}


element_t* EMultCElement(element_t el, mpz_t constant, PublicKey * pk){
        element_t * res = (element_t*)malloc(sizeof(element_t)); 
        element_init_same_as(*res, el);
        element_pow_mpz(*res, el, constant);
        
        if(pk->Deterministic){
                mpz_t r;
                element_t q; 
                mpz_init(r);
                mpz_set(r, *newCryptoRandom(pk->N));
                element_init_same_as(q, el);
                element_mul_mpz(q, pk->Q, r); 
                element_mul(*res, *res, q);
        }
        return res;
}
Ciphertext * EMult(Ciphertext * ct1, Ciphertext * ct2, PublicKey * pk, Ciphertext * cipher, SecretKey * sk){
        clock_t t;
        t = clock();
        int degree, i, k, index; degree = ct1->Degree + ct2->Degree;
      
        element_t *coeff;
        element_t* result = (element_t*)malloc(degree*sizeof(element_t)); 
        pairing_t pair;
        pairing_init_pbc_param(pair, params); 
	element_t temp;
        for(i = 0; i < degree; i++){
                element_init_GT(temp, pair);
                element_init_same_as(result[i], temp);
                element_set(result[i], temp);
        }
         
        element_t * tempMult; element_t *tempAdd;
        for(i = ct1->Degree - 1; i >= 0; i--){
                for(k = ct2->Degree - 1; k >= 0; k--){
                       index = i + k;
                       coeff = EMultElements(ct1->Coefficients[i],ct2->Coefficients[k],pk);
                       
                       element_set(result[index],*EAddL2Elements(result[index], *coeff,pk));
                       
                }
        }
        constructCiphertext(cipher, result, degree, ct1->ScaleFactor + ct2->ScaleFactor, true);
        t = clock() - t;
	element_t temp11;	
	element_init_same_as(temp11,cipher->Coefficients[0]);
	printf ("It took me EMULT (%f seconds).\n",((float)t)/CLOCKS_PER_SEC);
        for(int i = 0; i < degree; i++)
                element_printf("element[%d]: %B \n", i, cipher->Coefficients[i]); 
	Plaintext * pt;
       // pt = Decrypt(cipher, pk, sk);
        return cipher;
}
element_t* EMultCElementL2(element_t el, mpz_t constant, PublicKey *pk){
        element_t * res = (element_t*)malloc(sizeof(element_t));
        element_init_same_as(*res,el);
        element_pow_mpz(*res, el, constant);
        
        if(!pk->Deterministic){
                mpz_t r;
                mpz_init(r);
                mpz_set(r, *newCryptoRandom(pk->N));
                element_t pair;
                element_pow_mpz(pair,pair,r);
                element_mul(*res, *res, pair);
        }
        return res;
}
element_t* EMultElements(element_t el1, element_t el2, PublicKey* pk){
        
        element_t * res= (element_t*)malloc(sizeof(element_t));
        pairing_t pair;
        pairing_init_pbc_param(pair, params);
        element_init_GT(*res, pair);
        element_init_same_as(*res,*res);
        pairing_apply(*res,el1,el2, pair);

        
        if(!pk->Deterministic){
                mpz_t r;
                mpz_init(r);
                mpz_init_set(r, *newCryptoRandom(pk->N));
                element_t pairs;
                element_init_GT(pairs, pair);
                pairing_apply(pairs,pk->Q,pk->Q,pair);
                element_pow_mpz(pairs,pairs,r);
                element_mul(*res, *res,pairs); 
        }
        return res;
}

Ciphertext * MakeL2(Ciphertext *ct, PublicKey *pk){
        Ciphertext * one;
        mpf_t temp;
        mpf_init(temp);
        mpf_set_d(temp,1.0);
        one = Encrypt(NewPlaintext(pk,temp), pk);
        return one;
}
element_t* toL2Element(element_t el, PublicKey * pk){
        element_t * result = (element_t*)malloc(sizeof(element_t));
        element_t pair;
        mpz_t r;
        mpz_init(r);
        mpz_set_si(r, 1); 
        element_init_GT(*result, *pk->Pairing);
        element_init_same_as(*result, *result);
        pairing_apply(*result, el, *EncryptElement(r,pk), *pk->Pairing);

        mpz_set(r, *newCryptoRandom(pk->N));
        element_init_GT(pair, *pk->Pairing);
        pairing_apply(pair, pk->Q, pk->Q, *pk->Pairing);
        element_pow_mpz(pair, pair, r);
        element_mul(*result, *result, pair);
        return result;
}
element_t* ToDeterministicL2Element(element_t el, PublicKey *pk){
        element_t * result = (element_t*)malloc(sizeof(element_t));
        element_init_GT(*result, *pk->Pairing);
        element_init_same_as(*result, *result);
        mpz_t temp;
        mpz_init(temp);
        mpz_set_ui(temp, 1);
        pairing_apply(*result,el,*EncryptDeterministic(temp,pk), *pk->Pairing);
        return result;
}
element_t* EncryptDeterministic(mpz_t x, PublicKey * pk){
        element_t * G = (element_t*)malloc(sizeof(element_t));
        element_init_same_as(*G, pk->G1);
        element_pow_mpz(*G, pk->P, x);
        return G;
}

element_t* EncryptElement(mpz_t x, PublicKey *pk){
        element_t G, H;
        element_t * C = (element_t*)malloc(sizeof(element_t));
        element_init_same_as(G, pk->G1);
        element_pow_mpz(G, pk->P, x);
        
        mpz_t r; 
        mpz_init(r);
        mpz_set(r,*newCryptoRandom(pk->N));
        element_init_same_as(H, pk->G1);
        element_pow_mpz(H, pk->Q, r);

        element_init_same_as(*C, pk->G1);
        element_mul(*C, H, G);
        return C;
}

mpz_t* RecoverMessageWithDL(element_t gsk, element_t csk, bool l2, PublicKey *pk){
        element_t zero;
        element_init_same_as(zero, gsk);
        mpz_t * m = (mpz_t*)malloc(sizeof(mpz_t));
        mpz_init(*m);  
        
        if(!element_cmp(zero,csk)){
                return m;
        } 
        
	/*Function*/
	element_t g0;
	
  	element_init_same_as(g0, gsk);

  	element_set(g0, gsk);
  	mpz_set_ui(*m, 1);
  	while (element_cmp(g0, csk)) {
        element_mul(g0, g0, gsk);
    	mpz_add_ui(*m, *m, 1);
  	}
        
        return m;
}

element_t* ESubElements(element_t coeff1, element_t coeff2, PublicKey * pk){
        element_t * result;
        element_init_same_as(*result, pk->G1);
        element_div(*result, coeff1, coeff2);
        if(pk->Deterministic){
                return result;
        }
        element_t h1;
        mpz_t rand;
        mpz_init(rand);
        mpz_set(rand,*newCryptoRandom(pk->N));

        element_init_same_as(h1,pk->G1);
        element_pow_mpz(h1,pk->Q,rand);
        element_mul(*result, *result, h1);
        return result;
}
element_t* ESubL2Elements(element_t coeff1, element_t coeff2, PublicKey * pk){
        element_t * result = (element_t*)malloc(sizeof(element_t));
        element_init_GT(*result, *pk->Pairing);
        element_init_same_as(*result, *result);
        element_div(*result, coeff1, coeff2);
        if(pk->Deterministic){
                return result;
        }
        mpz_t r;
        mpz_init(r);
        mpz_set(r,*newCryptoRandom(pk->N));

        element_t pair;
  
        element_init_GT(pair, *pk->Pairing);
        pairing_apply(pair, pk->Q, pk->Q, *pk->Pairing);
        element_pow_mpz(pair, pair, r);
        
        /*
        This is conversion for pair.PowBig(pair, r); 
        and result Mult(result pair)
        */

        element_mul(*result, *result, pair);
        return result;
}
element_t* EAddElements(element_t coeff1, element_t coeff2, PublicKey *pk){
        element_t * result = (element_t*)malloc(sizeof(element_t));
        element_t h1;
        element_init_same_as(*result, pk->G1);
        element_set(*result, *result);
	element_mul(*result, coeff1, coeff2);
        if(pk->Deterministic){
                return result;
        }
        mpz_t rand; 
        mpz_init(rand);
        mpz_set(rand,*newCryptoRandom(pk->N));
        
        element_init_same_as(h1,pk->G1);
        element_pow_mpz(h1,pk->Q,rand);
        element_mul(*result, *result, h1);
        return result;
}
element_t* EAddL2Elements(element_t coeff1, element_t coeff2, PublicKey *pk){
        element_t * result = (element_t*)malloc(sizeof(element_t));
        pairing_t pairing;
        pairing_init_pbc_param(pairing, params);
        element_init_GT(*result, pairing);
        element_init_same_as(*result, *result);
        element_set(*result, *result);
        element_mul(*result, coeff1, coeff2);
        if(pk->Deterministic){
                return result;
        }
        mpz_t r;
        mpz_init(r);
        mpz_set(r,*newCryptoRandom(pk->N));
        
        element_t pair;
        element_init_GT(pair, pairing);
        pairing_apply(pair, pk->Q, pk->Q, pairing); 

        element_pow_mpz(pair, pair, r);
        
        element_mul(*result, *result, pair);
        return result;
        
}
element_t* EPolyEval(Ciphertext * ct, PublicKey * pk){
        element_t * acc = (element_t*)malloc(sizeof(acc));
        mpz_t x;
        mpz_init(x);
        element_set(*acc, *EncryptDeterministic(x,pk));
        mpz_set_ui(x, pk->PolyBase);
        int i;
        for(i = ct->Degree - 1; i >= 0; i--){
                element_set(*acc,*EMultCElement(*acc,x,pk));
                element_set(*acc, *EAddElements(*acc, ct->Coefficients[i],pk));
        }
        return acc;
                         
}       
void alignCiphertexts(Ciphertext * ct1, Ciphertext * ct2, bool level2, PublicKey * pk) {
        int diff;
        if(ct1->ScaleFactor > ct2->ScaleFactor){
                diff = ct1->ScaleFactor - ct2->ScaleFactor;
                mpf_t temp;
                mpf_init(temp);
                mpf_set_d(temp,pow((double)pk->FPScaleBase,(double)diff));
                ct2 = EMultC(ct2,temp,pk);
                ct2->ScaleFactor = ct1->ScaleFactor;
        }
        else if(ct1->ScaleFactor < ct2->ScaleFactor){
                alignCiphertexts(ct2,ct1,level2,pk);
        }  
}
element_t *encryptZero(PublicKey *pk){
        mpz_t temp;
        mpz_init(temp);
        return EncryptElement(temp,pk);
}
element_t *encryptZeroL2(PublicKey *pk){
        element_t zero; 
        element_t * result = (element_t*)malloc(sizeof(element_t));
        element_set(zero, *encryptZero(pk));
           
        element_init_GT(*result, *pk->Pairing);
        element_init_same_as(*result, *result); 
        pairing_apply(*result, zero, zero, *pk->Pairing);
        
        mpz_t * r;
 
        r = newCryptoRandom(pk->N);

        element_t pair;
        element_init_GT(pair, *pk->Pairing);
        pairing_apply(pair, pk->Q, pk->Q, *pk->Pairing);
        element_pow_mpz(pair, pair, *r);

        element_mul(*result,*result,pair);
        return result;
}
mpz_t * newCryptoRandom(mpz_t max){
        uint32_t randseed; /*generates new random number */
        mpz_t * rand = (mpz_t*)malloc(sizeof(mpz_t));
        mpz_init(*rand); 
        randseed = randombytes_random();
        gmp_randstate_t gmprand;
        gmp_randinit_mt(gmprand);
        gmp_randseed_ui(gmprand, randseed);
        mpz_urandomm(*rand, gmprand, max);
        gmp_randclear(gmprand);
        return rand;
}
mpz_t *parseLFromPBCParams(pbc_param_t params){
        /*
        This is workaround for parsing l from the PBC params. 
        Still looking for a more secure method, but this will work for 
        researching purposes.
        */
        freopen("output.txt", "a+", stdout);
        char find_l;
        FILE* outfile;
        pbc_param_out_str(stdout, params);
        freopen("/dev/tty", "w", stdout);
        outfile = fopen("output.txt", "r");
        int64_t l = 0;
        while(find_l != EOF){
                find_l = getc(outfile);
                if(find_l == 'l'){
                      fscanf(outfile, "%lld", &l);
                }
        }
        fclose(outfile);
        mpz_t * returnL = (mpz_t *)(malloc(sizeof(mpz_t)));
        mpz_init_set_d(*returnL, l);
        return returnL;
}
