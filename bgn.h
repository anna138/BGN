//
// bgn.h
//
#ifndef _BGN_H_
#define _BGN_H_

#include <gmp.h>
#include <sodium.h>
#include <pbc/pbc.h>

#include "ciphertext.h"

//pbc_param_t params;

typedef struct PKey {
    pairing_t *Pairing; /*pairing between G1 and G2*/
    element_t G1; /*G1 Group*/
    element_t P; /*generator of G1*/
    element_t Q; /*generator Q*/
    mpz_t N; /*product of two primes*/
    mpz_t T; /*message space T*/
    int PolyBase; /*ciphertext polynomial encoding base */
    int FPScaleBase; /*fixed point encoding scale base */
    double FPPrecision; /* min error tolerance for fixed point encoding*/
    bool Deterministic; /* whether or not the homomorphic operations                                                           are deterministic*/
} PublicKey;

typedef struct SKey {
    mpz_t Key; /*big integer for the key*/
    int PolyBase;
} SecretKey;

typedef struct Ptext {
    PublicKey *Pk;
    int64_t Coefficients[512];  //may need to change this to dynamic
    int Degree;
    int ScaleFactor;
} Plaintext;

void cleanParams();
void constructPlaintext(Plaintext *that, PublicKey *pk, int64_t *coeff,
                        int degree, int sfactor);
void constructPublicKey(PublicKey *that, pairing_t* par, element_t g1,
                        element_t p1, element_t q1, mpz_t n1, mpz_t t1, int PB,
                        int FPS, double FPP, bool DTRM);
void constructSecretKey(SecretKey *that, mpz_t s1, int PB);
void destructSecretKey(SecretKey *that);
void destructPublicKey(PublicKey *that);
double max(double x, double y);
void NewKeyGen(int keyBits, mpz_t T, int polyBase, int fpScaleBase,
               double fpPrecision, bool deterministic, PublicKey * pk,
               SecretKey * sk);
Ciphertext * Encrypt(Plaintext *pt, PublicKey * pk);
Ciphertext * AInv(Ciphertext *ct, PublicKey * pk);
element_t* AInvElement(element_t el, PublicKey * pk);
element_t* AInvElementL2(element_t el, PublicKey * pki);
Ciphertext * EAdd(Ciphertext * ciphertext1, Ciphertext * ciphertext2,
                  PublicKey * pk);
mpz_t* DecodeSign(mpz_t m, PublicKey * pk);
Plaintext * Decrypt(Ciphertext * ct, PublicKey * pk, SecretKey * sk);
Ciphertext * aInvL2(Ciphertext * ct, PublicKey * pk);
mpz_t* DecryptElement(element_t el, PublicKey * pk, bool failed, SecretKey *sk);
mpz_t* DecryptElementL2(element_t el, PublicKey * pk, bool failed,
                        SecretKey *sk);
Plaintext * decryptL2(Ciphertext * ct, PublicKey * pk, SecretKey * sk);
Ciphertext * EAddL2(Ciphertext * ciphertext1, Ciphertext * ciphertext2,
                    PublicKey *pk);
Ciphertext * EMultC(Ciphertext * ct, mpf_t constant, PublicKey * pk);
Ciphertext * eMultC(Ciphertext * ct, mpf_t constant, PublicKey *pk);
Ciphertext * eMultCL2(Ciphertext * ct, mpf_t constant, PublicKey *pk);
element_t* EMultCElement(element_t el, mpz_t constant, PublicKey * pk);
Ciphertext * EMult(Ciphertext * ct1, Ciphertext * ct2, PublicKey * pk,
                   Ciphertext * cipher, SecretKey *sk);
element_t* EMultCElementL2(element_t el1, mpz_t constant, PublicKey *pk);
element_t* EMultElements(element_t el1, element_t el2, PublicKey* pk);
Ciphertext * MakeL2(Ciphertext *ct, PublicKey *pk);
element_t* toL2Element(element_t el, PublicKey * pk);
element_t* ToDeterministicL2Element(element_t el, PublicKey *pk);
element_t* EncryptDeterministic(mpz_t x, PublicKey * pk);
element_t* EncryptElement(mpz_t x, PublicKey *pk);
mpz_t *RecoverMessageWithDL(element_t gsk, element_t csk, bool l2,
                            PublicKey *pk);
element_t* ESubElements(element_t coeff1, element_t coeff2, PublicKey * pk);
element_t* ESubL2Elements(element_t coeff1, element_t coeff2, PublicKey * pk);
element_t* EAddElements(element_t coeff1, element_t coeff2, PublicKey *pk);
element_t* EAddL2Elements(element_t coeff1, element_t coeff2, PublicKey *pk);
element_t* EPolyEval(Ciphertext * ct, PublicKey * pk);
void alignCiphertexts(Ciphertext * ct1, Ciphertext * ct2, bool level2,
                      PublicKey * pk);
element_t* encryptZero(PublicKey *pk);
element_t* encryptZeroL2(PublicKey *pk);
mpz_t *newCryptoRandom(mpz_t max);
mpz_t *parseLFromPBCParams(pbc_param_t params);

#endif // _BGN_H_
