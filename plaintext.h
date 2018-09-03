/*This is plaintext.h*/
#ifndef PLAINTEXT
#define PLAINTEXT
#include "bgn.h"

/*Prototypes in Code */
Plaintext * NewUnbalancedPlaintext(PublicKey *pk, float m);
Plaintext * NewPlaintext(PublicKey * pk, mpf_t m);
void computeEncodingTable(PublicKey * pk);
int degree(mpz_t target, mpz_t * sums, int bound, bool balanced);
int64_t * unbalancedEncode(mpz_t target, int base, mpz_t * degrees, mpz_t * sumDegrees, int * encodeInt);
int64_t * balancedEncode(mpz_t target, int base, mpz_t * degrees, mpz_t * sumDegrees, int * encodeInt);
int64_t * reverse(int64_t * numbers, int size);
void rationalize(double x, int base, double precision, int64_t * numb, int * pew);
mpf_t * PolyEval(Plaintext * p);
bool checkOverflow(mpz_t x);
std::string String(Plaintext * p);
void cleanTablrus();
#endif
