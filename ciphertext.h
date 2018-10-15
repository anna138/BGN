//
// ciphertext.h
//
#ifndef CIPHERTEXT
#define CIPHERTEXT

#include <gmp.h>
#include <sodium.h>
#include <pbc/pbc.h>

typedef struct Ctext {
    element_t Coefficients[512];  //may need change this to dynamic
    int Degree;
    int ScaleFactor;
    bool L2;
} Ciphertext;

void constructCiphertext(Ciphertext *that, element_t *coeff, int degree,
                         int sfactor, bool l2);
void destructCiphertext(Ciphertext * that);
Ciphertext * Copy(Ciphertext * ct);
Ciphertext * NewCiphertext(element_t * coefficients, int degree,
                           int scaleFactor, bool l2);

#endif
