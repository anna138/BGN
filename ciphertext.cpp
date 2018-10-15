//
// ciphertext.cpp
//
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>

#include "ciphertext.h"
#include "bgn.h"

#include <gmp.h>
#include <sodium.h>
#include <pbc/pbc.h>

void constructCiphertext(Ciphertext *that, element_t *coeff, int degree,
                         int sfactor, bool l2)
{
    int i;

    for(i = 0; i < degree; i++) {
        element_init_same_as(that->Coefficients[i], coeff[i]);
        element_set(that->Coefficients[i], coeff[i]);
    }  //may change this to dynamic

    that->Degree = degree;
    that->ScaleFactor = sfactor;
    that->L2 = l2;
}
void destructCiphertext(Ciphertext *that)
{
    int i;
    for(i = 0; i < that->Degree; i++) {
        element_clear(that->Coefficients[i]);
    }
}
Ciphertext * Copy(Ciphertext * ct)
{
    Ciphertext * c = (Ciphertext*) malloc(sizeof(Ciphertext));
    constructCiphertext(c, ct->Coefficients, ct->Degree, ct->ScaleFactor,
                        ct->L2);
    return c;
}
Ciphertext * NewCiphertext(element_t * coefficients, int degree,
                           int scaleFactor, bool l2)
{
    Ciphertext * c = (Ciphertext*) malloc(sizeof(Ciphertext));
    constructCiphertext(c, coefficients, degree, scaleFactor, l2);
    return c;
}

