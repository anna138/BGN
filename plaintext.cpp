/*This is plaintext.c*/
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <gmp.h>
#include <pbc.h>
#include <sodium.h>
#include <math.h>
#include <cstdlib>
#include <iostream>
#include "plaintext.h"

mpz_t * degreeSumTable;
mpz_t * degreeTable;
int computedBase;


const int degreeBound = 10;

/*
Will be adding NTRU polynomial conversion function on plaintext.c
*/
void cleanTablrus(){
        for(int i = 0; i < degreeBound; i++){

                mpz_clear(degreeSumTable[i]);
                mpz_clear(degreeTable[i]);
        }
        free(degreeSumTable);
        free(degreeTable);
}

Plaintext * NewUnbalancedPlaintext(PublicKey *pk, float m){
        if(degreeTable == 0){
                perror("Encoding tables are not computed!");
        }
        srand(time(NULL));
        double mFloat = (double)rand() / (double)RAND_MAX;
        int64_t numerator; int scaleFactor;
        int64_t * coeffs; int degree;

        Plaintext * pt = (Plaintext*)malloc(sizeof(Plaintext));
        /*
        This is only for floating point numbers. We don't need it for Genomic Data

        if (fmod(mFloat, 1.0) !=0.0){
                mFloat = (double)rand() / (double)RAND_MAX;
                rationalize(floor(mFloat), pk->FPScaleBase, pk->FPPrecision, &numerator, &scaleFactor);
                mpz_t mInt;  
               
                mpz_t temp;
                mpz_init_set_d(temp, numerator);
                mpz_add(mInt, mInt, temp);
                //gmp_printf("encoded %f as %s/%s\n", m, mInt 
                
                coeffs = unbalancedEncode(mInt, pk->PolyBase, degreeTable, degreeSumTable,  &degree);
                constructPlaintext(pt, pk, coeffs,degree,scaleFactor);
                return pt;
        }
        */
                mpz_t mInt;
                mpz_init(mInt);
                mpz_init_set_d(mInt, m);
                coeffs = unbalancedEncode(mInt, pk->PolyBase, degreeTable, degreeSumTable,  &degree);
                constructPlaintext(pt, pk, coeffs,degree, 0);
                return pt;

                
}

Plaintext * NewPlaintext(PublicKey * pk, mpf_t m){
        if(degreeTable == 0){
                perror("Encoding tables are not computed");
        }

        double mFloat = (double)rand() / (double)RAND_MAX;
        
        Plaintext * pt = (Plaintext*)malloc(sizeof(Plaintext));

        int64_t * coeffs; int degree;
        /*
        if (fmod(mFloat, 1.0) != 0.0){
                int64_t numerator; int scaleFactor;
                
                rationalize(floor(mFloat), pk->FPScaleBase, pk->FPPrecision, &numerator, &scaleFactor);
                mpz_t mInt;  
                mpz_t tmp;
                double tmp2 = pow((double)pk->FPScaleBase, (double)scaleFactor);
                mpz_init_set_d(tmp, tmp2);
                mpz_mul(mInt, mInt, tmp); 
                mpz_t temp;
                mpz_init_set_d(temp, numerator);
                mpz_add(mInt, mInt, temp);
                printf("encoded"); 
                
                coeffs = balancedEncode(mInt, pk->PolyBase, degreeTable, degreeSumTable, &degree);
                constructPlaintext(pt, pk, coeffs,degree,scaleFactor);
                return pt;
       }
       */

                mpz_t mInt;
                mpz_init(mInt);
                mpz_set_f(mInt, m);
                coeffs = balancedEncode(mInt, pk->PolyBase, degreeTable, degreeSumTable, &degree);
                //printf("DEGREE: %d\n",degree);
                constructPlaintext(pt, pk, coeffs,degree, 0);
                return pt;
}
void computeEncodingTable(PublicKey * pk){
        mpz_t base;
        mpz_init_set_si(base, pk->PolyBase);
        int bound = degreeBound;
        degreeTable = (mpz_t*)malloc(bound*(sizeof(mpz_t)));
        degreeSumTable = (mpz_t*)malloc(bound*(sizeof(mpz_t)));
        mpz_t sum;
        mpz_init_set_si(sum, 1);
        mpz_init_set_si(degreeSumTable[0], 1);
        //gmp_printf("\n106- DegreeSum: %Zd \n", degreeSumTable[0]);
        mpz_init_set_si(degreeTable[0], 1);
        int i;
        mpz_t result;
        for(i = 1; i < bound; i++){
                mpz_init(result);
                /*
                temp2 was to be mpz value of 0, 
                find base to power of i mod 0 but now removing temp2
                
                bc mpz_powm(result, base, temp, temp2): changed, 
                bc result := big.NewInt(0).Exp(base, big.NewInt(int64(i)), nil) 
                is always nil so x **y
                */
                mpz_pow_ui(result, base, i);
                mpz_add(sum, sum, result);
                mpz_init_set(degreeTable[i], result);
                mpz_init_set_si(degreeSumTable[i], 0);
                mpz_init_set(degreeSumTable[i], sum);
                //gmp_printf("125- DegreeSum: %Zd \n", degreeSumTable[i]);
        }
        mpz_clear(base);
        mpz_clear(sum);
        mpz_clear(result);
        //gmp_printf("This is the end of the function");
        
}
int degree(mpz_t target, mpz_t * sums, int bound, bool balanced){
       
        if(mpz_cmp_ui(target, 1) == 0){
                return 0;
        }
        if(balanced){
                int i;
                for(i = 1; i <=bound; i++){
                        //gmp_printf("The degreeSumTable[i]: %Zd \n Target: %Zd \n", degreeSumTable[i], target);
                        if(mpz_cmp(degreeSumTable[i], target) >= 0){
                                return i;
                        }
                }
        }
        else{
                int i;
                for(i = 1; i <=bound; i++){
                        if(mpz_cmp(degreeTable[i], target) >= 1){
                                return i-1;
                        }
                }
        }
        //printf("MINUS ONE Line 152\n");
        return -1;
}
int64_t * unbalancedEncode(mpz_t target, int base, mpz_t * degrees, mpz_t * sumDegrees, int * encodeInt){ 
        /* original function return two pm:  array and encodeInt*/
        mpz_t zero;
        mpz_init(zero);/* This is the Special Case*/
        if((mpz_cmp(target, zero)) == 0){
                int64_t * coefficients = (int64_t*)malloc(sizeof(int64_t));
                coefficients[0] = 0;
                *encodeInt = 1;
                return coefficients;
        }
        if((mpz_cmp(target, zero)) < 0){
                perror("Negative encoding not support");
                //exit(1);
        }
        if(!sumDegrees){        
                perror("No precomputed degree table!");
        }
        
        int64_t * coefficients = (int64_t*)malloc(degreeBound*sizeof(int64_t));
        int bound = (degreeBound); /*what is the length of mpz_t*/
        int lastDegree = degreeBound;
        int index = 0;
        while(1){
                //gmp_printf("degree(%Zd, sumDegrees, %d, false);", target, lastDegree);
                index = degree(target, sumDegrees, lastDegree, false);
                //gmp_printf("This is index after degree function: %d", index);
                lastDegree = index +1;
                //gmp_printf("This is last degree: %d", lastDegree);
                if(bound == (degreeBound)){ 
                        /*
                        SumDegrees was declared to be the the length 
                        is the degreeBound declared in compute
                        */
                        bound = index + 1;
                }
                mpz_t value;
                mpz_init(value);
                mpz_t value2;
                mpz_init(value2);
                mpz_init_set(value, degrees[index]);
                mpz_t two;
                mpz_init(two);
                mpz_init_set_si(two, 2);
                mpz_t degIndx;
                mpz_init(degIndx);

                mpz_init_set(value, degrees[index]);
                mpz_mul(value2, degIndx, two);

                if(mpz_cmp(value, target)){
                        mpz_set(value, value2);
                }
                else{
                        coefficients[index] = 1;
                }
                if(mpz_cmp(value, target)){
                        *encodeInt = bound + 1;
                        return &coefficients[*encodeInt];
                }
                mpz_sub(target, target, value);
        }
}
int64_t * balancedEncode(mpz_t target, int base, mpz_t * degrees, mpz_t * sumDegrees, int * encodeInt){
        mpz_t zero;
        mpz_init(zero);
        mpz_t negone;
        mpz_init_set_si(negone, -1);
        //gmp_printf("Plaintext.c Line:221 This is target: %Zd \nThis is zero %Zd \n", target, zero);
        /* This is the Special Case*/
        if((mpz_cmp(target, zero)) == 0){
                int64_t * coefficients = (int64_t*)malloc(sizeof(int64_t));
                coefficients[0] = 0;
                *encodeInt = 1;
                return coefficients;
        }
        /*End of Special Case*/
        bool isNegative = (mpz_cmp(zero, target) > 0);
        if (isNegative){
                mpz_mul(target,target, negone);
                
        }
        if(!sumDegrees){        
                perror("No precomputed degree table!");
                exit(1);
        }
        
        int64_t * coefficients = (int64_t*)malloc(degreeBound*sizeof(int64_t));
        int i;
        for(i = 0; i < degreeBound; i++){
                coefficients[i] = 0;
        }
        int bound = degreeBound; 
        int lastIndex = degreeBound;
        int index;
        bool nextNegative = false;
        while(1){
                //gmp_printf("Plaintext.c Line 245: Degree(%Zd, sumDegrees, %d, true);\n\n", target, lastIndex);
                index = degree(target, sumDegrees, lastIndex, true);
                //gmp_printf("This is index after degree function: %d \n\n", index);
                lastIndex = index;
                //gmp_printf("LastIndex: %d \n\n\n", lastIndex);
                if(bound == (degreeBound)){ 
                        /*sumDegrees is the length of degreeBound*/
                        bound = index;
                }
                //gmp_printf("Is this index %d \n\n", index);
                coefficients[index] = 1;
                if(nextNegative){
                        coefficients[index] *= -1;
                }
                if((mpz_cmp(degrees[index], target)) == 0){
                        if(isNegative){
                                for(int i = 0; i <= bound; i++){
                                        coefficients[i] *= -1;
                                        //gmp_printf("Coefficient: %lld", coefficients[i]);
                                }
                        }
                        *encodeInt = bound + 1;
                        return coefficients;//&coefficients[*encodeInt]; commenting out to return full coefficients
                }
                if(mpz_cmp(degrees[index], target) >= 1){
                        //gmp_printf("this is target: %Zd \n this is index: %d \n\n this is degrees %Zd \n\n\n", target, index, degrees[index]);
                        nextNegative = !nextNegative;
                        mpz_sub(target, degrees[index], target);
                }
                else{
                        //gmp_printf("this is target: %Zd \n this is index: %d \n\n this is degrees %Zd \n\n\n", target, index, degrees[index]);
                        mpz_sub(target, target, degrees[index]);
                }
                //gmp_printf("balanced Encode: this is target: %Zd \n this is degrees %Zd \n this is index: %d \n\n\n", target, degrees[index], index);
        }

}

int64_t * reverse(int64_t * numbers, int size){
        int j;
        for(int i = 0; i < size; i++){ /*length of numbers*/                
                j = size - i - 1;
                numbers[i] = numbers[j];
                numbers[j] = numbers[i];
        }
        return numbers;
}

/* This is only used for floating points. Since we use genomic data, no need*/

void rationalize(double x, int base, double precision, int64_t * numb, int * pew){
        int factor = floor(x);

        x = 1.0 + fmod(x, 1.0);

        if(abs(x) > 1.0){
                x += 1.0;
        }
        if(x >= 0.0){
                x -= (double)((int)x);
        }
        else if(x <= 0.0){
                x += (double)((int)x);
        }
        //printf("\n x is: %f \n Precizion is: %f \n\n", x, precision);
        double num = 1.0;
        double temp = 1.0; 
        /*this is pew double to later be casted as int, 
        pew variable is used as work around to return*/
        x = 1; //testing
        double qmin = x - precision;
        double qmax = x + precision;
         
        double denom = 0.0;
        double rat = 0.0;
        int loop = 0;

        //printf("QMin: %f \n QMax: %f \n Demon: %f \n Rat: %f \n", qmin,qmax,denom,rat); 
        while(1){
                denom = pow((double)(base), temp);
                rat = num/denom;

                //printf("QMin: %f \n QMax: %f \n Demon: %f \n Rat: %f \n", qmin,qmax,denom,rat);
                if(rat <= qmax && rat >= qmin){
                        while(((int)(num) % base) == 0){
                                num = num/(double)(base);
                                temp--;
                        }
                        denom = pow((double)base, temp);
                        *pew = (int)temp;
                        *numb = (int64_t)(factor*denom + num);
                        return;
                }
                //*pew = (int)temp;
                if(num+1 >= denom){
                        num = (double)1;
                        temp = temp + 1.0;
                        //*pew++;
                }
                loop++;
                num++; //num++
                //printf("Loop: %d \n\n QMin: %f \n QMax: %f \n Demon: %f \n Rat: %f \n Num: %f \n Base: %d \n Temp: %f, Pew: %d", loop,qmin,qmax,denom,rat,num, base, temp, *pew); 
        }
        //gmp_printf("This is the end of Rationalize");
}
mpf_t * PolyEval(Plaintext * p){
        mpf_t * acc = (mpf_t*)malloc(sizeof(mpf_t));
        mpf_init(*acc);
        mpf_t x;
        mpf_init(x);
        mpf_set_si(x, (double)p->Pk->PolyBase);
        double var1 = 0;
        mpf_t temp;
        mpf_init(temp);
        int i;
        for(i = p->Degree-1; i >= 0; i--){
                mpf_mul(*acc, *acc, x);
                mpf_set_d(temp, p->Coefficients[i]);
                mpf_add(*acc, *acc, temp);
        }
        //workaround for trailing zeros problem
        for(i = p->Degree-1; i >=0; i--){
                var1 = (var1 * p->Pk->PolyBase) + p->Coefficients[i];
        } 
        std::cout << "This is var: " << var1 << std::endl;

        return acc;
}

bool checkOverflow(mpz_t x){
        mpz_t max;
        mpz_init_set_d(max, 9223372036854775807); 
        /*boolean test for overflow*/
        return (mpz_cmp(x, max) > 0);
}

std::string String(Plaintext * p){
        mp_exp_t exponent;
        std::string str;
        str = mpf_get_str(NULL, &exponent, 10, 0, *PolyEval(p));   
        return str;
}

