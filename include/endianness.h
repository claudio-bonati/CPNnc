#ifndef ENDIANNESS_H
#define ENDIANNESS_H

#include<stdio.h>

int endian(void); // return 0 if little endian
void SwapBytesInt(void *pv);
void SwapBytesFloat(void *pv);
void SwapBytesDouble(void *pv);

int print_on_binary_file_bigen_double(FILE *fp, double r);
int read_from_binary_file_bigen_double(FILE *fp, double *r);

#endif
