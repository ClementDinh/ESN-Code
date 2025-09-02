#pragma once
#include <cstdio>

typedef struct{
  double val;
  int x, y, left, right, bottom, top, topLeft, topRight, bottomLeft, bottomRight;
} Index;

int mod(int x, int m);
void fillIndex(Index *lat, int L);
// void printLat(Index *lat, int L);
void writeLat(Index *lat, int L, FILE *data);
double corr(int x, Index* lat, int L); // calculate and add correlation func value to an arr
double corr2(int x, Index* lat, int L);
void updateLat(Index *lat, int L);
int* neighbor(Index* lat, int L, int point);

