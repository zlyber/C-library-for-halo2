#ifndef __SETUP_H__
#define __SETUP_H__

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdint.h>
#include"curve.h"

/// Initializes parameters for the curve, draws toxic secret from given rng.
/// MUST NOT be used in production.
void setup(int k);

#endif