#ifndef JLMD_CELL_H
#define JLMD_CELL_H

#include "ljmd.h"
#include <stdlib.h>
#include <vector>

// cell structure
struct _mdcell {
    std::vector<int> idx;
};

typedef struct _mdcell cell_t;

int index3d(mdsys_t *sys, int i, int j, int k);

void sort(mdsys_t *sys, cell_t *cel);

#endif
