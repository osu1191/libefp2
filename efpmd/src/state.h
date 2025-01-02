#ifndef STATE_H
#define STATE_H

#include "../torch/c_libtorch.h"

typedef struct {
    ANIModel* model;
} ANIState;

extern ANIState global_state;

#endif // STATE_H
