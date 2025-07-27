//
//  element.h
//  graphcluster
//
//  Created by Martin Steinegger on 21.05.13.
//  Copyright (c) 2013 Martin Steinegger. All rights reserved.
//

#ifndef __graphcluster__element__
#define __graphcluster__element__
#include <iostream>
#include "SetElement.h"

struct element_meta_data{
    int element_id;
    int count;
};


struct set_ {
    int set_id;
    unsigned short weight;
    struct element {
        unsigned int element_id;
        unsigned short weight;
        element(unsigned int element_id,unsigned short weight) : element_id(element_id), weight(weight){};
        set_ * parent_set;
        element * last;
        element * next;
    };
    element * elements;
    set_ * last;
    set_ * next;
};



#endif /* defined(__graphcluster__element__) */
