//
//  SetCover.h
//  graphcluster
//
//  Created by Martin Steinegger on 21.05.13.
//  Copyright (c) 2013 Martin Steinegger. All rights reserved.
//

#ifndef __graphcluster__SetCover__
#define __graphcluster__SetCover__

#include <iostream>
#include <list>
#include "SetElement.h"
#include "LinearMultiArray.h"
class set_cover {
public:
    set_cover(unsigned int set_size,
              unsigned int element_size,
              unsigned int weight_range,
              unsigned int all_element_count,
              unsigned int * element_size_lookup);
    ~set_cover();

    void add_set(const int set_id, const int set_weight,
                 const unsigned int * element_ids,
                 const unsigned short * weights,
                 const int element_size);
    std::list<set_ *> execute_set_cover();
/*
    get_highest_weighted_set 
    input 
            int start_pos is the search start position 
    output
            returns the highest sets from a certaint position
*/
    set_* get_highest_weighted_set(int start_pos);
private:
    unsigned int add_position;
    int element_size;
    int set_size;
    int weight_range;
    set_** ordered_by_score_set;
    set_::element * set_elements;
    set_* sets;

    linear_multi_array<set_::element *> * element_lookup;
    linear_multi_array<set_::element *> * set_element_lookup;
    // methodes
    void removeSet(set_* s);
    set_::element * unplug_element(set_::element * element_to_unplug, set_::element * first_element);
    void unplug_set(set_* set_to_remove);
    set_* create_set_at_weight_position(int weight, set_* set_to_add);

};
#endif /* defined(__graphcluster__SetCover__) */
