//
//  SetCover.cpp
//  graphcluster
//
//  Created by Martin Steinegger on 21.05.13.
//  Copyright (c) 2013 Martin Steinegger. All rights reserved.
//

#include "SetCover.h"
#include <string.h> // memset


set_cover::set_cover(unsigned int set_size,
                     unsigned int uniqu_element_size,
                     unsigned int weight_range,
                     unsigned int all_element_count,
                     unsigned int * element_size_lookup){
    this->set_size = set_size;
    this->element_size = uniqu_element_size;
    this->weight_range = weight_range;
    this->ordered_by_score_set = new set_ *[weight_range+1]; // score range
    memset(this->ordered_by_score_set, 0, sizeof(set_*)*(weight_range+1)); // init with 0 (no element is set)
    
    this->element_lookup     = new linear_multi_array<set_::element*>(element_size_lookup,uniqu_element_size,all_element_count);
    this->set_elements = (set_::element *) malloc(sizeof(set_::element)*all_element_count);
    
    this->sets = (set_*) malloc(sizeof(set_)*(set_size+1));
    this->add_position = 0;
}

set_cover::~set_cover(){
    free(set_elements);
    free(sets);
    delete element_lookup;
    delete[] ordered_by_score_set;
}


set_* set_cover::create_set_at_weight_position(int weight, set_* set_to_add){
    set_* weighted_position_start_set=this->ordered_by_score_set[weight];
    
    if(weighted_position_start_set == 0) { // first element is not yet set
        set_to_add->next = NULL;
        set_to_add->last = NULL;
    }else{  // first element is already setted
        weighted_position_start_set->last = set_to_add;
        set_to_add->next = weighted_position_start_set;
    }
    this->ordered_by_score_set[weight] = set_to_add;
    return set_to_add;
}

void set_cover::add_set(const int set_id, const int set_weight,
                        const unsigned int * element_ids,
                        const unsigned short * weights,
                        const int element_size){
    
    set_::element * element_last_ptr = NULL;
    set_::element * element_first_ptr = NULL;
    set_* curr_set = &this->sets[set_id];
    curr_set->next = NULL;
    curr_set->last = NULL;
    curr_set->elements = NULL;		
    curr_set =create_set_at_weight_position(set_weight,curr_set);
    curr_set->set_id = set_id;
    curr_set->weight = set_weight;
    
    // set up doubled linked list + fill element_lookup
    for(int i = 0 ; i < element_size; i++) {
        // init elemnt with id, weight information
        set_::element * curr_element_ptr=&set_elements[add_position+i];
        curr_element_ptr->element_id=element_ids[i];
        curr_element_ptr->weight=weights[i];
        if(element_first_ptr == NULL) // first ptr is not yet set
            element_first_ptr = curr_element_ptr;
        // navigation (double linked list, parent)
        curr_element_ptr->parent_set = curr_set;
        curr_element_ptr->last = element_last_ptr;
        if(element_last_ptr != NULL) // not in first iteration
            element_last_ptr->next = curr_element_ptr;
        element_last_ptr = curr_element_ptr;
        // element_lookup fill up
        int element_id=curr_element_ptr->element_id;
        element_lookup->add_value_at(element_id,curr_element_ptr);
    }
    add_position += element_size;
    if(element_last_ptr!=NULL)
        element_last_ptr->next = NULL; // last element points to NULL
    curr_set->elements = element_first_ptr; // set element pointer to current_set

}



void set_cover::removeSet(set_* s){
    set_::element * element=s->elements;
    int s_set_id = s->set_id;
    unplug_set(s);
    do{ // for(Element element in elements
        set_::element * element_to_remove;
        int element_id=element->element_id;
        
        std::pair<set_::element**,int> element_lookup_structur=element_lookup->get_array(element_id);
        set_::element ** element_lookup_array = element_lookup_structur.first;
        int array_size = element_lookup_structur.second;
        for(int i =0; i < array_size;i++){
            element_to_remove = element_lookup_array[i];
            set_* parent_set = element_to_remove->parent_set;
            if(parent_set!=NULL && parent_set->set_id!=s_set_id){
                unplug_set(parent_set);
                parent_set->weight  -= element_to_remove->weight;
                parent_set->elements = unplug_element(element_to_remove,parent_set->elements);
                create_set_at_weight_position(parent_set->weight,parent_set);

            }
            //delete element_to_remove;
        }

    }while((element=element->next) != NULL);
}

set_::element * set_cover::unplug_element(set_::element * element_to_unplug, set_::element * first_element) {
    set_::element * last_element=element_to_unplug->last;
    set_::element * next_element=element_to_unplug->next;
    element_to_unplug->last       = NULL;
    element_to_unplug->next       = NULL;
    element_to_unplug->parent_set = NULL;

    if(last_element == NULL && next_element==NULL){
        return NULL;
    }if(last_element == NULL){ // first element
        next_element->last = NULL;
        return next_element;
    } else if (next_element == NULL) { // end of list
        last_element->next = NULL;
    } else { // middle of the list
        last_element->next = next_element;
        next_element->last = last_element;
    }
    return first_element;
}


void set_cover::unplug_set(set_* set_to_unplug){
    set_* last_set = set_to_unplug->last;
    set_* next_set = set_to_unplug->next;
    set_to_unplug->next = NULL;
    set_to_unplug->last = NULL;
    if(last_set == NULL && next_set == NULL){ // only one element left
        ordered_by_score_set[set_to_unplug->weight] = NULL;
    }else if(last_set == NULL){ // first set
        next_set->last = NULL;
        ordered_by_score_set[set_to_unplug->weight] = next_set;
    } else if (next_set == NULL) { // last set
        last_set->next = NULL;
    } else { // set in the middle
        last_set->next = next_set;
        next_set->last = last_set;
    }
}



std::list<set_*> set_cover::execute_set_cover(){
    //covered_vertex is a bit set which holds the covered states for all vertexes
    std::list<set_*> result;
    for(int i = weight_range; i > 0;) {
        while(this->ordered_by_score_set[i] == NULL && i > 0){
            i--;
        }
        if(i==0){ // if i == 0 than score = 0 -> no element in set anymore
            break;
        }
        set_* max_set = this->ordered_by_score_set[i];
        
        removeSet(max_set);
        result.push_back(max_set); // O(1)
    }
    return result;
}
