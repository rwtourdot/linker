#include "read_tree.h"

void read_tree::set_name( std::string name ) { readname = name; };

void read_tree::add_connection( std::string hetname ) {
    connected_strings.push_back( hetname );
    num_cnx += 1;
};

void read_tree::add_connection_int( int index, std::string chr ) {
    connected_ints.push_back( index );
    connected_chr.push_back( chr );
    num_cnx += 1;
};

