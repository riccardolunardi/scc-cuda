#ifndef IS_CHECKED_CPP
#define IS_CHECKED_CPP

#include <iostream>
using namespace std;
 
// Usiamo un vettore di char
// ogni char è composto da 8 bit così assegnati:
//  - is_fw_visited     1° bit da destra
//  - is_bw_visited     2° bit da destra
//  - is_eliminated     3° bit da destra
//  - is_fw_expanded    4° bit da destra
//  - is_bw_expanded    5° bit da destra
//  - is_u              6° bit da destra
//  - is_scc            7° bit da destra

bool get_is_fw_visited(char value)  { return value & 1; }
bool get_is_bw_visited(char value)  { return value & 2; }
bool get_is_eliminated(char value)  { return value & 4; }
bool get_is_fw_expanded(char value) { return value & 8; }
bool get_is_bw_expanded(char value) { return value & 16; }
bool get_is_u(char value)           { return value & 32; }
bool get_is_scc(char value)         { return value & 64; }

void set_is_fw_visited(char & value)    { value |= 1; }
void set_is_bw_visited(char & value)    { value |= 2; }
void set_is_eliminated(char & value)    { value |= 4; }
void set_is_fw_expanded(char & value)   { value |= 8; }
void set_is_bw_expanded(char & value)   { value |= 16; }
void set_is_u(char & value)             { value |= 32; }
void set_is_scc(char & value)           { value |= 64; }

void set_not_is_fw_visited(char & value)    { value &= 254; }
void set_not_is_bw_visited(char & value)    { value &= 253; }
void set_not_is_eliminated(char & value)    { value &= 251; }
void set_not_is_fw_expanded(char & value)   { value &= 247; }
void set_not_is_bw_expanded(char & value)   { value &= 239; }
void set_not_is_u(char & value)             { value &= 223; }
void set_not_is_scc(char & value)           { value &= 191; }

#endif