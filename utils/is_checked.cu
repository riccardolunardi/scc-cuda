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

// Alcune di queste funzioni verranno passate come parametri. Per permettere questo, in CUDA,
// vanno creati dei tipi appositi per le funzioni, che verranno utilizzati per inizializzare
// delle variabili che punteranno al blocco della funzione
typedef bool (*get_status)(char * value);
typedef void (*set_status)(char * value);

__host__ __device__ bool get_is_d_fw_visited(char * value)  { return *value & 1; }
__host__ __device__ bool get_is_d_bw_visited(char * value)  { return *value & 2; }
__host__ __device__ bool get_is_d_eliminated(char * value)  { return *value & 4; }
__host__ __device__ bool get_is_d_fw_expanded(char * value) { return *value & 8; }
__host__ __device__ bool get_is_d_bw_expanded(char * value) { return *value & 16; }
__host__ __device__ bool get_is_d_u(char * value)           { return *value & 32; }
__host__ __device__ bool get_is_d_scc(char * value)         { return *value & 64; }

__host__ __device__ void set_is_d_fw_visited(char * value)      { *value |= 1; }
__host__ __device__ void set_is_d_bw_visited(char * value)      { *value |= 2; }
__host__ __device__ void set_is_d_bw_fw_visited(char * value)   { *value |= 3; }
__host__ __device__ void set_is_d_eliminated(char * value)      { *value |= 4; }
__host__ __device__ void set_is_d_fw_expanded(char * value)     { *value |= 8; }
__host__ __device__ void set_is_d_bw_expanded(char * value)     { *value |= 16; }
__host__ __device__ void set_is_d_u(char * value)               { *value |= 32; }
__host__ __device__ void set_is_d_scc(char * value)             { *value |= 64; }

__host__ __device__ void set_not_is_d_fw_visited(char * value)      { *value &= 254; }
__host__ __device__ void set_not_is_d_bw_visited(char * value)      { *value &= 253; }
__host__ __device__ void set_not_is_d_eliminated(char * value)      { *value &= 251; }
__host__ __device__ void set_not_is_d_fw_expanded(char * value)     { *value &= 247; }
__host__ __device__ void set_not_is_d_bw_expanded(char * value)     { *value &= 239; }
__host__ __device__ void set_not_is_d_u(char * value)               { *value &= 223; }
__host__ __device__ void set_not_is_d_scc(char & value)             { *value &= 191; }

// Inizializzazione delle variabili che puntano alle funzioni.
// Da qeusti puntatori l'API riuscirà a risalire al codice della funzione
__device__ const get_status dev_get_fw_visited = get_is_d_fw_visited;
__device__ const get_status dev_get_bw_visited = get_is_d_bw_visited;
__device__ const get_status dev_get_fw_expanded = get_is_d_fw_expanded;
__device__ const get_status dev_get_bw_expanded = get_is_d_bw_expanded;

__device__ const set_status dev_set_fw_visited = set_is_d_fw_visited;
__device__ const set_status dev_set_bw_visited = set_is_d_bw_visited;
__device__ const set_status dev_set_fw_expanded = set_is_d_fw_expanded;
__device__ const set_status dev_set_bw_expanded = set_is_d_bw_expanded;