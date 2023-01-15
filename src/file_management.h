#ifndef FILE_MANAGEMENT
#define FILE_MANAGEMENT

double fast_atof (const char *p);

char* read_chains(char* sub_name, char* chain);

char* read_binding_sites(char* sub_name);

void create_sites(char* sub_name, char* chain_ab, char* chain_ag, char* write_ab, char* write_ag);

void string_split(char *string, char sep, char ***r_array_string, int *r_size);

#endif