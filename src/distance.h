#ifndef DISTANCE
#define DISTANCE

double threshold;

char* get_chains(FILE* ptr, char* file_name);

void get_residue(char* str, char** residue);

void get_coordinate(char* str, double** coordinate);

double euclidienne_distance(double* coordinates_antibody, double* coordinates_antigen);

void get_alpha_carbon(char** array, int size_array, char* residue, int first_pos, char** buffer);

void evaluate_distance(char* file_name, char* chains, int is_light);

#endif