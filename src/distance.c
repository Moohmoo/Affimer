#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <omp.h>
#include "file_management.h"
#include "distance.h"

double threshold;

char* get_chains(FILE* ptr, char* file_name) {
    char* chains = malloc(sizeof(char) * 100);
    strcpy(chains, "HL");

    int bufferLength = 255;
    char buffer[bufferLength];

    FILE* f_light;
    FILE* f_heavy;
    FILE* f_chain = NULL;
    FILE* f_bis = NULL;
    FILE* f_ter = NULL;

    char* path = "../data/chains/";
    char sub_name[(strlen(file_name) - 4) + 1];
    memcpy(sub_name, &file_name[0], strlen(file_name) - 4);
    sub_name[strlen(file_name) - 4] = '\0';

    char path_file_light[100];
    strcpy(path_file_light, path);
    strcat(path_file_light, sub_name);
    strcat(path_file_light, "_L.pdb");
    f_light = fopen(path_file_light, "w+");

    char path_file_heavy[100];
    strcpy(path_file_heavy, path);
    strcat(path_file_heavy, sub_name);
    strcat(path_file_heavy, "_H.pdb");
    f_heavy = fopen(path_file_heavy, "w+");

    int flag_chain = 0;
    int flag_bis = 0;
    int flag_ter = 0;

    while(fgets(buffer, bufferLength, ptr)) {
        if (strstr(buffer, "ATOM") != NULL) {
            char chain_name[2];
            memcpy(chain_name, &buffer[21], 1);
            chain_name[1] = '\0';
            if (strcmp(chain_name, "H") == 0) {
                fprintf(f_heavy, "%s", buffer);
            } else if (strcmp(chain_name, "L") == 0) {
                fprintf(f_light, "%s", buffer);
            } else {
                if (strstr(chains, chain_name) != NULL) {
                    if (buffer[21] == chains[2]) {
                        fprintf(f_chain, "%s", buffer);
                    } else if (buffer[21] == chains[3]) {
                        fprintf(f_bis, "%s", buffer);
                    } else {
                        fprintf(f_ter, "%s", buffer);
                    }
                } else if (flag_chain == 0) {
                    strcat(chains, chain_name);

                    char path_file_chain[100];
                    strcpy(path_file_chain, path);
                    strcat(path_file_chain, sub_name);
                    strcat(path_file_chain, "_");
                    strcat(path_file_chain, chain_name);
                    strcat(path_file_chain, ".pdb");
                    f_chain = fopen(path_file_chain, "w+");

                    fprintf(f_chain, "%s", buffer);
                    flag_chain = 1;
                } else if (flag_bis == 0) {
                    strcat(chains, chain_name);

                    char path_file_bis[100];
                    strcpy(path_file_bis, path);
                    strcat(path_file_bis, sub_name);
                    strcat(path_file_bis, "_");
                    strcat(path_file_bis, chain_name);
                    strcat(path_file_bis, ".pdb");
                    f_bis = fopen(path_file_bis, "w+");

                    fprintf(f_bis, "%s", buffer);
                    flag_bis = 1;
                } else if (flag_ter == 0) {
                    strcat(chains, chain_name);

                    char path_file_ter[100];
                    strcpy(path_file_ter, path);
                    strcat(path_file_ter, sub_name);
                    strcat(path_file_ter, "_");
                    strcat(path_file_ter, chain_name);
                    strcat(path_file_ter, ".pdb");
                    f_ter = fopen(path_file_ter, "w+");

                    fprintf(f_ter, "%s", buffer);
                    flag_ter = 1;
                }
            }
        }
    }

    fclose(f_heavy);
    fclose(f_light);

    if (f_chain != NULL)
        fclose(f_chain);
    if (f_bis != NULL)
        fclose(f_bis);
    if (f_ter != NULL)
        fclose(f_ter);

    return chains;
}

void get_residue(char* str, char** residue) {
    char* tmp = malloc(sizeof(char) * strlen(str) + 1);
    strcpy(tmp, str + 22);
    tmp = strtok(tmp, " ");
    *residue = tmp;
}

void get_coordinate(char* str, double** coordinate) {
    (*coordinate)[0] = fast_atof(str+30);
    (*coordinate)[1] = fast_atof(str+38);
    (*coordinate)[2] = fast_atof(str+46);
}

double euclidienne_distance(double* coordinates_antibody, double* coordinates_antigen) {
    double distance = 0;
    for (int j = 0; j < 3; j++) {
        distance += (coordinates_antigen[j] - coordinates_antibody[j]) 
        * (coordinates_antigen[j] - coordinates_antibody[j]);
    }
    distance = sqrt(distance);
    return distance;
}

void get_alpha_carbon(char** array, int size_array, char* residue, int first_pos, char** buffer) {
    char* current_residue;
    get_residue(array[first_pos], &current_residue);
    for (int i = first_pos; strcmp(current_residue, residue) == 0 && i < size_array - 1; i++) {
        if (strstr(array[i], "CA") != NULL) {
            *buffer = array[i];
            return;
        }
    }
    return;
}

void evaluate_distance(char* file_name, char* chains, int is_light) {
    char* buffer_antibody = 0;
    char* buffer_antigen = 0;

    char sub_name[(strlen(file_name) - 4) + 1];
    memcpy(sub_name, &file_name[0], strlen(file_name) - 4);
    sub_name[strlen(file_name) - 4] = '\0';
    char chain_name[2];

    buffer_antibody = read_chains(sub_name, (is_light)? "L" : "H");

    char** array_antibody;
	int size_antibody;
    char** array_antigen;
	int size_antigen;

    string_split(buffer_antibody, '\n', &array_antibody, &size_antibody);

    int max_j = (size_antibody - 1); 
    int max_l;

    double* coordinate_antibody = malloc(sizeof(double) * 3);
    double* coordinate_antigen = malloc(sizeof(double) * 3);

    char* residue_antibody = malloc(sizeof(char) * 10);
    char* residue_antigen = malloc(sizeof(char) * 10);

    char* current_residue_antibody = "1";
    char* current_residue_antigen = "1";

    int first_pos_antibody = 0;
    int first_pos_antigen = 0;

    double distance = 6.0;

    char* alpha_carbon_antibody;
    char* alpha_carbon_antigen;

    alpha_carbon_antibody = NULL;

    char* write_antibody = malloc(sizeof(char) * 100000000);
    char* write_antigen = malloc(sizeof(char) * 100000000);

    for (int i = 2; i < strlen(chains); i++) {
        strcpy(write_antibody, "");
        strcpy(write_antigen, "");

        memcpy(chain_name, &chains[i], 1);
        chain_name[1] = '\0';
        buffer_antigen = read_chains(sub_name, chain_name);

        string_split(buffer_antigen, '\n', &array_antigen, &size_antigen);

        alpha_carbon_antigen = NULL;
        max_l = (size_antigen - 1);

        //#pragma omp parallel for collapse(2)
        for (int j = 0; j < max_j; j++) {
            get_coordinate(array_antibody[j], &coordinate_antibody);
            get_residue(array_antibody[j], &residue_antibody);

            if (strcmp(current_residue_antibody, residue_antibody) != 0) {
                current_residue_antibody = residue_antibody;
                first_pos_antibody = j;
            }

            for (int l = 0; l < max_l; l++) {
                get_coordinate(array_antigen[l], &coordinate_antigen);
                get_residue(array_antigen[l], &residue_antigen);

                if (strcmp(current_residue_antigen, residue_antigen) != 0) {
                    current_residue_antigen = residue_antigen;
                    first_pos_antigen = l;
                }

                distance = euclidienne_distance(coordinate_antibody, coordinate_antigen);
                if (distance <= threshold) {
                    get_alpha_carbon(array_antibody, size_antibody, residue_antibody, first_pos_antibody, &alpha_carbon_antibody);
                    get_alpha_carbon(array_antigen, size_antigen, residue_antigen, first_pos_antigen, &alpha_carbon_antigen);

                    if (alpha_carbon_antibody != NULL && strstr(write_antibody, alpha_carbon_antibody) == NULL) {
                        strcat(write_antibody, alpha_carbon_antibody);
                        strcat(write_antibody, "\n");
                    }
                    if (strstr(write_antibody, array_antibody[j]) == NULL) {
                        strcat(write_antibody, array_antibody[j]);
                        strcat(write_antibody, "\n");
                    }

                    if (alpha_carbon_antigen != NULL && strstr(write_antigen, alpha_carbon_antigen) == NULL) {
                        strcat(write_antigen, alpha_carbon_antigen);
                        strcat(write_antigen, "\n");
                    }
                    if (strstr(write_antigen, array_antigen[l]) == NULL) {
                        strcat(write_antigen, array_antigen[l]);
                        strcat(write_antigen, "\n");
                    }
                }
            }
        }
        
        if (strlen(write_antibody) != 0 && strlen(write_antigen) != 0) {
            create_sites(sub_name, (is_light) ? "L" : "H", chain_name, write_antibody, write_antigen);
        }
    }
}

int main(int argc, char** argv)
{
    FILE* ptr;

    char* path = "../data/dataset/";

    if (argc == 4) {
        
        double input_threshold = fast_atof(argv[3]);

        if ((strcmp(argv[2], "H") == 0 || strcmp(argv[2], "L") == 0) && input_threshold != 0.0) {
            threshold = input_threshold;

            char path_file[200];
            strcpy(path_file, path);
            strcat(path_file, argv[1]);
            if (access(path_file, F_OK) == 0) {
                ptr = fopen(path_file, "r");
                if (ptr == NULL) {
                    fprintf(stderr, "[ERROR] File can't be opened.\n");
                    exit(1);
                } else {
                    char* chains = get_chains(ptr, argv[1]);
                    evaluate_distance(argv[1], chains, (strcmp(argv[2], "H") == 0) ? 0 : 1);
                }
            } else {
                fprintf(stderr, "[ERROR] This file does not exist or is not in the dedicated directory.\n");
                exit(1);
            }
        } else {
            fprintf(stderr, "[ERROR] Invalid arguments.\n");
            exit(1);
        }
    } else {
        fprintf(stderr, "[ERROR] Invalid arguments.\n");
        exit(1);
    }
    return 0;
}