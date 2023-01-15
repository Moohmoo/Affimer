#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "file_management.h"

double threshold = 16;
char* path = "Patch_Search/Sites_Modifies/";

void get_residus_implicated() {
    return;
}

void evaluate_score(char* sub_name_target) {
    FILE* f_input;
    FILE* f_save;

    long length;
    char* buffer = 0;

    char* buffer_pdb = 0;

    char input_file[20];
    strcpy(input_file, sub_name_target);
    strcat(input_file, ".out");

    f_input = fopen(input_file, "rb");

    if (f_input) {
        fseek(f_input, 0, SEEK_END);
        length = ftell(f_input);
        fseek(f_input, 0, SEEK_SET);
        buffer = malloc(length * sizeof(char) + 1);
        if (buffer) {
            fread(buffer, 1, length, f_input);
        }
        fclose(f_input);
    }

    char** array_input_file;
	int size_input_file;

    string_split(buffer, '\n', &array_input_file, &size_input_file);

    int max_i = (size_input_file - 1);

    char** array;
    int size;

    char* score_char;
    double score;

    char tmp[80];

    char sub_name[20];
    char chain[2];

    char** array_pdb;
    int size_pdb;

    int max_j;

    char residu[20];
    char* residus_infos = malloc(sizeof(char) * 100000000);

    char* write = malloc(sizeof(char) * 100000000);
    strcpy(write, "Score\n");

    for (int i = 1; i < max_i; i++) {
        string_split(array_input_file[i], ' ', &array, &size);
        score_char = array[size-1];
        score = fast_atof(score_char);
        if (score != 0.0) {
            strcat(write, score_char);
            strcat(write, "\n");
        }
        
        if (score > threshold) {
            strcpy(tmp, array[0] + strlen(path) + 1);
            strncpy(sub_name, tmp, strlen(tmp) - 10);
            sub_name[strlen(tmp) - 10] = '\0';
            chain[0] = sub_name[strlen(sub_name) - 3];
            chain[1] = '\0';
            
            buffer_pdb = read_binding_sites(sub_name);
            array_pdb = 0;
            string_split(buffer_pdb, '\n', &array_pdb, &size_pdb);

            max_j = (size_pdb - 1);
            strcpy(residus_infos, "");
            for(int j = 0; j < max_j; j++) {
                strcpy(tmp, array_pdb[j] + 22);
                strncpy(residu, tmp, strlen(tmp) - 50);
                residu[strlen(tmp) - 50] = '\0';
                if (array_pdb[j][21] == chain[0] && strstr(residus_infos, residu) == NULL) {
                    strcat(residus_infos, residu);
                    strcat(residus_infos, "\n");
                }
            }
            printf("%s\n", residus_infos);
        }
    }

    free(residus_infos);

    f_save = fopen("score.csv", "w");
    if (f_save) {
        fprintf(f_save, "%s", write);
        fclose(f_save);    
    }
    printf("Done !\n");
}


int main(int argc, char** argv) 
{
    evaluate_score("5GS6");
    return 0;
}