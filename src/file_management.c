#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "file_management.h"

#define white_space(c) ((c) == ' ' || (c) == '\t')
#define valid_digit(c) ((c) >= '0' && (c) <= '9')

double fast_atof (const char *p)
{
    int frac;
    double sign, value, scale;

    // Skip leading white space, if any.

    while (white_space(*p) ) {
        p += 1;
    }

    // Get sign, if any.

    sign = 1.0;
    if (*p == '-') {
        sign = -1.0;
        p += 1;

    } else if (*p == '+') {
        p += 1;
    }

    // Get digits before decimal point or exponent, if any.

    for (value = 0.0; valid_digit(*p); p += 1) {
        value = value * 10.0 + (*p - '0');
    }

    // Get digits after decimal point, if any.

    if (*p == '.') {
        double pow10 = 10.0;
        p += 1;
        while (valid_digit(*p)) {
            value += (*p - '0') / pow10;
            pow10 *= 10.0;
            p += 1;
        }
    }

    // Handle exponent, if any.

    frac = 0;
    scale = 1.0;
    if ((*p == 'e') || (*p == 'E')) {
        unsigned int expon;

        // Get sign of exponent, if any.

        p += 1;
        if (*p == '-') {
            frac = 1;
            p += 1;

        } else if (*p == '+') {
            p += 1;
        }

        // Get digits of exponent, if any.

        for (expon = 0; valid_digit(*p); p += 1) {
            expon = expon * 10 + (*p - '0');
        }
        if (expon > 308) expon = 308;

        // Calculate scaling factor.

        while (expon >= 50) { scale *= 1E50; expon -= 50; }
        while (expon >=  8) { scale *= 1E8;  expon -=  8; }
        while (expon >   0) { scale *= 10.0; expon -=  1; }
    }

    // Return signed and scaled floating point result.

    return sign * (frac ? (value / scale) : (value * scale));
}


char* read_chains(char* sub_name, char* chain) {
    FILE* f_pdb;
    long length;
    char* buffer = 0;

    char path_file[100];
    char* path = "../data/chains/";

    strcpy(path_file, path);
    strcat(path_file, sub_name);
    strcat(path_file, "_");
    strcat(path_file, chain);
    strcat(path_file, ".pdb");
    f_pdb = fopen(path_file, "rb");

    if (f_pdb) {
        fseek(f_pdb, 0, SEEK_END);
        length = ftell(f_pdb);
        fseek(f_pdb, 0, SEEK_SET);
        buffer = malloc(length * sizeof(char) + 1);
        if (buffer) {
            fread(buffer, 1, length, f_pdb);
        }
        fclose(f_pdb);
    }

    return buffer;
}

char* read_binding_sites(char* sub_name) {
    FILE* f_pdb;
    long length;
    char* buffer = 0;

    char path_file[100];
    char* path = "../data/binding_sites/";

    strcpy(path_file, path);
    strcat(path_file, sub_name);
    strcat(path_file, ".pdb");
    f_pdb = fopen(path_file, "rb");
    if (f_pdb) {
        fseek(f_pdb, 0, SEEK_END);
        length = ftell(f_pdb);
        fseek(f_pdb, 0, SEEK_SET);
        buffer = malloc(length * sizeof(char) + 1);
        if (buffer) {
            fread(buffer, 1, length, f_pdb);
        }
    }
    fclose(f_pdb);
    return buffer;
}

void create_sites(char* sub_name, char* chain_ab, char* chain_ag, char* write_ab, char* write_ag) {
    FILE* f_out;

    char path_file_out[100];
    char* path = "../data/binding_sites/";

    strcpy(path_file_out, path);
    strcat(path_file_out, sub_name);
    strcat(path_file_out, "_");
    strcat(path_file_out, chain_ag);
    strcat(path_file_out, "-");
    strcat(path_file_out, chain_ab);
    strcat(path_file_out, ".pdb");
    f_out = fopen(path_file_out, "w");

    if (f_out) {
        fprintf(f_out, "%s", write_ab);
        fprintf(f_out, "%s", write_ag);
        fclose(f_out);    
    }
}


void string_split(char *string, char sep, char ***r_array_string, int *r_size) {
    int i, k, len, size;
	char** array_string;
    int is_continuous = 0;
    char* tmp;

    size = 0, len = strlen(string);
	for(i = 0; i < len - 1; i++) {
		if(string[i] != sep && !is_continuous) {
			size++;
            is_continuous = 1;
		} else if (string[i] == sep) {
            is_continuous = 0;
        }
	}
    array_string = malloc(size * sizeof(char*));

    tmp = malloc(sizeof(char) * strlen(string) + 1);
    strcpy(tmp, string);

    i=0, k=0;
	array_string[k++] = tmp;

	while(k < size) {
		if(tmp[i++] == sep) {
            while(tmp[i] == sep)  {
                i++;
            }
			tmp[i-1] = '\0'; // Set end of substring
			array_string[k++] = (tmp+i); // Save the next substring pointer
		}
	}
	*r_array_string = array_string;
	*r_size = size;
	return;
}

/*void string_split(char *string, char* sep, char ***r_array_string, int *r_size) {
    char* tmp = malloc(sizeof(char) * strlen(string) + 1);
    strcpy(tmp, string);

    char** array_string = malloc(sizeof(char*));

    char * token = strtok(tmp, sep);
    int size = 0;

    while(token != NULL) {
        array_string = realloc(array_string, (size + 1) * sizeof(char*));
        array_string[size++] = token;
        token = strtok(NULL, sep);
    }

    *r_array_string = array_string;
	*r_size = size;
}*/

