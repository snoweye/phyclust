/* This file contains functions to read input file. */

#include <stdlib.h>
#include <stdio.h>
#include "phyclust_constant.h"
#include "phyclust_file_input.h"
#include "phyclust_tool.h"


input_struct* initialize_input_struct(int code_type, int N_X_org, int L){
	input_struct *ins;

	ins = (input_struct*) malloc(sizeof(input_struct));
	ins->code_type = code_type;
	ins->ncode = NCODE[code_type];
	ins->N_X_org = N_X_org;
	ins->L = L;
	ins->X_org = allocate_int_RT(N_X_org, L);
	ins->X_name = allocate_char_RT(N_X_org, NAME_LENGTH);

	return(ins);
} /* End of initialize_input_struct(). */


void free_input_struct(input_struct *ins){
	free_int_RT(ins->N_X_org, ins->X_org);
	free_char_RT(ins->N_X_org, ins->X_name);
	free(ins);
} /* End of free_input_struct(). */




int nucleotide_to_id(char x){
	int ret;

	switch(x){
		case 'A': case 'a':
			ret = A;
			break;
		case 'G': case 'g':
			ret = G;
			break;
		case 'C': case 'c':
			ret = C;
			break;
		case 'T': case 't': case 'U': case 'u':
			ret = T;
			break;
		case '-':
			ret = GAP;
			break;
		default:
			ret = MISSING_ALLELE;
			break;
	}
	return(ret);
} /* End of nucleotide_to_id(). */


int is_nucleotide(char X_org){
	if(X_org == 'A' || X_org == 'G' || X_org == 'C' || X_org == 'T' || X_org == 'U' ||
		X_org == 'a' || X_org == 'g' || X_org == 'c' || X_org == 't' || X_org == 'u' ||
		X_org == '-' || X_org == '.' ||
		X_org == 'K' || X_org == 'M' || X_org == 'N' || X_org == 'R' || X_org == 'W' || X_org == 'Y' ||
		X_org == 'k' || X_org == 'm' || X_org == 'n' || X_org == 'r' || X_org == 'w' || X_org == 'y'
		){
		return(1);
	} else{
		return(0);
	}
} /* End of is_nucleotide(). */


/* Read PHYLIP file. */
input_struct* read_input_phylip(char *file_name){
	input_struct *ins = NULL;
	FILE *fp;
	char X_org;
	int i, j=0, N_X_org, L, count_L;

	fp = fopen(file_name, "r");

	if(fp == NULL){
		fprintf_stderr("PE: can't open file \"%s\".\n", file_name);
		exit(1);
	} else{
		if(!fscanf(fp, "%d %d", &N_X_org, &L)) {
			fprintf_stderr("PE: invalid PHYLIP format in file \"%s\".\n", file_name);
			exit(1);
		}
		printf("Read PHYLIP(%s): N_X_org=%d L=%d code_type=%s\n", file_name, N_X_org, L, CODE_TYPE[NUCLEOTIDE]);
		do{
			X_org = (char) fgetc(fp);
		} while(X_org != '\n');

		ins = initialize_input_struct(NUCLEOTIDE, N_X_org, L);
		count_L = 0;

		for(i = 0; i < N_X_org; i++){
			for(j = 0; j < NAME_LENGTH; j++){
				ins->X_name[i][j] = (char) fgetc(fp);
			}

			j = count_L;
			do{
				X_org = (char) fgetc(fp);
				if(is_nucleotide(X_org)){
					ins->X_org[i][j++] = nucleotide_to_id(X_org);
				}
			} while(X_org != '\n');
		}
		count_L = j;

		while(count_L < L){
			do{
				X_org = (char) fgetc(fp);	/* Skip a new line. */
			} while(X_org != '\n');
			for(i = 0; i < N_X_org; i++){
				j = count_L;
				do{
					X_org = (char) fgetc(fp);
					if(is_nucleotide(X_org)){
						ins->X_org[i][j++] = nucleotide_to_id(X_org);
					}
				} while(X_org != '\n');
			}
			count_L = j;
		}
	}

	fclose(fp);

	return(ins);
} /* End of read_input_phylip(). */


/* Read FASTA file. */
input_struct* read_input_fasta(char *file_name){
	input_struct *ins = NULL;
	FILE *fp;
	char X_org;
	int i, j, N_X_org, L;

	fp = fopen(file_name, "r");

	if(fp == NULL){
		fprintf_stderr("PE: can't open file \"%s\".\n", file_name);
		exit(1);
	} else{
		/* Count total sequences. */
		N_X_org = 0;
		L = 0;
		while(! feof(fp)){
			X_org = (char) fgetc(fp);
			if(X_org == '\r'){
				continue;
			}
			if(X_org == '>'){
				N_X_org++;
				while(! feof(fp)){
					X_org = (char) fgetc(fp);
					if(X_org == '\r'){
						continue;
					}
					if(X_org == '\n'){
						break;
					}
				}
				continue;
			}
			if(N_X_org == 1 && X_org != '\n'){
				L++;
			}
		}

		fseek(fp, 0, SEEK_SET);
		printf("Read FASTA(%s): N_X_org=%d L=%d code_type=%s\n", file_name, N_X_org, L, CODE_TYPE[NUCLEOTIDE]);

		ins = initialize_input_struct(NUCLEOTIDE, N_X_org, L);

		/* Read sequences. */
		i = -1;
		j = 0;
		while(! feof(fp)){
			X_org = (char) fgetc(fp);
			if(X_org == '\r'){
				continue;
			}
			if(X_org == '>'){
				j = 0;
				while(! feof(fp)){
					X_org = (char) fgetc(fp);
					if(X_org == '\r'){
						continue;
					}
					if(X_org == '\n'){
						break;
					}
					if(j < NAME_LENGTH){
						ins->X_name[i + 1][j++] = X_org;
					}
				}
				ins->X_name[i + 1][j] = '\0';
				i++;
				j = 0;
				continue;
			}
		       
			if(is_nucleotide(X_org)){
				ins->X_org[i][j++] = nucleotide_to_id(X_org);
			}
		}
	}

	fclose(fp);

	return(ins);
} /* End of read_input_fasta(). */




int snp_to_id(char x){
	int ret;

	switch(x){
		case '1':
			ret = NL;
			break;
		case '2':
			ret = SU;
			break;
		case '-':
			ret = UN;
			break;
		default:
			ret = MISSING_ALLELE;
			break;
	}
	return(ret);
} /* End of snp_to_id(). */


int is_snp(char X_org){
	if(X_org == '1' || X_org == '2' || X_org == '-' || X_org == '.'){
		return(1);
	} else{
		return(0);
	}
} /* End of is_snp(). */


/* Read SNP file. The same format as PHYLIP, but only for SNP. */
input_struct* read_input_snp(char *file_name){
	input_struct *ins = NULL;
	FILE *fp;
	char X_org;
	int i, j=0, N_X_org, L, count_L;

	fp = fopen(file_name, "r");

	if(fp == NULL){
		fprintf_stderr("PE: can't open file \"%s\".\n", file_name);
		exit(1);
	} else{
		if(!fscanf(fp, "%d %d", &N_X_org, &L)) {
			fprintf_stderr("PE: invalid PHYLIP format in file \"%s\".\n", file_name);
			exit(1);
		}
		printf("Read SNP(%s): N_X_org=%d L=%d code_type=%s\n", file_name, N_X_org, L, CODE_TYPE[SNP]);
		do{
			X_org = (char) fgetc(fp);
		} while(X_org != '\n');

		ins = initialize_input_struct(SNP, N_X_org, L);
		count_L = 0;

		for(i = 0; i < N_X_org; i++){
			for(j = 0; j < NAME_LENGTH; j++){
				ins->X_name[i][j] = (char) fgetc(fp);
			}

			j = count_L;
			do{
				X_org = (char) fgetc(fp);
				if(is_snp(X_org)){
					ins->X_org[i][j++] = snp_to_id(X_org);
				}
			} while(X_org != '\n');
		}
		count_L = j;

		while(count_L < L){
			do{
				X_org = (char) fgetc(fp);	/* Skip a new line. */
			} while(X_org != '\n');
			for(i = 0; i < N_X_org; i++){
				j = count_L;
				do{
					X_org = (char) fgetc(fp);
					if(is_snp(X_org)){
						ins->X_org[i][j++] = snp_to_id(X_org);
					}
				} while(X_org != '\n');
			}
			count_L = j;
		}
	}

	fclose(fp);

	return(ins);
} /* End of read_input_snp(). */

