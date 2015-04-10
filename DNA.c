#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

/*Implementation of the Needleman-Wunch algorithm for finding the DNA alignment that create the highest score using this scoring system:
Match: +3
Mismatch: -3
Indel start: -8
Indel exstension: -1

The program only print the highest score
TODO: add the directions to backtrack the solutions
TODO: print the 2 sequences of DNA aligned.
*/
#define MAX_LEN 400

const short ind_start = -8;
const short ind_ext = -1;

typedef struct matrix_t {
	int len;
}Matr;

int max_val(int,...);
void print(Matr *m,int r,int c){
int i,j;
for (i=0;i<r;i++){
	for(j=0;j<c;j++)
		printf(" %d ",m[j+i*c].len);
 printf("\n");
}
};
int max_val(int n_args,...){
	Matr next;
	register int i;
	int max,a;
	va_list ap;
	va_start(ap,n_args);
	max = va_arg(ap,int);
	for(i=2;i<=n_args;i++){
		if((a =	va_arg(ap,int)) > max)
			max=a;
	}	
	va_end(ap);	
	next.len = max;
	return next.len;
};

int isequal(char a,char b)
{
	if(a == b)
		return 3;
	else
		return -3;
}
//Needlman-Wunsch algorithm based on the scoring system defined above
Matr *DNA_Align(char *A,char *B){

int lenA = strlen(A);
int lenB = strlen(B);

Matr *m = (Matr *) malloc ((lenA)*(lenB)*sizeof(Matr));
Matr *D = (Matr *) malloc ((lenA)*(lenB)*sizeof(Matr));
Matr *H = (Matr *) malloc ((lenA)*(lenB)*sizeof(Matr));
Matr *V = (Matr *) malloc ((lenA)*(lenB)*sizeof(Matr));


int i;
m[0].len = 0;
for(i = 1; i < lenB;i++){
	m[i].len = V[i].len = ind_start + (i-1)*ind_ext;
}
int j;
for (j = 1; j < lenA;j++){
	m[j*lenB].len = H[j*lenB].len = ind_start + (j-1)*ind_ext;
}

for(i = 1; i < lenA;i++){
    for (j = 1; j < lenB;j++){

		D[j+i*lenB].len = m[(j-1)+(i-1)*lenB].len + isequal(A[i],B[j]);
		H[j+i*lenB].len = max_val(2,H[(j-1)+i*lenB].len + ind_ext,m[(j-1)+i*lenB].len + ind_start); 
		V[j+i*lenB].len = max_val(2,V[j+(i-1)*lenB].len + ind_ext,m[j+(i-1)*lenB].len + ind_start);
		m[j + i*lenB].len = max_val(3,D[j+i*lenB].len,V[j+i*lenB].len,H[j+i*lenB].len);		
	
	}
	
}
						
//print(m,lenA,lenB);
printf("%d\n",m[(j-1)+(i-1)*lenB]);
return m;
};

int main(int argc,char **argv){

FILE *file = fopen(argv[1], "r");
char line[1024];
char *A;
char *B;
int i,endA,endB;
Matr *dna;
while (fgets(line, 1024, file)) {
	if (line[0] == '\n')
		continue;
	else {
		endA = 0;
		while (line[endA] != ' ')endA++;
		A = (char *) malloc ((endA+1) *sizeof(char));
		for (i = 0; i < endA;i++){
			if(line[i] != ' ')
				A[i+1] = line[i];
		}
		A[i+1] = '\0';
		A[0] = '0';
		endB = 0;
		while (line[endB+endA+3] != '\n') endB++;
		B = (char *) malloc ((endB+1) *sizeof(char));
		for (i = 0;i < endB;i++)
				B[i+1] = line [endA+3 + i];
		B[i+1] = '\0';
		B[0] = '0';
		if(strlen(A) != MAX_LEN && strlen(B) != MAX_LEN)
 			dna = DNA_Align(B,A);
		}
}
return 0;
}
