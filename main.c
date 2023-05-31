#include <math.h>
#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#define MAX_TEXT_SIZE 255

bool is_string_valid (char *);
bool is_consonant (char );
void print_matrix (gsl_matrix *, char *);
void encryption (char *);
void decryption (int *, int);

// declaring global variables
int n = 0;
int s = 0;
gsl_matrix * A;
gsl_matrix * B;
gsl_matrix * N;
gsl_matrix * K;
gsl_matrix * Q;
gsl_matrix * C1;
gsl_matrix * C2;

int main()
{
    char plain_text[MAX_TEXT_SIZE];

    // Gets the plain text
    printf("Type the message to be ciphered\n");
    fgets(plain_text, MAX_TEXT_SIZE, stdin);

    // Check if the plain text is valid
    if (!is_string_valid(plain_text)){
        printf("The plain text only allows letters A through Z, space and dot\n");
        return 1;
    }

    char non_repeated_str[29];
    int repeated_str[504];

    int nr_size = 0;
    int  r_size = 0;
    for ( ;n < MAX_TEXT_SIZE && plain_text[n] != '\n'; n++){
        bool is_repeated = false;
        for(int j = 0; j < nr_size; j++){
            if(non_repeated_str[j] == plain_text[n]){
                is_repeated = true;
                break;
            }
        }

        if (is_repeated){
            repeated_str[r_size] = n;
            repeated_str[r_size +1] = plain_text[n] +n;
            r_size += 2;
        } else {
            non_repeated_str[nr_size] = plain_text[n];
            nr_size++;
        }
    }
    non_repeated_str[nr_size] = '\n';
    n = nr_size;

    encryption(non_repeated_str);

    gsl_matrix_set_all(C1, 0);
    gsl_matrix_set_all(N, 0);
    gsl_matrix_set_all(B, 0);

    decryption(repeated_str, r_size/2);

    gsl_matrix_free(A);
    gsl_matrix_free(B);
    gsl_matrix_free(N);
    gsl_matrix_free(K);
    gsl_matrix_free(C1);
    gsl_matrix_free(Q);
    gsl_matrix_free(C2);
    return 0;
}


bool is_string_valid (char * plain_text){

    for (char *sub_string = plain_text; sub_string[0] != '\n'; sub_string++){
        if(!isalpha(sub_string[0]) && sub_string[0] != ' ' && sub_string[0] != '.'){
            return false;
        }

        // Turn to upperCase every letter
        sub_string[0] = toupper(sub_string[0]);
    }
    return true;
}


bool is_consonant (char letter){

    for (char* vowels= "AEIOU"; vowels[0] != '\0'; vowels++)
        if(vowels[0] == letter)
            return false;

    return true;
}

void print_matrix (gsl_matrix * matrix, char * title){

    printf("\n\n %s\n", title);

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            printf("%g ", gsl_matrix_get(matrix, i, j));
        }
        printf("\n");
    }
}


void encryption (char plain_text[MAX_TEXT_SIZE]){

    // Check if the message has 3 distinc characters minimum
    if (n < 3){
        printf("The plain text should have at least 3 distinc characters\n");
        exit(2);
    }

    // Translate the plain_text to an string of numbers
    int coded_plain_text[MAX_TEXT_SIZE];
    for(int i = 0; i < n; i++){

        char letter = plain_text[i];

        // Handle corner cases when the character is space or an dot
        if(letter == ' '){
            coded_plain_text[i] = 105;
            continue;
        } else if(letter == '.'){
            coded_plain_text[i] = 106;
            continue;
        }

        // Get the position of the letter in the alphabet
        int letter_num = letter- 64;

        // Calculate the row and column of the letter according to the especified table
        int row = ceil((double)letter_num /7) +6;
        int column = letter_num -(row -7)*7 -1;

        // Set the final value of the translation [A-Z] -> number
        if (is_consonant(letter)){
            coded_plain_text[i] = column*10 + row;
        } else{
            coded_plain_text[i] = row*10 + column;
        }
    }

    // Generate matrix A
    A = gsl_matrix_alloc(n, n);

    for(int i = 0; i < n; i++){
        for(int j = 0; j < i; j++){
            int value = abs(coded_plain_text[i]-coded_plain_text[j]);
            gsl_matrix_set(A, i, j, value);
            gsl_matrix_set(A, j, i, value);
        }
    }

    print_matrix(A, "matriz A >>>");

    // Calculate the matrix B
    B = gsl_matrix_alloc(n, n);

    for(int i = 0; i < n; i++){
        int y = i-1;
        int x = i+1;

        if (y == -1){
            y = n -1;
        } else if ( x == n ){
            x = 0;
        }

        gsl_matrix_set(B, i, y, gsl_matrix_get(A, i, y));
        gsl_matrix_set(B, i, i, plain_text[i]);
        gsl_matrix_set(B, i, x, gsl_matrix_get(A, i, x));
    }

    print_matrix(B, "matriz B >>>");

    // Calculate the matrix N
    N = gsl_matrix_alloc(n, n);
    gsl_blas_dsymm(CblasLeft, CblasUpper, 1.0, A, B, 0.0, N);

    print_matrix(N, "N =>>>");

    // Calculate the key matrix K
    K = gsl_matrix_alloc(n, n);

    for (int i = 0; i < n; i++){
        for (int j = 0, max = n -i; j < max; j++){
            gsl_matrix_set(K, i, i+j, j+1);
        }
    }

    print_matrix(K, "matriz K >>>");

    // Calculate the matrix C1
    C1 = gsl_matrix_alloc(n, n);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, N, K, 0.0, C1);


    print_matrix(C1, "matriz C1 >>>");

    // Calculate S
    int hamilton_cycle[n], index = 0;
    hamilton_cycle[0] = 0;

    while (index < n-1){

        int vertexes_weight[n];
        int idx_vertex_lowest_weight = -1;

        // Calculate the weight of the direct transitive closure from the last vertex in hamilton_cycle
        for (int i = 0; i < n; i++){

            // The matrix diagonal is discarted
            if (hamilton_cycle[index] == i){
                vertexes_weight[i] = -1;
                continue;
            }

            vertexes_weight[i] = gsl_matrix_get(A, hamilton_cycle[index],i);
            if (idx_vertex_lowest_weight == -1 || vertexes_weight[idx_vertex_lowest_weight] > vertexes_weight[i]){
                idx_vertex_lowest_weight = i;
            }
        }

        // Ensure the vertex is not repeated in the hamilton_cycle
        while(true){
            bool is_repeated = false;

            for (int i = 0; i <= index; i++){
                if (hamilton_cycle[i] == idx_vertex_lowest_weight ){
                    is_repeated = true;
                    break;
                }
            }

            if (!is_repeated)
                break;

            // Discards idx_vertex_lowest_weight in vertexes_weight because it is repeated in hamilton_cycle
            vertexes_weight[idx_vertex_lowest_weight] = -1;
            idx_vertex_lowest_weight = -1;

            // Finds the next vertex with the lowest weight
            for (int i = 0; i < n; i++)
                if ((idx_vertex_lowest_weight == -1 || vertexes_weight[idx_vertex_lowest_weight] > vertexes_weight[i]) && vertexes_weight[i] != -1)
                    idx_vertex_lowest_weight = i;
        }

        s += vertexes_weight[idx_vertex_lowest_weight];         // Sum the vertex add in the hamilton cycle
        hamilton_cycle[++index] = idx_vertex_lowest_weight;    // add the Vertex to the hamilton cycle
    }

    s += gsl_matrix_get(A, 0, hamilton_cycle[index]);       // sum the vertex of the first and last vertex to calculate the final value of s
    printf("\nvalor de s: %d\n",s);

    // Calculate C2 and Q
    C2 = gsl_matrix_alloc(n, n);
    Q  = gsl_matrix_alloc(n, n);

    printf("\nCifra final em forma linear: ");


    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            int C1_value = gsl_matrix_get(C1, i, j);
            gsl_matrix_set(C2, i, j, C1_value % s);

            printf("%g", gsl_matrix_get(C2, i, j));         // Final cypher printed in linear form

            gsl_matrix_set( Q, i, j, floor(C1_value/s));
        }
    }

    print_matrix(C2, "matriz C2 >>>");
    print_matrix(Q, "matriz Q >>>");
}


void decryption (int * repeated_letters, int size){

    // Obtein matrix C1
    C1 = gsl_matrix_alloc(n, n);

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            int remainder = gsl_matrix_get(C2, i, j);
            int quocient = gsl_matrix_get(Q, i, j);
            gsl_matrix_set(C1, i, j, quocient *s + remainder);
        }
    }

    print_matrix(C1, "matrix C1 reconstruída >>>");

    // Calculate the inverse matrix of K
    gsl_matrix * inverse_K = gsl_matrix_alloc(n, n);
    for(int i = 0; i < n; i++){
        gsl_matrix_set(inverse_K, i, i, 1);
        if (i+1 < n){
            gsl_matrix_set(inverse_K, i, i+1, -2);
        }
        if (i+2 < n){
            gsl_matrix_set(inverse_K, i, i+2, 1);
        }
    }


    print_matrix(inverse_K, "inverso da matriz K >>>");

    // Calculate the matrix N
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, C1, inverse_K, 0.0, N);

    print_matrix(N, "matriz N reconstruída");

    // Calculate the inverse matrix of A
    gsl_matrix * inverse_A = gsl_matrix_alloc(n, n);

    int signum;
    gsl_permutation * p = gsl_permutation_alloc(n);

    gsl_linalg_LU_decomp (A, p, &signum);
    gsl_linalg_LU_invert (A, p, inverse_A);

    gsl_permutation_free(p);

    print_matrix(inverse_A, "inverso da matriz A >>>");


    // Calculate the matrix B
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, inverse_A, N, 0.0, B);

    print_matrix(B, "matriz B reconstruída");

    printf("\nMensagem decodificada: ");

    int count_repeated_letters = 0;
    // Print the decoded message
    for(int i = 0, max = n + size; i < max; i++){

        int idx_value = repeated_letters[count_repeated_letters*2];
        if (i == idx_value){
            printf("%c", repeated_letters[(count_repeated_letters*2)+1] -idx_value);
            count_repeated_letters++;
        }
        else
            printf("%c", (int)round(gsl_matrix_get(B, i-count_repeated_letters, i-count_repeated_letters)));
    }

    printf("\n\n\n");

    gsl_matrix_free(inverse_K);
    gsl_matrix_free(inverse_A);
}

