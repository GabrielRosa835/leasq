#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

typedef struct 
{
    float* coeficientes;
    int tamanho;
    float termo_independente;
} linha;

typedef struct 
{
    linha* linhas;
    int tamanho;
} sistema;

sistema criar_sistema(int size, float coeficientes[][size + 1]) 
{
    linha* linhas = (linha*) malloc(size * sizeof(linha));

    for(int j = 0; j < size; j++) 
    {
        float* coefs = (float*) malloc(size * sizeof(float));

        for(int i = 0; i < size; i++) 
        {
            *(coefs+i) = coeficientes[j][i];
        }

        linha lin = {
            .coeficientes = coefs,
            .tamanho = size,
            .termo_independente = coeficientes[j][size],
        };

        *(linhas+j) = lin;
    }

    sistema sis = {
        .linhas = linhas,
        .tamanho = size,
    };

    return sis;
}

void print_sistema(sistema* s)
{
    for(int j = 0; j < s->tamanho; j++) 
    {
        linha l = *(s->linhas+j);
        printf("%i: ", j);
        for(int i = 0; i < s->tamanho; i++) 
        {
            printf("%6.3f, ", *(l.coeficientes+i));
        }
        printf("%6.3f\n", l.termo_independente);
    }
}

void aplicar_gauss(sistema* s) 
{
    for(int k = 0; k < s->tamanho - 1; k++) 
    {
        linha linha_k = *(s->linhas+k);

        float a_kk = *(linha_k.coeficientes + k);
        
        for(int i = k + 1; i < s->tamanho; i++) 
        {
            linha linha_i = *(s->linhas+i);

            float a_ik = *(linha_i.coeficientes + k);
            float m_ik = a_ik / a_kk;

            for(int j = 0; j < s->tamanho; j++) 
            {
                float L_kj = *(linha_k.coeficientes + j);
                float* L_ij = (linha_i.coeficientes + j);
                *L_ij = *L_ij - (m_ik * L_kj);
            }

            float b_i = linha_i.termo_independente - (m_ik * linha_k.termo_independente);
            (s->linhas+i)->termo_independente = b_i;
        }
    }
}

float* resolver_escalonado(sistema* s) 
{
    float* resultados = (float*) malloc(s->tamanho * sizeof(float));
    
    for(int i = s->tamanho - 1; i >= 0; i--) 
    {
        linha L_i = *(s->linhas+i);
        float* x_i = resultados+i;

        float a_nn = *(L_i.coeficientes+i);
        float b_nn = L_i.termo_independente;

        float sum = 0;
        for (int j = i + 1; j < s->tamanho; j++) 
        {
            float x_j = *(resultados+j);
            float a_ij = *(L_i.coeficientes+j);
            sum += x_j * a_ij;
        }

        *(x_i) = (b_nn - sum) / a_nn;
    }

    return resultados;
}

int main() 
{
    /* 
     *  2x + 1y - 1z =  3
     *  4x - 2y + 5z = 16
     * -2x + 1y + 3z =  1
     */

    // int size = 3;
    // float coeficientes[3][4] = 
    // {
    //     {  2,  1, -1,  3 },
    //     {  4, -2,  5, 16 },
    //     { -2,  1,  3,  1 }
    // };

    /* 
     * 2w +  2x +  1y +  1z =  7
     * 1w + -1x +  1y + -1z =  1
     * 3w +  2x + -3y + -2z =  4
     * 4w +  3x +  2y +  1z = 12
     */

    int size = 4;
    float coeficientes[4][5] = 
    {
        {  2,  2,  1,  1,  7 },
        {  1, -1,  1, -1,  1 },
        {  3,  2, -3, -2,  4 },
        {  4,  3,  2,  1, 12 }
    };

    sistema s = criar_sistema(size, coeficientes);
    printf("SISTEMA ANTES:\n");
    print_sistema(&s);

    aplicar_gauss(&s);

    printf("SISTEMA DEPOIS:\n");
    print_sistema(&s);

    float* resultados = resolver_escalonado(&s);

    printf("RESULTADOS:\n");
    for(int i = 0; i < size; i++) 
    {
        printf("x_%i = %8.3f\n", i, *(resultados+i));
    }

    return 0;
}