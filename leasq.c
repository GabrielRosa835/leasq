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

float** construir_matriz(int tamanho);
float** adaptar_matriz(int tamanho, float matriz[][tamanho+1]);
sistema criar_sistema(int size, float** matriz_determinante);
void imprimir_sistema(sistema* s);
void aplicar_gauss(sistema* s);
float* resolver_escalonado(sistema* s);

// ==================================================================

int main() 
{
    int escolha = -1;

    printf("ESCOLHA UMA OPÇÃO:\n");
    printf("Sistema 3x3 de teste -> 0\n");
    printf("Sistema 4x4 de teste -> 1\n");
    printf("Sistema customizado  -> 2...6\n");

    while(escolha < 0 || escolha > 6) 
    {
        printf(" - ");
        scanf("%i", &escolha);
    }

    float** linhas = NULL;
    int tamanho = 0;

    if (escolha == 0) 
    {
        float coeficientes[3][3+1] = 
        {
            {  2,  1, -1,  3 },  //   2x + 1y - 1z =  3
            {  4, -2,  5, 16 },  //   4x - 2y + 5z = 16
            { -2,  1,  3,  1 }   //  -2x + 1y + 3z =  1
        };
        tamanho = 3;
        linhas = adaptar_matriz(tamanho, coeficientes);
    } 
    else if (escolha == 1) 
    {
        float coeficientes[4][4+1] = 
        {
            {  2,  2,  1,  1,  7 },  //  2w + 2x + 1y + 1z =  7
            {  1, -1,  1, -1,  1 },  //  1w - 1x + 1y - 1z =  1
            {  3,  2, -3, -2,  4 },  //  3w + 2x - 3y - 2z =  4
            {  4,  3,  2,  1, 12 }   //  4w + 3x + 2y + 1z = 12
        };
        tamanho = 4;
        linhas = adaptar_matriz(tamanho, coeficientes);
    }
    else
    {
        tamanho = escolha;
        linhas = construir_matriz(tamanho);
    }

    sistema s = criar_sistema(tamanho, linhas);

    printf("SISTEMA ANTES:\n");
    imprimir_sistema(&s);

    aplicar_gauss(&s);

    printf("SISTEMA DEPOIS:\n");
    imprimir_sistema(&s);

    float* resultados = resolver_escalonado(&s);

    printf("RESULTADOS:\n");
    for(int i = 0; i < tamanho; i++) 
    {
        printf("x_%i = %8.3f\n", i+1, *(resultados+i));
    }

    return 0;
}

// ==================================================================

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

// ==================================================================

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

// ==================================================================

void imprimir_sistema(sistema* s)
{
    for(int j = 0; j < s->tamanho; j++) 
    {
        linha l = *(s->linhas+j);
        printf("%i: ", j+1);

        for(int i = 0; i < s->tamanho; i++) 
        {
            if (i > 0) printf(" + ");
            printf("%7.3f*x_%i", *(l.coeficientes+i), i+1);
        }
        
        printf(" = %7.3f", l.termo_independente);
        printf("\n");
    }
}

// ==================================================================

sistema criar_sistema(int size, float** matriz_determinante)
{
    linha* linhas = (linha*) malloc(size * sizeof(linha));

    for(int j = 0; j < size; j++) 
    {
        float* coefs = (float*) malloc(size * sizeof(float));
        float* coefs_matriz = *(matriz_determinante+j);

        for(int i = 0; i < size; i++) 
        {
            *(coefs+i) = *(coefs_matriz+i);
        }

        linha lin = {
            .coeficientes = coefs,
            .tamanho = size,
            .termo_independente = *(coefs_matriz+size),
        };

        *(linhas+j) = lin;
    }

    sistema sis = {
        .linhas = linhas,
        .tamanho = size,
    };

    return sis;
}

// ==================================================================

float** construir_matriz(int tamanho)
{
    printf("TEMPLATE:\n");
    for(int j = 0; j < tamanho; j++) 
    {
        printf("%i: ", j);
        for(int i = 0; i < tamanho; i++) 
        {
            if (i > 0) printf(" + ");
            printf("a_%i%i*x_%i", j+1, i+1, i+1);
        }
        printf(" = b_%i", j+1);
        printf("\n");
    }

    float** sistema = (float**) malloc(tamanho * sizeof(float*));
    printf("VALORES:\n");

    for(int j = 0; j < tamanho; j++) 
    {
        float* linha = (float*) malloc((tamanho+1) * sizeof(float));
        *(sistema+j) = linha;

        for(int i = 0; i < tamanho; i++) 
        {
            printf("a_%i%i: ", j+1, i+1);
            scanf("%f", linha+i);
        }
        printf("b_%i: ", j+1);
        scanf("%f", linha+tamanho);
    }

    return sistema;
}

// ==================================================================

float** adaptar_matriz(int tamanho, float matriz[][tamanho+1]) 
{
    float** linhas = (float**) malloc(tamanho * sizeof(float*));
    for(int i = 0; i < tamanho; i++) 
    {
        float* linha = (float*) malloc((tamanho + 1) * sizeof(float));
        *(linhas+i) = linha;
        for (int j = 0; j < tamanho; j++) 
        {
            *(linha+j) = matriz[i][j];
        }
        *(linha+tamanho) = matriz[i][tamanho];
    }
    return linhas;
}