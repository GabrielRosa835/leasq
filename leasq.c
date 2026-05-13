#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

typedef struct {
    float* coeficientes;
    int tamanho;
    float termo_independente;
} linha_t;

typedef struct {
    linha_t* linhas;
    int tamanho;
} sistema_t;

typedef struct {
    float x;
    float f_x;
} ponto_t;

typedef struct {
    ponto_t* pontos;
    int tamanho;
} tabela_t;

tabela_t adaptar_tabela(int tamanho, float valores[][2]);
float** criar_matriz(int grau, const tabela_t* tabela);
tabela_t aquisitar_tabela();

float** aquisitar_sistema(int tamanho);
float** adaptar_valores(int tamanho, float matriz[][tamanho+1]);
sistema_t criar_sistema(int tamanho, float** matriz_determinante);
void imprimir_sistema(sistema_t* s);
void aplicar_gauss(sistema_t* s);
float* resolver_escalonado(sistema_t* s);

// ==================================================================



// ==================================================================

int main() 
{
    int escolha_metodo = -1;
    int escolha_opcao = -1;
    float** linhas = NULL;
    int tamanho = 0;

    printf("ESCOLHA UMA OPÇÃO:\n");
    printf("Eliminação de Gauss -> 0\n");
    printf("Mínimos Quadrados   -> 1\n");

    while(escolha_metodo < 0 || escolha_metodo > 1) 
    {
        printf("> ");
        scanf("%i", &escolha_metodo);
    }

    if (escolha_metodo == 0) 
    {
        printf("ESCOLHA UMA OPÇÃO:\n");
        printf("Sistema 3x3 de teste -> 0\n");
        printf("Sistema 4x4 de teste -> 1\n");
        printf("Sistema customizado  -> 2...6\n");
    
        while(escolha_opcao < 0 || escolha_opcao > 6) 
        {
            printf("> ");
            scanf("%i", &escolha_opcao);
        }

        if (escolha_opcao == 0) 
        {
            float coeficientes[3][3+1] = 
            {
                {  2,  1, -1,  3 },  //   2x + 1y - 1z =  3
                {  4, -2,  5, 16 },  //   4x - 2y + 5z = 16
                { -2,  1,  3,  1 }   //  -2x + 1y + 3z =  1
            };
            tamanho = 3;
            linhas = adaptar_valores(tamanho, coeficientes);
        } 
        else if (escolha_opcao == 1) 
        {
            float coeficientes[4][4+1] = 
            {
                {  2,  2,  1,  1,  7 },  //  2w + 2x + 1y + 1z =  7
                {  1, -1,  1, -1,  1 },  //  1w - 1x + 1y - 1z =  1
                {  3,  2, -3, -2,  4 },  //  3w + 2x - 3y - 2z =  4
                {  4,  3,  2,  1, 12 }   //  4w + 3x + 2y + 1z = 12
            };
            tamanho = 4;
            linhas = adaptar_valores(tamanho, coeficientes);
        }
        else
        {
            tamanho = escolha_opcao;
            linhas = aquisitar_sistema(tamanho);
        }
    } 
    else 
    {
        printf("ESCOLHA UMA OPÇÃO:\n");
        printf("Dados de teste 1   -> 0\n");
        printf("Dados de teste 2   -> 1\n");
        printf("Dados customizados -> 2\n");
    
        while(escolha_opcao < 0 || escolha_opcao > 2) 
        {
            printf("> ");
            scanf("%i", &escolha_opcao);
        }

        if (escolha_opcao == 0)
        {
            float valores[4][2] = 
            {
                { -1,  0 },
                {  0, -1 },
                {  1,  0 },
                {  2,  7 },
            };
            tabela_t tabela = adaptar_tabela(4, valores);
            int grau = 3;
            tamanho = grau;
            linhas = criar_matriz(grau, &tabela);
        }
        else if (escolha_opcao == 1)
        {
            float valores[4][2] = 
            {
                { -1,  0 },
                {  0, -1 },
                {  1,  0 },
                {  2,  7 },
            };
            tabela_t tabela = adaptar_tabela(4, valores);
            int grau = 5;
            tamanho = grau;
            linhas = criar_matriz(grau, &tabela);
        }
        else 
        {
            tabela_t tabela = aquisitar_tabela();
            int grau = 0;

            printf("Digite o grau de aproximação (1 < grau < 6):\n");

            while(grau < 1 || grau > 6) 
            {
                printf("> ");
                scanf("%i", &grau);
            }

            tamanho = grau + 1;
            linhas = criar_matriz(tamanho, &tabela);
        }
    }

    sistema_t s = criar_sistema(tamanho, linhas);

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

float* resolver_escalonado(sistema_t* s) 
{
    float* resultados = (float*) malloc(s->tamanho * sizeof(float));
    
    for(int i = s->tamanho - 1; i >= 0; i--) 
    {
        linha_t L_i = *(s->linhas+i);
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

void aplicar_gauss(sistema_t* s) 
{
    for(int k = 0; k < s->tamanho - 1; k++) 
    {
        linha_t linha_k = *(s->linhas+k);

        float a_kk = *(linha_k.coeficientes + k);
        
        for(int i = k + 1; i < s->tamanho; i++) 
        {
            linha_t linha_i = *(s->linhas+i);

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

void imprimir_sistema(sistema_t* s)
{
    for(int j = 0; j < s->tamanho; j++) 
    {
        linha_t l = *(s->linhas+j);
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

sistema_t criar_sistema(int size, float** matriz_determinante)
{
    linha_t* linhas = (linha_t*) malloc(size * sizeof(linha_t));

    for(int j = 0; j < size; j++) 
    {
        float* coefs = (float*) malloc(size * sizeof(float));
        float* coefs_matriz = *(matriz_determinante+j);

        for(int i = 0; i < size; i++) 
        {
            *(coefs+i) = *(coefs_matriz+i);
        }

        linha_t lin = {
            .coeficientes = coefs,
            .tamanho = size,
            .termo_independente = *(coefs_matriz+size),
        };

        *(linhas+j) = lin;
    }

    sistema_t sis = {
        .linhas = linhas,
        .tamanho = size,
    };

    return sis;
}

// ==================================================================

float** aquisitar_sistema(int tamanho)
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
            printf("> a_%i%i: ", j+1, i+1);
            scanf("%f", linha+i);
        }
        printf("> b_%i: ", j+1);
        scanf("%f", linha+tamanho);
    }

    return sistema;
}

// ==================================================================

float** adaptar_valores(int tamanho, float matriz[][tamanho+1]) 
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

// ==================================================================

tabela_t aquisitar_tabela()
{
    int tamanho;

    do 
    {
        printf("Digite a quantidade de pontos presentes (0 < quantidade < 7): \n");
        printf("> ");
        scanf("%i", &tamanho);
    }
    while (tamanho <= 0);

    tabela_t t = 
    {
        .pontos = (ponto_t*) malloc(tamanho * sizeof(ponto_t)),
        .tamanho = tamanho,
    };

    printf("Insira os pontos tabelados, digitando o valor de xi e f(xi) alternadamente:\n");

    for(int i = 0; i < t.tamanho; i++)
    {
        if (i != 0) printf("---+-------------\n");
        ponto_t p = *(t.pontos+i);
        printf(" %d | x    = ", i+1);
        scanf("%f", &(p.x));
        printf("   | f(x) = ");
        scanf("%f", &(p.f_x));
    }

    return t;
}

// ==================================================================

float** criar_matriz(int grau, const tabela_t* tabela)
{
    float** linhas = (float**) malloc(grau * sizeof(float*));

    // Sistema
    for(int i = 0; i < grau; i++) 
    {
        float* linha = (float*) malloc((grau + 1) * sizeof(float));
        *(linhas+i) = linha;

        // Equação
        for (int j = 0; j < grau; j++) 
        {
            float coeficiente = 0;

            // Somatória
            for (int k = 0; k < tabela->tamanho; k++) 
            {
                float x = ((tabela->pontos)+k)->x;
                coeficiente += pow(x, i + j);
            }
            *(linha+j) = coeficiente;
        }

        float termo_independente = 0;
        for(int k = 0; k < tabela->tamanho; k++)
        {
            float x = ((tabela->pontos)+k)->x;
            float f_x = ((tabela->pontos)+k)->f_x;
            termo_independente += f_x * pow(x, i);
        }
        *(linha+grau) = termo_independente;
    }

    return linhas;
}

// ==================================================================

tabela_t adaptar_tabela(int tamanho, float valores[][2])
{
    tabela_t t = 
    {
        .pontos = (ponto_t*) malloc(tamanho * sizeof(ponto_t)),
        .tamanho = tamanho,
    };

    for(int i = 0; i < tamanho; i++) 
    {
        ponto_t* p = t.pontos+i;
        p->x = valores[i][0];
        p->f_x = valores[i][1];
    }

    return t;
}