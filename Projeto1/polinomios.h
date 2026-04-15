#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>

// Função auxiliar para o caso base do Karatsuba (n <= 16)
uint32_t *multiplica_ingenua(uint32_t grau, const uint32_t *coef1, const uint32_t *coef2)
{
    uint32_t tamanho_resultado = 2 * grau + 1;
    uint32_t *resultado = (uint32_t *)calloc(tamanho_resultado, sizeof(uint32_t));
    
    for (uint32_t i = 0; i <= grau; i++) {
        for (uint32_t j = 0; j <= grau; j++) {
            resultado[i + j] += coef1[i] * coef2[j];
        }
    }
    
    return resultado;
}

uint32_t *multiplica_polinomio(uint32_t grau, const uint32_t *coef1, const uint32_t *coef2)
{
    // Caso base otimizado para n <= 16
    if (grau <= 16) {
        return multiplica_ingenua(grau, coef1, coef2);
    }

    uint32_t tamanho = grau + 1;
    uint32_t m = tamanho / 2; 
    
    uint32_t grau_alto = m - 1;
    uint32_t grau_baixo = (tamanho - m) - 1;

    // Alocação das metades
    uint32_t *A1 = (uint32_t *)malloc((grau_alto + 1) * sizeof(uint32_t));
    uint32_t *B1 = (uint32_t *)malloc((grau_alto + 1) * sizeof(uint32_t));
    uint32_t *A0 = (uint32_t *)malloc((grau_baixo + 1) * sizeof(uint32_t));
    uint32_t *B0 = (uint32_t *)malloc((grau_baixo + 1) * sizeof(uint32_t));

    for (uint32_t i = 0; i <= grau_alto; i++) {
        A1[i] = coef1[i];
        B1[i] = coef2[i];
    }
    for (uint32_t i = 0; i <= grau_baixo; i++) {
        A0[i] = coef1[m + i];
        B0[i] = coef2[m + i];
    }

    // Multiplica metades
    uint32_t *Z2 = multiplica_polinomio(grau_alto, A1, B1);
    uint32_t *Z0 = multiplica_polinomio(grau_baixo, A0, B0);

    // Soma das metades
    uint32_t grau_soma = grau_baixo;
    uint32_t *SomaA = (uint32_t *)calloc(grau_soma + 1, sizeof(uint32_t));
    uint32_t *SomaB = (uint32_t *)calloc(grau_soma + 1, sizeof(uint32_t));

    for (uint32_t i = 0; i <= grau_alto; i++) {
        SomaA[i + (grau_soma - grau_alto)] += A1[i];
        SomaB[i + (grau_soma - grau_alto)] += B1[i];
    }
    for (uint32_t i = 0; i <= grau_baixo; i++) {
        SomaA[i] += A0[i];
        SomaB[i] += B0[i];
    }

    // Multiplica as somas
    uint32_t *Z1_temp = multiplica_polinomio(grau_soma, SomaA, SomaB);

    // Montagem do resultado
    uint32_t *resultado = (uint32_t *)calloc(2 * grau + 1, sizeof(uint32_t));

    for (uint32_t i = 0; i <= 2 * grau_alto; i++) {
        resultado[i] += Z2[i];
    }

    for (uint32_t i = 0; i <= 2 * grau_baixo; i++) {
        resultado[2 * grau - 2 * grau_baixo + i] += Z0[i];
    }

    // Liberação de memória
    free(A1); free(B1); free(A0); free(B0);
    free(SomaA); free(SomaB);
    free(Z2); free(Z0); free(Z1_temp);

    return resultado;
}

// Verifica erros
int32_t avalia_polinomio(int32_t x, uint32_t grau, const uint32_t *coef)
{
    int32_t resultado = coef[0];
    for (uint32_t i = 1; i <= grau; i++)
    {
        resultado = (resultado * x) + (int32_t)coef[i];
    }
    return resultado;
}

// Decodifica
uint32_t *divide_por_binomio(uint32_t grau, uint32_t *coef, int32_t raiz)
{
    // Divide o polinomio de grau 'grau' por um binomio do tipo (x - raiz), onde raiz é uma raiz do polinomio, ou seja, P(raiz) == 0
    // Dividir por um binomio vai diminuir um grau
    uint32_t *resultado = (uint32_t *)calloc(grau, sizeof(uint32_t));
    resultado[0] = coef[0]; // O coeficiente de maior grau do resultado é igual ao coeficiente de maior grau do polinômio original
    for(uint32_t i = 1; i < grau; i++)
    {
        resultado[i] = coef[i] + raiz * resultado[i - 1];
    }
    return resultado;
}

uint32_t *decodifica_mensagem(uint32_t grau, uint32_t *coef, uint32_t n_raizes, int32_t *raizes)
{
    // Decodifica a mensagem original M(X) a partir do polinômio codificado P(X) e das raízes fornecidas
    uint32_t *coeficientes = (uint32_t *)calloc(grau + 1, sizeof(uint32_t));
    // Copia P
    for (uint32_t i = 0; i <= grau; i++)
    {
        coeficientes[i] = coef[i];
    }
    for (uint32_t i = 0; i < n_raizes; i++)
    {
        uint32_t *temp = coeficientes;
        coeficientes = divide_por_binomio(grau - i, temp, raizes[i]);
        free(temp); // Libera a memória da iteração anterior para evitar memory leak
    }
    return coeficientes;
}

_Bool verifica_erros(uint32_t grau, const uint32_t *coef, uint32_t n_raizes, const int32_t *raizes)
{
    for (uint32_t i = 0; i < n_raizes; i++)
    {
        if (avalia_polinomio(raizes[i], grau, coef) != 0)
        {
            return 1; // erro detectado
        }
    }
    return 0;
}

// Adiciona ruido
void adiciona_ruido(uint32_t grau, uint32_t *coef, uint32_t m)
{
    // Adiciona ruído em m coeficientes aleatórios do polinômio
    for (uint32_t i = 0; i < m; i++)
    {
        uint32_t indice = rand() % (grau + 1);
        uint32_t ruido = (rand() % 10) + 1; // Ruído entre 1 e 10
        coef[indice] += ruido;
    }
}