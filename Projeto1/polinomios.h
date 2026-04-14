#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>



uint32_t *multiplica_polinomio(uint32_t grau, const uint32_t *coef1, const uint32_t *coef2)
{
    // Você deve implementar essa função
    // Os dois polinômios de entrada têm grau "grau", ou seja, cada um tem "grau + 1" coeficientes.
    // O resultado da multiplicação terá grau no máximo 2*grau, ou seja, terá no máximo 2*grau + 1 coeficientes.
    // Você deve implementar a multiplicação usando o método de Karatsuba, que é mais eficiente do que o método ingênuo para polinômios de grau grande.
    // Você pode criar sub-funções auxiliares
}

// Verifica erros
int32_t avalia_polinomio(int32_t x, uint32_t grau, const uint32_t *coef)
{
    // Você deve implementar essa função
    // Avalia o polinômio de grau "grau" com coeficientes "coef" no ponto "x"
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