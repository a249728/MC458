#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

void imprime_array_ui32(uint32_t n, const uint32_t *array)
{   
    printf("[");
    for (uint32_t i = 0; i < n; i++)
    {
        printf("%u, ", array[i]);
    }
    printf("]\n");
}

void imprime_array_i32(uint32_t n, const int32_t *array)
{
    printf("[");
    for (uint32_t i = 0; i < n; i++)
    {
        printf("%d, ", array[i]);
    }
    printf("]\n");
}

void le_linha_ui32(FILE *fptr, char *buffer, size_t buffer_size, uint32_t n, uint32_t *out)
{
    for (uint32_t i = 0; i < n; i++)
    {
        if (fgets(buffer, (int)buffer_size, fptr) == NULL)
        {
            fprintf(stderr, "Erro ao ler linha de uint32.\n");
            exit(EXIT_FAILURE);
        }
        out[i] = (uint32_t)strtoul(buffer, NULL, 10);
    }
}

void le_linha_i32(FILE *fptr, char *buffer, size_t buffer_size, uint32_t n, int32_t *out)
{
    for (uint32_t i = 0; i < n; i++)
    {
        if (fgets(buffer, (int)buffer_size, fptr) == NULL)
        {
            fprintf(stderr, "Erro ao ler linha de int32.\n");
            exit(EXIT_FAILURE);
        }
        out[i] = (int32_t)strtol(buffer, NULL, 10);
    }
}

void imprime_polinomio(uint32_t grau, const uint32_t *coef)
{
    // Considera que é um polinomio de grau 'grau', e que coef[0] é o maior coeficiente
    for (int32_t i = 0; i < (int32_t)grau; i++)
    {
        printf("%uX^%d + ", coef[i], grau - i);
    }
    printf("%u\n", coef[grau]);
}

