#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include "io.h"
#include "polinomios.h"


int main(int argc, char *argv[])
{   
    clock_t start, end;
    double cpu_time_used;

    FILE *fptr;
    fptr = fopen("input.txt", "r");
    char buffer[100];

    // Pega o número de casos de teste
    int num_testes;
    fgets(buffer, sizeof(buffer), fptr);
    num_testes = atoi(buffer);


    for (int i = 0; i < num_testes; i++)
    {
        int grau;
        int k;
        // Lê o valor de grau dos polinomios
        fgets(buffer, sizeof(buffer), fptr);
        grau = atoi(buffer);

        // Lê as raízes de G(X), lê com sinal
        int32_t *raizes = (int32_t *)malloc((grau) * sizeof(int32_t));
        le_linha_i32(fptr, buffer, sizeof(buffer), grau, raizes);

        // Lê os coeficientes de G(x)
        uint32_t *G = (uint32_t *)malloc((grau + 1) * sizeof(uint32_t));
        le_linha_ui32(fptr, buffer, sizeof(buffer), grau + 1, G);

        // Lê o número de mensagens a serem codificadas
        fgets(buffer, sizeof(buffer), fptr);
        k = atoi(buffer);


        for (int j = 0; j < k; j++)
        {
            // Lê a mensagem M(X)
            uint32_t *mensagem = (uint32_t *)malloc((grau + 1) * sizeof(uint32_t));
            le_linha_ui32(fptr, buffer, sizeof(buffer), grau + 1, mensagem);

            // Mede o tempo de multiplicar os polinômios
            start = clock();
            uint32_t *P = multiplica_polinomio(grau, G, mensagem);
            end = clock();
            cpu_time_used = cpu_time_used + ((double)(end - start));
            
            // Verifica erros antes de adicionar ruído
            if(verifica_erros(grau * 2, P, grau, raizes)){
                printf("Erro detectado antes de adicionar ruído\n");
            };

            // Exibe o polinômio resultante da multiplicação
            imprime_polinomio(grau * 2, P);

            adiciona_ruido(grau * 2, P, 5); // Adiciona ruído em 5 coeficientes aleatórios
            if(verifica_erros(grau * 2, P, grau, raizes)){
                printf("Ruído identificado\n");
            }; 
            free(mensagem);
            free(P);
        }
        cpu_time_used = cpu_time_used / CLOCKS_PER_SEC;
        printf("Tempo gasto para grau %d: %f segundos\n", grau, cpu_time_used);
        free(raizes);
        free(G);
    }
    fclose(fptr);
}
