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
    fptr = fopen("input_teste.txt", "r");
    char buffer[100];

    // Pega o número de casos de teste
    int num_testes;
    fgets(buffer, sizeof(buffer), fptr);
    num_testes = atoi(buffer);


    for (int i = 0; i < num_testes; i++)
    {
        printf("Teste %d:\n\n", i + 1);
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

        if (grau < 50){
            printf("Polinômio G(X): \n");
            imprime_polinomio(grau, G);
        }

        // Lê o número de mensagens a serem codificadas
        fgets(buffer, sizeof(buffer), fptr);
        k = atoi(buffer);


        for (int j = 0; j < k; j++)
        {
            uint32_t *mensagem = (uint32_t *)malloc((grau + 1) * sizeof(uint32_t));
            le_linha_ui32(fptr, buffer, sizeof(buffer), grau + 1, mensagem);
            printf("\n");
            if (grau < 50){
                printf("Mensagem M(X): \n");
                imprime_polinomio(grau, mensagem);
            }

            uint32_t *P = multiplica_polinomio(grau, G, mensagem);

            if (grau < 50){
                printf("Polinômio codificado P(X): \n");
                imprime_polinomio(2*grau, P);
            }

            if (verifica_erros(grau * 2, P, grau, raizes)){
                printf("Erros detectados!\n");
            } else {
                printf("Nenhum erro detectado.\n");
            }

            uint32_t *M = decodifica_mensagem( 2 * grau, P, grau, raizes);

            if (grau < 50){
                printf("Mensagem decodificada M'(X): \n");
                imprime_polinomio(grau, M);
            }


            adiciona_ruido(grau * 2, P, 2); // Adiciona um erro

            if (grau < 50){
                printf("Polinômio codificado P(X) com ruído em 2 dos coeficientes: \n");
                imprime_polinomio(2*grau, P);
            }

            if (verifica_erros(grau * 2, P, grau, raizes)){
                printf("Erros detectados após adicionar ruído!\n");
            } else {
                printf("Nenhum erro detectado após adicionar ruído.\n");
            }

            free(mensagem);
            free(P);
        }
        free(raizes);
        free(G);
        printf("\n");
    }
    fclose(fptr);
}
