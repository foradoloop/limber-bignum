#include <stdio.h>
#include <stdlib.h>
#include "limber.h" // Certifique-se de que o nome do seu header está correto

/**
 * @brief Função de depuração para imprimir o estado interno de um Limber.
 * Essencial para testar sem uma função to_string.
 */
/*static void limber_debug_print(const char *name, Limber l) {
    printf("--- DEBUG: %s ---\n", name);
    if (l == NULL || l->limbs == NULL) {
        printf("  -> Limber nulo ou nao inicializado.\n");
        return;
    }
    printf("  -> sign:  %d\n", l->sign);
    printf("  -> size:  %d\n", l->size);
    printf("  -> limbs[0] (hex): 0x%lX\n", l->limbs[0]);
    printf("  -> limbs[0] (dec): %lu\n", l->limbs[0]);
    printf("------------------------\n\n");
}
*/

int main(void) {
    printf("Iniciando testes simples com Limbers...\n\n");

    // --- 1. Inicialização ---
    // Sempre inicialize todas as variáveis.
    Limber zero, a, resultado;
    limber_init(zero);
    limber_init(a);
    limber_init(resultado);
    
    // O 'zero' será nosso ponto de partida para criar outros números.
    limber_debug_print("Variavel 'zero' inicializada", zero);

    // --- 2. Criando um número (usando seu truque) ---
    // Vamos criar o número 100 em 'a'.
    // a = zero + 100
    printf("=== Etapa 1: Criando o numero 100 ===\n");
    limber_add_limb(a, zero, 100);
    limber_debug_print("Variavel 'a'", a);
    printf(">>> VERIFIQUE: 'a' deve ter o valor 100 (ou 0x64 em hex).\n\n");

    // --- 3. Testando Adição ---
    // Vamos calcular: resultado = a + 50
    printf("=== Etapa 2: Testando adicao (100 + 50) ===\n");
    limber_add_limb(resultado, a, 50);
    limber_debug_print("Variavel 'resultado' apos adicao", resultado);
    printf(">>> VERIFIQUE: 'resultado' deve ter o valor 150 (ou 0x96 em hex).\n\n");

    // --- 4. Testando Multiplicação (com aliasing) ---
    // Vamos calcular: resultado = resultado * 10
    // Este é um teste crucial para sua refatoração anti-aliasing!
    printf("=== Etapa 3: Testando multiplicacao (150 * 10) ===\n");
    limber_mul_limb(resultado, resultado, 10);
    limber_debug_print("Variavel 'resultado' apos multiplicacao", resultado);
    printf(">>> VERIFIQUE: 'resultado' deve ter o valor 1500 (ou 0x5DC em hex).\n\n");

    // --- 5. Testando Divisão ---
    // Vamos calcular: resultado = 1507 / 100
    printf("=== Etapa 4: Testando divisao (1507 / 100) ===\n");
    // Primeiro, criamos o número 1507 em 'a'.
    limber_add_limb(a, zero, 1507);
    limber_debug_print("Novo valor de 'a' para a divisao", a);
    
    limb_t resto = 0;
    limber_div_limb(resultado, &resto, a, 100);

    printf("Resultados da divisao:\n");
    limber_debug_print("Quociente", resultado);
    printf("Resto: %lu\n", resto);
    printf(">>> VERIFIQUE: Quociente deve ser 15 (0xF) e Resto deve ser 7.\n\n");

    // --- 6. Limpeza Final ---
    printf("Liberando memoria...\n");
    limber_clear(zero);
    limber_clear(a);
    limber_clear(resultado);

    printf("Testes concluidos!\n");

    return 0;
}
