#ifndef LIMBER_H
#define LIMBER_H

typedef unsigned long int limb_t;

struct limber {
	int size;
	int alloc;
	int sign;
	limb_t *limbs;
};
typedef struct limber Limber[1];

void limber_init(Limber l);
void limber_clear(Limber l);

char *limber_to_str(Limber l, int base);
void limber_set_from_str(Limber rop, const char *str, int base);
void limber_debug_print(const char *name, Limber l);

void limber_add_limb(Limber rop, Limber op1, limb_t op2);
void limber_sub_limb(Limber rop, Limber op1, limb_t op2);
void limber_mul_limb(Limber rop, Limber op1, limb_t op2);
void limber_div_limb(Limber q, limb_t *r, Limber n, limb_t d);

void limber_add(Limber rop, Limber op1, Limber op2);
void limber_sub(Limber rop, Limber op1, Limber op2);
void limber_mul(Limber rop, Limber op1, Limber op2);

#endif
