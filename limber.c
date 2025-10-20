#include "limber.h"
#include <stdlib.h>
#include <stdbool.h>
#include <ctype.h>
#include <string.h>
#include <stdio.h>

#define LIMB_MAX ((limb_t)(0 - 1))
#define BITS_PER_LIMB ((sizeof(limb_t)) * 8)
#define HALF_BITS (BITS_PER_LIMB / 2)
#define LOWER_MASK (LIMB_MAX >> HALF_BITS)

const static struct limber limber_null = { .size = 0, .alloc = 0, .sign = 0, .limbs = NULL };
#define LIMBER_NULL limber_null

void limber_init(Limber l)
{
	l->size = 0;
	l->alloc = 0;
	l->sign = 0;
	l->limbs = NULL;
}

void limber_clear(Limber l)
{
	free(l->limbs);
	*l = LIMBER_NULL;
}

static void limber_swap(Limber a, Limber b)
{
	struct limber temp;

	temp = *a;
	*a = *b;
	*b = temp;
}

bool limber_is_null(Limber l)
{
	if (l == NULL) {
		return true;
	}

	return l->size == 0 && l->alloc == 0 && l->sign == 0 && l->limbs == NULL;
}

#define LIMBER_GROW_FACTOR 2
static void limber_realloc(Limber l, int min_size)
{
	int new_size = min_size * LIMBER_GROW_FACTOR;
	limb_t *temp = realloc(l->limbs, sizeof(limb_t) * new_size);
	if (temp) {
		l->alloc = new_size;
		l->limbs = temp;
	} else {
		// TODO(): Handle realloc fail
	}
}

static bool limber_is_full(Limber l)
{
	return (l->size == l->alloc);
}

void limber_set(Limber rop, Limber op)
{
	limber_clear(rop);
	if (op == NULL || limber_is_null(op)) {
		return;
	}

	rop->limbs = malloc(sizeof(limb_t) * op->alloc);
	if (rop->limbs) {
		rop->size = op->size;
		rop->alloc = op->alloc;
		rop->sign = op->sign;
		memcpy(rop->limbs, op->limbs, sizeof(limb_t) * rop->size);
	} else {
		*rop = LIMBER_NULL;
	}
}

static void limber_normalize(Limber l)
{
	while (l->size > 1 && l->limbs[l->size - 1] == 0) {
		l->size--;
	}

	if (l->size == 1 && l->limbs[0] == 0) {
		l->sign = 0;
	}
}

void limber_add_limb(Limber rop, Limber op1, limb_t op2)
{
	Limber temp_rop;
	limber_init(temp_rop);

	limber_set(temp_rop, op1);

	limb_t carry = op2;
	for (int i = 0; i < temp_rop->size; i++) {
		limb_t sum = temp_rop->limbs[i] + carry;
		carry = sum < temp_rop->limbs[i];
		temp_rop->limbs[i] = sum;
	}

	if (carry > 0) {
		if (limber_is_full(temp_rop)) {
			limber_realloc(temp_rop, temp_rop->size + 1);
			// TODO(): Handle realloc fail
		}
		
		temp_rop->limbs[temp_rop->size++] = carry;
	}

	limber_normalize(temp_rop);
	limber_swap(temp_rop, rop);
	limber_clear(temp_rop);
}

void limber_sub_limb(Limber rop, Limber op1, limb_t op2)
{
	Limber temp_rop;
	limber_set(temp_rop, op1);
	if (limber_is_null(temp_rop)) {
		return;
	}

	limb_t borrow = op2;
	for (int i = 0; i < temp_rop->size; i++) {
		limb_t sub = temp_rop->limbs[i] - borrow;
		borrow = sub > temp_rop->limbs[i];
		temp_rop->limbs[i] = sub;
	}

	if (borrow > 0) {
		// TODO(): Handle the case when op1 < op2
		// This will only happen if op1 has just one limb
		// and its value is less than op2
	}

	limber_normalize(temp_rop);
	limber_swap(temp_rop, rop);
	limber_clear(temp_rop);
}

static void limber_mul_single(limb_t *mul_hi, limb_t *mul_lo, limb_t op1, limb_t op2)
{
	limb_t op1_lo = op1 & LOWER_MASK;
	limb_t op1_hi = op1 >> HALF_BITS;
	limb_t op2_lo = op2 & LOWER_MASK;
	limb_t op2_hi = op2 >> HALF_BITS;

	limb_t p0 = op1_lo * op2_lo;
	limb_t p1 = op1_lo * op2_hi;
	limb_t p2 = op1_hi * op2_lo;
	limb_t p3 = op1_hi * op2_hi;

	limb_t c0 = (p0 & LOWER_MASK);
	limb_t c1 = (p0 >> HALF_BITS) + (p1 & LOWER_MASK) + (p2 & LOWER_MASK);
	limb_t c2 = (p1 >> HALF_BITS) + (p2 >> HALF_BITS) + (p3 & LOWER_MASK) + (c1 >> HALF_BITS);
	limb_t c3 = (p3 >> HALF_BITS) + (c2 >> HALF_BITS);

	*mul_lo = ((c1 & LOWER_MASK) << HALF_BITS) | (c0 & LOWER_MASK);
	*mul_hi = ((c3 & LOWER_MASK) << HALF_BITS) | (c2 & LOWER_MASK);
}

void limber_mul_limb(Limber rop, Limber op1, limb_t op2)
{
	Limber temp_rop;
	limber_init(temp_rop);

	limber_set(temp_rop, op1);
	if (limber_is_null(temp_rop)) {
		return;
	}

	limb_t carry = 0;
	for (int i = 0; i < temp_rop->size; i++) {
		limb_t current_limb = temp_rop->limbs[i];
		limb_t mul_hi = 0;
		limb_t mul_lo = 0;

		limber_mul_single(&mul_hi, &mul_lo, current_limb, op2);

		temp_rop->limbs[i] = mul_lo + carry;
		carry = mul_hi + (temp_rop->limbs[i] < mul_lo);
	}

	if (carry > 0) {
		if (limber_is_full(temp_rop)) {
			limber_realloc(temp_rop, temp_rop->size + 1);
			// TODO(): Handle realloc fail
		}
		
		temp_rop->limbs[temp_rop->size++] = carry;
	}

	limber_normalize(temp_rop);
	limber_swap(temp_rop, rop);
	limber_clear(temp_rop);
}


static bool is_greater_than(limb_t op1_hi, limb_t op1_lo, limb_t op2_hi, limb_t op2_lo)
{
	if (op1_hi > op2_hi) {
		return true;
	}

	if (op1_hi < op2_hi) {
		return false;
	}

	return op1_lo > op2_lo;
}

static int clz(limb_t n)
{
	if (n == 0) {
		return BITS_PER_LIMB;
	}

	int count = 0;
	limb_t mask = (limb_t) 1 << (BITS_PER_LIMB - 1);
	while ((n & mask) == 0) {
		count++;
		mask >>= 1;
	}

	return count;
}

static void divrem_2by1(limb_t *q_ptr, limb_t *r_ptr, limb_t p_hi, limb_t p_lo, limb_t q)
{
	if (p_hi == 0) {
		if (q_ptr) {
			*q_ptr = p_lo / q;
		}

		if (r_ptr) {
			*r_ptr = p_lo % q;
		}
		return;
	}

	int s = clz(q);
	limb_t qn = q << s;
	limb_t qn_hi = qn >> HALF_BITS;

	limb_t pn_hi = (p_hi << s) | (p_lo >> (BITS_PER_LIMB - s));
	limb_t pn_lo = p_lo << s;

	limb_t est = 0;
	if (p_hi >= q) {
		est = LIMB_MAX;
	} else {
		if (pn_hi >= qn_hi) {
			est = LIMB_MAX;
		} else {
			est = pn_hi / qn_hi;
		}
	}

	limb_t prod_hi = 0;
	limb_t prod_lo = 0;
	limber_mul_single(&prod_hi, &prod_lo, qn, est);

	while (is_greater_than(prod_hi, prod_lo, pn_hi, pn_lo)) {
		est--;
		limber_mul_single(&prod_hi, &prod_lo, qn, est);
	}

	limb_t rem_hi = pn_hi - prod_hi - (pn_lo < prod_lo);
	limb_t rem_lo = pn_lo - prod_lo;

	if (r_ptr) {
		*r_ptr = (rem_hi << (BITS_PER_LIMB - s)) | (rem_lo >> s);
	}

	if (q_ptr) {
		*q_ptr = est;
	}
}


void limber_div_limb(Limber rop, limb_t *remainder, Limber dividend, limb_t divider)
{
	// TODO: divider can't be zero
	Limber temp_rop;
	limber_set(temp_rop, dividend);
	if (limber_is_null(temp_rop)) {
		if (remainder) {
			*remainder = 0;
		}

		return;
	}

	limb_t partial_remainder = 0;
	for (int i = temp_rop->size - 1; i >= 0; i--) {
		limb_t dividend_hi = partial_remainder;
		limb_t dividend_lo = temp_rop->limbs[i];
		
		divrem_2by1(&temp_rop->limbs[i], &partial_remainder, dividend_hi, dividend_lo, divider);
	}

	if (remainder) {
		*remainder = partial_remainder;
	}

	limber_normalize(temp_rop);
	limber_swap(temp_rop, rop);
	limber_clear(temp_rop);
}

bool limber_is_zero(Limber op)
{
	return op->size == 1 && op->limbs[0] == 0 && op->sign == 0;
}

void limber_set_from_str(Limber rop, const char *str, int base)
{
	while (isspace(*str)) {
		str++;
	}

	int sign = 1;
	if (*str == '-') {
		sign = -1;
		str++;
	} else if (*str == '+') {
		str++;
	}

	while (*str) {
		int digit = *str - '0';

		limber_mul_limb(rop, rop, base);
		limber_add_limb(rop, rop, digit);

		str++;
	}

	if (limber_is_zero(rop)) {
		rop->sign = 0;
	} else {
		rop->sign = sign;
	}
}

void limber_debug_print(const char *name, Limber l)
{
	printf("DEBUG: %s\n", name);

	if (limber_is_null(l)) {
		printf("	-> NULL LIMBER\n");
		return;
	}

	printf("	-> SIGN:  %d\n", l->sign);
	printf("	-> SIZE:  %d\n", l->size);
	printf("	-> ALLOC: %d\n", l->alloc);

	if (l->size == 0) {
		printf("	-> LIMBS: [EMPTY]\n");
		return;
	}

	printf("	-> LIMBS: [");

	for (int i = l->size - 1; i >= 0; i--) {
		printf(" %016lX", l->limbs[i]);
	}
	printf("]\n");
}

