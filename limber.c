#include "limber.h"
#include <stdlib.h>
#include <stdbool.h>
#include <ctype.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#define LIMB_MAX ((limb_t)(0 - 1))
#define BITS_PER_LIMB ((sizeof(limb_t)) * 8)
#define HALF_BITS (BITS_PER_LIMB / 2)
#define LOWER_MASK (LIMB_MAX >> HALF_BITS)

static const struct limber limber_null = { .size = 0, .alloc = 0, .sign = 0, .limbs = NULL };
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
	if (limber_is_null(rop)) {
		rop->limbs = malloc(sizeof(limb_t) * op->alloc);
	} else {
		if (rop->alloc < op->alloc) {
			limb_t *temp = realloc(rop->limbs, op->alloc);
			if (!temp) {
				// realoc fail
				limber_clear(rop);
				*rop = LIMBER_NULL;
			}
			rop->limbs = temp;
		}
	}

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

static int limber_clz(limb_t n)
{
	if (n == 0) {
		return BITS_PER_LIMB;
	}

	int count = 0;
	limb_t mask = (limb_t) 1 << (BITS_PER_LIMB - 1);
	while (mask != 0 && (n & mask) == 0) {
		count++;
		mask >>= 1;
	}

	return count;
}

static void divmnu(limb_t high, limb_t low, limb_t divisor, limb_t *q_hi_ptr, limb_t *q_lo_ptr, limb_t *r_ptr)
{
	limb_t quotient_high = 0;
	limb_t quotient_low  = 0;
	limb_t rem = 0;

	for (int i = BITS_PER_LIMB - 1; i >= 0; i--) {
		bool overflow = (rem >> (BITS_PER_LIMB - 1)) & 1;

		unsigned char current_bit = (high >> i) & 1;
		rem = (rem << 1) | current_bit;

		if (overflow || rem >= divisor) {
			rem -= divisor;
			quotient_high |= ((limb_t)1 << i);
		}
	}

	for (int i = BITS_PER_LIMB - 1; i >= 0; i--) {
		bool overflow = (rem >> (BITS_PER_LIMB - 1)) & 1;

		unsigned char current_bit = (low >> i) & 1;
		rem = (rem << 1) | current_bit;

		if (overflow || rem >= divisor) {
			rem -= divisor;
			quotient_low |= ((limb_t)1 << i);
		}
	}

	if (q_hi_ptr) {
		*q_hi_ptr = quotient_high;
	}

	if (q_lo_ptr) {
		*q_lo_ptr = quotient_low;
	}

	if (r_ptr) {
		*r_ptr = rem;
	}
}
/*
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

	int s = limber_clz(q);
	limb_t qn = q << s;
	limb_t qn_hi = qn >> HALF_BITS;

	limb_t pn_hi = (p_hi << s) | (p_lo >> (BITS_PER_LIMB - s));
	limb_t pn_lo = p_lo << s;

	limb_t est = 0;
	limb_t pn_hi_hi = pn_hi >> HALF_BITS;
	limb_t pn_hi_lo = pn_hi & LOWER_MASK;

	if ((pn_hi >> HALF_BITS) >= qn_hi) {
		est = LIMB_MAX;
	} else {
		est = ( (pn_hi_hi << HALF_BITS) | pn_hi_lo ) / qn_hi;
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

	while (rem_hi > 0 || rem_lo >= qn) {
		est++;

		limb_t old_rem_lo = rem_lo;
		rem_lo -= qn;
		if (rem_lo > old_rem_lo) {
			rem_hi--;
		}
	}

	if (r_ptr) {
		*r_ptr = (rem_hi << (BITS_PER_LIMB - s)) | (rem_lo >> s);
	}

	if (q_ptr) {
		*q_ptr = est;
	}
}
*/
bool limber_is_zero(Limber op)
{
        return op->size == 1 && op->limbs[0] == 0 && op->sign == 0;
}


void limber_div_limb(Limber rop, limb_t *remainder, Limber dividend, limb_t divider)
{
	// TODO: divider can't be zero
	Limber temp_rop;
	limber_init(temp_rop);
	limber_set(temp_rop, dividend);
	if (limber_is_null(temp_rop) || limber_is_zero(temp_rop)) {
		if (remainder) {
			*remainder = 0;
		}

		limber_clear(temp_rop);
		limber_realloc(temp_rop, 1);
		temp_rop->limbs[0] = 0;
		temp_rop->size = 1;
		temp_rop->sign = 0;
	} else {
		limb_t partial_remainder = 0;
		for (int i = temp_rop->size - 1; i >= 0; i--) {
			limb_t dividend_hi = partial_remainder;
			limb_t dividend_lo = temp_rop->limbs[i];
			divmnu(dividend_hi, dividend_lo, divider, NULL, &temp_rop->limbs[i], &partial_remainder);
			//divrem_2by1(&temp_rop->limbs[i], &partial_remainder, dividend_hi, dividend_lo, divider);
		}

		if (remainder) {
			*remainder = partial_remainder;
		}
	}

	limber_normalize(temp_rop);
	limber_swap(temp_rop, rop);
	limber_clear(temp_rop);
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

int limber_cmp(Limber op1, Limber op2)
{
	if (op1->sign != op2->sign) {
		return op1->sign > op2->sign ? 1 : -1;
	}

	if (op1->sign == 0) {
		return 0;
	}

	int sign = op1->sign;

	if (op1->size != op2->size) {
		return (op1->size - op2->size) * sign;
	}

	for (int i = op1->size - 1; i >= 0; i--) {
		if (op1->limbs[i] > op2->limbs[i]) {
			return 1 * sign;
		} else if (op1->limbs[i] < op2->limbs[i]) {
			return -1 * sign;
		}
	}

	return 0;
}

int limber_cmp_mag(Limber op1, Limber op2)
{
	if (op1->size > op2->size ) {
		return 1;
	} else if (op1->size < op2->size) {
		return -1;
	} else {
		int i = op1->size - 1;
		while (i >= 0) {
			if (op1->limbs[i] > op2->limbs[i]) {
				return 1;
			} else if (op1->limbs[i] < op2->limbs [i]) {
				return -1;
			}
			i--;
		}

		return 0;
	}
}

static void limber_add_mag(Limber rop, Limber op1, Limber op2)
{
	limb_t carry = 0;
	for (int i = 0; i < op2->size; i++) {
		limb_t sum = op1->limbs[i] + op2->limbs[i];
		limb_t carry_out = sum < op1->limbs[i];
		sum += carry;
		carry_out = carry_out || (sum < op2->limbs[i]);
		carry = carry_out;
		rop->limbs[i] = sum;
	}

	if (carry == 0) {
		memcpy(&rop->limbs[op2->size], &op1->limbs[op2->size], (op1->size - op2->size) * sizeof(limb_t));
		rop->size = op1->size;
		return;
	}

	for (int i = op2->size; i < op1->size; i++) {
		limb_t sum = op1->limbs[i] + carry;
		carry = sum < op1->limbs[i];
		rop->limbs[i] = sum;
	}

	if (carry) {
		rop->limbs[op1->size] = carry;
		rop->size = op1->size + 1;
	} else {
		rop->size = op1->size;
	}
}

static void limber_sub_mag(Limber rop, Limber op1, Limber op2)
{
	limb_t borrow = 0;
	for (int i = 0; i < op2->size; i++) {
		limb_t borrow_in = borrow;
		limb_t borrow_out = op1->limbs[i] < borrow_in;

		limb_t a = op1->limbs[i] - borrow_in;
		limb_t sub = a - op2->limbs[i];
		borrow_out |= a < op2->limbs[i];
		rop->limbs[i] = sub;

		borrow = borrow_out;
	}

	if (borrow == 0) {
		memcpy(&rop->limbs[op2->size], &op1->limbs[op2->size], (op1->size - op2->size) * sizeof(limb_t));
		rop->size = op1->size;
		return;
	}

	for (int i = op2->size; i < op1->size; i++) {
		limb_t sub = op1->limbs[i] - borrow;
		borrow = op1->limbs[i] < borrow;
		rop->limbs[i] = sub;
	}

	rop->size = op1->size;
}

void limber_add(Limber rop, Limber op1, Limber op2)
{
	Limber temp;
	limber_init(temp);

	int r = limber_cmp_mag(op1, op2);
	int needed_size = (op1->size > op2->size) ? op1->size : op2->size;
	needed_size++;
	limber_realloc(temp, needed_size);

	if (op1->sign == op2->sign) {
		temp->sign = op1->sign;

		if (r >= 0) {
			limber_add_mag(temp, op1, op2);
		} else {
			limber_add_mag(temp, op2, op1);
		}
	} else {
		if (r == 1) {
			limber_sub_mag(temp, op1, op2);
			temp->sign = op1->sign;
		} else if (r == -1) {
			limber_sub_mag(temp, op2, op1);
			temp->sign = op2->sign;
		} else if (r == 0) {
			//TODO(): Rop is zero
		}
	}

	limber_normalize(temp);
	limber_swap(rop, temp);
	limber_clear(temp);
}

void limber_sub(Limber rop, Limber op1, Limber op2)
{
	Limber temp;
	limber_init(temp);

	int needed_size = (op1->size > op2->size) ? op1->size : op2->size;
	needed_size++;
	limber_realloc(temp, needed_size);

	int r = limber_cmp_mag(op1, op2);

	if (op1->sign == op2->sign) {
		if (r == 1) {
			limber_sub_mag(temp, op1, op2);
			temp->sign = op1->sign;
		} else if (r == -1) {
			limber_sub_mag(temp, op2, op1);
			temp->sign = -(op1->sign);
		} else {
			//TODO(): Rop is zero
		}
	} else {
		temp->sign = op1->sign;

		if (r >= 0) {
			limber_add_mag(temp, op1, op2);
		} else {
			limber_add_mag(temp, op2, op1);
		}
	}

	limber_normalize(temp);
	limber_swap(rop, temp);
	limber_clear(temp);
}

void limber_mul(Limber rop, Limber op1, Limber op2)
{
	Limber temp;
	limber_init(temp);
	limber_realloc(temp, op1->size + op2->size + 1);
	temp->sign = op1->sign * op2->sign;

	for (int i = 0; i < op1->size; i++) {
		for (int j = 0; j < op2->size; j++) {
			limb_t mul_hi = 0;
			limb_t mul_lo = 0;

			limber_mul_single(&mul_hi, &mul_lo, op1->limbs[i], op2->limbs[i]);
			temp->limbs[i + j + 1] = mul_hi;
			temp->limbs[i + j] = mul_lo;
		}
	}

	limber_normalize(temp);
	limber_swap(rop, temp);
	limber_clear(temp);
}

int limber_sizeinbase(Limber l, int base)
{
	int last_limb_index = l->size - 1;
	limb_t last_limb = l->limbs[last_limb_index];
	int bits_in_last_limb = BITS_PER_LIMB - limber_clz(last_limb);
	int total_bits = (last_limb_index * BITS_PER_LIMB) + bits_in_last_limb;

	double log2_base = log2((double)base);
	double estimated_digits = (double)total_bits / log2_base;

	return ceil(estimated_digits);
}

static const char *LIMBER_DIGITS_MAP = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";

static void limber_str_reverse(char *start, char *end)
{
	char temp;
	while (start < end) {
		temp = *start;
		*start = *end;
		*end = temp;
		start++;
		end--;
	}
}


char *limber_to_str(Limber l, int base)
{
	int digits_estimated = limber_sizeinbase(l, base);
	int buffer_size = digits_estimated + 2;
	char *str = malloc(sizeof(char) * buffer_size);
	if (str == NULL) {
		return NULL;
	}

	if (limber_is_zero(l)) {
		str[0] = '0';
		str[1] = '\0';
		return str;
	}

	Limber temp;
	limber_init(temp);
	limber_set(temp, l);

	char *p = str;
	char *digits_start = str;

	if (temp->sign == -1) {
		*p++ = '-';
		digits_start++;
	}

	limb_t remainder;
	while (!limber_is_zero(temp)) {
		limber_div_limb(temp, &remainder, temp, (limb_t)base);
		*p++ = LIMBER_DIGITS_MAP[remainder];
	}
	*p = '\0';

	limber_str_reverse(digits_start, p - 1);
	limber_clear(temp);
	return str;
}
