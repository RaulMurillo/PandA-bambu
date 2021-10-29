#include <sollya.h>
#include <inttypes.h>

/* All warning messages must occur twice in a row */
int message_occurrence_number = 0;
int last_message = 0;

/* Notice: we assume that a call to
   sollya_lib_get_constant_as_[mpz|uint64_array] emits at most one
   message. This is indeed the case at the time when these tests are written. A
   warning is emitted if this assumption is not true anymore. In such a case,
   the callback would have to be changed to a more sophisticated one, in order
   to stack the messages and compare the stacks of both calls.
*/
int catch_message_and_check(sollya_msg_t msg, void *data) {
  (void)data;
  message_occurrence_number++;
  switch ( message_occurrence_number ) {
  case 1:
    last_message = sollya_lib_get_msg_id(msg);
    break;
  case 2:
    if (sollya_lib_get_msg_id(msg) != last_message) {
      sollya_lib_printf("The previous message id was %d and the current message id is %d\n", last_message, sollya_lib_get_msg_id(msg));
    }
    break;
  default: sollya_lib_printf("Warning: more than one warning message is emitted by one of the sollya_lib_get_constant_as_* functions!\n");
  }

  return 0;
}

static inline void mpz_set_ui64(mpz_t b, uint64_t a) {
  mpz_import(b, 1, 1, sizeof(a), 0, 0, &a); return;
}

static inline void mpz_add_ui64(mpz_t res, mpz_t a, uint64_t b) {
  mpz_t tmp;
  mpz_init(tmp);
  mpz_set_ui64(tmp, b);
  mpz_add(res, a, tmp);
  mpz_clear(tmp);
  return;
}

void check_specif(sollya_obj_t arg) {
  int s = -17;
  uint64_t *value = NULL;
  size_t length = 17;
  int test1, test2;
  mpz_t val, tmp, B;
  int same_as_mpz = 1;
  int i;

  mpz_init(val);
  mpz_init(tmp);
  mpz_init(B);

  message_occurrence_number = 0;
  if (!sollya_lib_install_msg_callback(catch_message_and_check, NULL))
    printf("An error happened when installing the message callback\n");
  test1 = sollya_lib_get_constant_as_uint64_array(&s, &value, &length, arg);
  test2 = sollya_lib_get_constant_as_mpz(val, arg);
  sollya_lib_uninstall_msg_callback();

  if (!test1) {
    if (s != -17) sollya_lib_printf("Warning: sollya_lib_get_constant_as_uint64_array assigned a sign despite a failure!\n");
    if (length != 17) sollya_lib_printf("Warning: sollya_lib_get_constant_as_uint64_array assigned a length despite a failure!\n");
    if (value != NULL) sollya_lib_printf("Warning: sollya_lib_get_constant_as_uint64_array assigned an array despite a failure!\n");
  }
  else {
    if ( (s != -1) && (s != 0) && (s != 1) ) sollya_lib_printf("Warning: the returned sign is %d instead of an element of {-1, 0, 1}\n", s);
    if (length < 1) sollya_lib_printf("Warning: the returned length is %ud but should be greater or equal to 1\n", length);

    if (s) {
      if (value[length-1]==0) sollya_lib_printf("Warning: the last element of the array is 0\n");
      mpz_set_ui(B, 1); mpz_mul_2exp(B, B, 64);
      mpz_set_ui64(tmp, value[length-1]);
      for(i=length-2;i>=0;i--) {
        mpz_mul(tmp, tmp, B);
        mpz_add_ui64(tmp, tmp, value[i]);
      }
      if (s<0) mpz_neg(tmp, tmp);
    }
    else { mpz_set_ui(tmp, 0); }
  }

  if ((test1 && !test2) || (!test1 && test2)) {
    same_as_mpz = 0;
    sollya_lib_printf("The boolean value returned by sollya_lib_get_constant_as_uint64_array is %d whereas the value returned by sollya_lib_get_constant_as_mpz is %d\n", test1, test2);
  }
  if (mpz_cmp(tmp, val)) {
    same_as_mpz = 0;
    sollya_lib_printf("The value returned by sollya_lib_get_constant_as_uint64_array is ");
    gmp_printf("%Zd", tmp);
    sollya_lib_printf(" whereas the value returned by sollya_lib_get_constant_as_mpz is ");
    gmp_printf("%Zd", val);
    sollya_lib_printf(".\n");
  }

  if (same_as_mpz) {
    sollya_lib_printf("%b is handled similarly by sollya_lib_get_constant_as_uint64_array and sollya_lib_get_constant_as_mpz\n", arg);
  }
  sollya_lib_printf("\n");

  mpz_clear(val);
  mpz_clear(tmp);
  mpz_clear(B);
  if (test1) sollya_lib_free(value);
  return;
}

int main(void) {
  sollya_obj_t a, prec;
  int i;

  sollya_lib_init();

  prec = SOLLYA_CONST(20);
  sollya_lib_set_prec_and_print(prec);
  sollya_lib_clear_obj(prec);

  /* something that is obviously not a constant */
  a = sollya_lib_parse_string("[1;2]");
  check_specif(a);
  sollya_lib_clear_obj(a);


  /* something that could be interpreted as constant but that is not a */
  /* constant, strictly speaking */
  a = sollya_lib_parse_string("[1;1]");
  check_specif(a);
  sollya_lib_clear_obj(a);


  /* something that is constant, but it is not obvious */
  a = sollya_lib_parse_string("3*cos(2*x)/(2*sin(x)*cos(x))");
  check_specif(a);
  sollya_lib_clear_obj(a);


  /* An obvious constant */
  a = SOLLYA_CONST(3);
  check_specif(a);
  sollya_lib_clear_obj(a);

  /* A constant, but that does not fit on 20 bits but fits on 32 bits */
  a = SOLLYA_CONST(1073741824);
  check_specif(a);
  sollya_lib_clear_obj(a);

  /* A constant, that does not fit on 32 bits but fits on 64 bits */
  a = SOLLYA_CONST(1);
  for (i=1;i<=35;i++) a = SOLLYA_MUL(a, SOLLYA_CONST(2));
  a = SOLLYA_ADD(a, SOLLYA_CONST(17));
  check_specif(a);
  sollya_lib_clear_obj(a);

  /* A negative constant */
  a = SOLLYA_CONST(-3);
  check_specif(a);
  sollya_lib_clear_obj(a);

  /* A constant that does not fit in a integer */
  a = SOLLYA_CONST(1);
  for (i=1;i<=90;i++) a = SOLLYA_MUL(a, SOLLYA_CONST(2));
  check_specif(a);
  sollya_lib_clear_obj(a);


  /* The same, but negative */
  a = SOLLYA_CONST(-1);
  for (i=1;i<=90;i++) a = SOLLYA_MUL(a, SOLLYA_CONST(2));
  check_specif(a);
  sollya_lib_clear_obj(a);


  /* A constant expression exactly representable as an int but it cannot be decided. */
  a = sollya_lib_parse_string("(1b200+1)-1b200*(log2(3)/log2(7) - log(3)/log(7) + 1)");
  check_specif(a);
  sollya_lib_clear_obj(a);


  /* A constant expression very close to the middle between two integers
   and it cannot be decided. */
  a = sollya_lib_parse_string("1 + 1b-400 + 0.5*(log2(3)/log2(7) - log(3)/log(7) + 1)");
  check_specif(a);
  sollya_lib_clear_obj(a);


  /* A constant expression very close to the middle between two doubles */
  a = sollya_lib_parse_string("1 + 1b-400 + 0.5");
  check_specif(a);
  sollya_lib_clear_obj(a);


  /* A constant expression exactly at the middle between two doubles, but
     it cannot be decided. */
  a = sollya_lib_parse_string("1 - 0.5*(log2(3)/log2(7) - log(3)/log(7) + 1)");
  check_specif(a);
  sollya_lib_clear_obj(a);


  /* The same constant but decidable */
  a = sollya_lib_parse_string("1 - 0.5");
  check_specif(a);
  sollya_lib_clear_obj(a);


  /* Another constant expression exactly at the middle between two doubles, but
     it cannot be decided. */
  a = sollya_lib_parse_string("1 + 1.5*(log2(3)/log2(7) - log(3)/log(7) + 1)");
  check_specif(a);
  sollya_lib_clear_obj(a);


  /* The same constant but decidable. */
  a = sollya_lib_parse_string("1 + 1.5");
  check_specif(a);
  sollya_lib_clear_obj(a);


  /* A transcendantal constant. */
  a = sollya_lib_parse_string("exp(pi) + log(2)");
  check_specif(a);
  sollya_lib_clear_obj(a);


  /* A constant hard to evaluate because exactly zero but undecidable. */
  a = sollya_lib_parse_string("log10(2)/log10(3) - log(2)/log(3)");
  check_specif(a);
  sollya_lib_clear_obj(a);


  /* A constant very hard to evaluate (cf. tevaluate_function_at_constant_expression). */
  a = sollya_lib_parse_string("(sin((pi) / 3) - sqrt(3) / 2) * (1 * 2^(100000)) + 3");
  check_specif(a);
  sollya_lib_clear_obj(a);

  
  /* Trying -inf */
  a = sollya_lib_parse_string("-@Inf@");
  check_specif(a);
  sollya_lib_clear_obj(a);


/* Trying NaN */
  a = sollya_lib_parse_string("@NaN@");
  check_specif(a);
  sollya_lib_clear_obj(a);

  sollya_lib_close();
  return 0;
}
