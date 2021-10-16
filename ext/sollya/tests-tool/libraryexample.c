#include <mpfi.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

int foo(mpfi_t rop, mpfi_t op, int n) {
  (void) n;
  (void) op;
  mpfi_interv_ui(rop, 0, 1);
  return 0;
}

/* Externalproc designed to test non-regression of bug #21283 */
int takeyourtime() {
  if (fork() == 0) {
    sleep(5);
    printf("After\n");
    exit(0);
  }
  return 1;
}
