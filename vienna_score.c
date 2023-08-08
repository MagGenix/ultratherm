#include <stdlib.h>
#include <stdio.h>
#include <string.h>
 
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/part_func.h>

float score(char *seq, unsigned char score_region[]) {
  /* allocate memory for pairing propensity string (length + 1) */
  char      *propensity = (char *)vrna_alloc(sizeof(char) * (strlen(seq) + 1));
 
  /* pointers for storing and navigating through base pair probabilities */
  vrna_ep_t *ptr, *pair_probabilities = NULL;
 
  float     en = vrna_pf_fold(seq, propensity, &pair_probabilities);
 
  /* print sequence, pairing propensity string and ensemble free energy */
  printf("%s\n%s [ %6.2f ]\n", seq, propensity, en);
  printf("%s\n", score_region);
  /* print all base pairs with probability above 50% */
  for (ptr = pair_probabilities; ptr->i != 0; ptr++)
    if (ptr->p > 0.5)
      printf("p(%d, %d) = %g\n", ptr->i, ptr->j, ptr->p);
  
  /* cleanup memory */
  free(pair_probabilities);
  free(propensity);
  return 0;

} 

int main() {
  char *seq = "GATTACA";
  unsigned char score_region[] = "0110101";
  score(seq, score_region);
}