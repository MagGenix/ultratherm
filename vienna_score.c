#include <stdlib.h>
#include <stdio.h>
#include <string.h>
 
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/utils/structures.h>
#include <ViennaRNA/dp_matrices.h>

float score(char *seq, unsigned char score_region[]) {

  /* allocate memory for structure */
  char *structure = (char *)space(sizeof(char) * (strlen(seq) + 1));
 
  /* create a new model details structure to store the Model Settings */
  vrna_md_t md;
  //float *pf;
  //FLT_OR_DBL *bppm;
  /* ALWAYS set default model settings first! */
  vrna_md_set_default(&md);

  /* change temperature and activate G-Quadruplex prediction */
  md.temperature  = 25.0; /* 25 Deg Celcius */
  md.gquad        = 1;    /* Turn-on G-Quadruples support */
  md.compute_bpp  = 1;    /* Compute base pair probabilities */
  md.dangles      = 2;    /* Compute stabilizing co-axial stacking */
  
  //NOTE: can pass in multiple sequences! must be separated by '&'
  vrna_fold_compound_t *fc = vrna_fold_compound(seq, &md, VRNA_OPTION_PF);
  vrna_mx_pf_add(fc, VRNA_MX_DEFAULT, VRNA_OPTION_HYBRID); // Add auxillary DP array for RNA:RNA hybrid fold

  FLT_OR_DBL en = vrna_pf(fc, structure);

  int* arr = vrna_idx_row_wise(strlen(seq));
  //int* arr2 = vrna_mx_pf_t_probs_get(fc->matrices);

  /* get the base pair probability matrix for the previous run of pf_fold() */

  /* print sequence, pairing propensity string and ensemble free energy */
  printf("%s\t%6.2f\n", seq, en);
  printf("%s\n", score_region);
  printf("%s\n", structure);
  //printf("%f\n", fc->matrices->type);
  //printf("%f\n", sizeof(fc->exp_matrices));
  /* print all base pairs with probability above 50% */
  /* cleanup memory */
  free(structure);
  //free(pf);
  //free(propensity);
  return 0;

} 

int main(void) {
  char *seq = "AAAGACUUCCUAAUAAGGAAAUCACAUUCGUGGCUUAAGGAGGUUCACCAUG";
  unsigned char score_region[] = "0110101";
  score(seq, score_region);
}
