/* 
* Inference of empirical DDF using Hill function.
* ------------------------------------------------------------------------------
* Author: Griffin Chure
* License: MIT
*
* Description
* ------------------------------------------------------------------------------
* This model infers the difference in free energy of binding by computing the 
* empirical delta F. The induction profile is fit using a Hill function from
* which the fold-change is used to compute the Delta F. The reason for this is
* that the Hill function is restricted to [0, 1], unlike the experimental data,
* and (hopefully) gives a more realistic picture of the empirical delta F.
*/
