#ifndef __c1_Ali_model_h__
#define __c1_Ali_model_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc1_Ali_modelInstanceStruct
#define typedef_SFc1_Ali_modelInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c1_sfEvent;
  boolean_T c1_isStable;
  boolean_T c1_doneDoubleBufferReInit;
  uint8_T c1_is_active_c1_Ali_model;
  real_T (*c1_kp)[16];
  real_T (*c1_kd)[16];
  real_T (*c1_e)[4];
  real_T (*c1_e_dot)[4];
  real_T (*c1_q)[8];
  real_T (*c1_Tau)[4];
  real_T (*c1_q_dot)[8];
} SFc1_Ali_modelInstanceStruct;

#endif                                 /*typedef_SFc1_Ali_modelInstanceStruct*/

/* Named Constants */

/* Variable Declarations */
extern struct SfDebugInstanceStruct *sfGlobalDebugInstanceStruct;

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c1_Ali_model_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c1_Ali_model_get_check_sum(mxArray *plhs[]);
extern void c1_Ali_model_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
