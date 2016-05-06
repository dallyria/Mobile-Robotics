/* Include files */

#include <stddef.h>
#include "blas.h"
#include "Ali_model_sfun.h"
#include "c3_Ali_model.h"
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "Ali_model_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static real_T _sfTime_;
static const char * c3_debug_family_names[14] = { "tau_l1", "tau_r1", "tau_l2",
  "tau_r2", "theta_1", "theta_2", "theta_dot1", "theta_dot2", "qdoubledot",
  "nargin", "nargout", "q", "tau", "qdot" };

static const char * c3_b_debug_family_names[83] = { "t2", "t3", "t4", "t5", "t6",
  "t7", "t8", "t9", "t10", "t11", "t12", "t13", "t14", "t15", "t16", "t17",
  "t18", "t19", "t25", "t26", "t27", "t28", "t20", "t21", "t22", "t23", "t24",
  "t29", "t30", "t31", "t32", "t33", "t61", "t62", "t34", "t35", "t36", "t52",
  "t37", "t38", "t39", "t40", "t41", "t42", "t43", "t44", "t56", "t57", "t45",
  "t46", "t47", "t48", "t49", "t50", "t63", "t64", "t51", "t53", "t54", "t55",
  "t58", "t59", "t60", "t65", "t66", "t67", "t68", "t69", "t70", "t71", "t72",
  "t73", "nargin", "nargout", "tau_l1", "tau_l2", "tau_r1", "tau_r2", "theta_1",
  "theta_2", "theta_dot1", "theta_dot2", "qdotdot" };

/* Function Declarations */
static void initialize_c3_Ali_model(SFc3_Ali_modelInstanceStruct *chartInstance);
static void initialize_params_c3_Ali_model(SFc3_Ali_modelInstanceStruct
  *chartInstance);
static void enable_c3_Ali_model(SFc3_Ali_modelInstanceStruct *chartInstance);
static void disable_c3_Ali_model(SFc3_Ali_modelInstanceStruct *chartInstance);
static void c3_update_debugger_state_c3_Ali_model(SFc3_Ali_modelInstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c3_Ali_model(SFc3_Ali_modelInstanceStruct
  *chartInstance);
static void set_sim_state_c3_Ali_model(SFc3_Ali_modelInstanceStruct
  *chartInstance, const mxArray *c3_st);
static void finalize_c3_Ali_model(SFc3_Ali_modelInstanceStruct *chartInstance);
static void sf_gateway_c3_Ali_model(SFc3_Ali_modelInstanceStruct *chartInstance);
static void mdl_start_c3_Ali_model(SFc3_Ali_modelInstanceStruct *chartInstance);
static void initSimStructsc3_Ali_model(SFc3_Ali_modelInstanceStruct
  *chartInstance);
static void c3_Qdotdot(SFc3_Ali_modelInstanceStruct *chartInstance, real_T
  c3_tau_l1, real_T c3_tau_l2, real_T c3_tau_r1, real_T c3_tau_r2, real_T
  c3_theta_1, real_T c3_theta_2, real_T c3_theta_dot1, real_T c3_theta_dot2,
  real_T c3_qdotdot[4]);
static void init_script_number_translation(uint32_T c3_machineNumber, uint32_T
  c3_chartNumber, uint32_T c3_instanceNumber);
static const mxArray *c3_sf_marshallOut(void *chartInstanceVoid, void *c3_inData);
static void c3_emlrt_marshallIn(SFc3_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c3_b_qdot, const char_T *c3_identifier, real_T c3_y[8]);
static void c3_b_emlrt_marshallIn(SFc3_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[8]);
static void c3_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_b_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static const mxArray *c3_c_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static real_T c3_c_emlrt_marshallIn(SFc3_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId);
static void c3_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static void c3_d_emlrt_marshallIn(SFc3_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[4]);
static void c3_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static void c3_info_helper(const mxArray **c3_info);
static const mxArray *c3_emlrt_marshallOut(const char * c3_u);
static const mxArray *c3_b_emlrt_marshallOut(const uint32_T c3_u);
static real_T c3_power(SFc3_Ali_modelInstanceStruct *chartInstance, real_T c3_a);
static void c3_eml_scalar_eg(SFc3_Ali_modelInstanceStruct *chartInstance);
static real_T c3_rdivide(SFc3_Ali_modelInstanceStruct *chartInstance, real_T
  c3_x, real_T c3_y);
static const mxArray *c3_d_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static int32_T c3_e_emlrt_marshallIn(SFc3_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId);
static void c3_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static uint8_T c3_f_emlrt_marshallIn(SFc3_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c3_b_is_active_c3_Ali_model, const char_T *c3_identifier);
static uint8_T c3_g_emlrt_marshallIn(SFc3_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId);
static void init_dsm_address_info(SFc3_Ali_modelInstanceStruct *chartInstance);
static void init_simulink_io_address(SFc3_Ali_modelInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c3_Ali_model(SFc3_Ali_modelInstanceStruct *chartInstance)
{
  chartInstance->c3_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c3_is_active_c3_Ali_model = 0U;
}

static void initialize_params_c3_Ali_model(SFc3_Ali_modelInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void enable_c3_Ali_model(SFc3_Ali_modelInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c3_Ali_model(SFc3_Ali_modelInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c3_update_debugger_state_c3_Ali_model(SFc3_Ali_modelInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c3_Ali_model(SFc3_Ali_modelInstanceStruct
  *chartInstance)
{
  const mxArray *c3_st;
  const mxArray *c3_y = NULL;
  int32_T c3_i0;
  real_T c3_u[8];
  const mxArray *c3_b_y = NULL;
  uint8_T c3_hoistedGlobal;
  uint8_T c3_b_u;
  const mxArray *c3_c_y = NULL;
  c3_st = NULL;
  c3_st = NULL;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_createcellmatrix(2, 1), false);
  for (c3_i0 = 0; c3_i0 < 8; c3_i0++) {
    c3_u[c3_i0] = (*chartInstance->c3_qdot)[c3_i0];
  }

  c3_b_y = NULL;
  sf_mex_assign(&c3_b_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 1, 8), false);
  sf_mex_setcell(c3_y, 0, c3_b_y);
  c3_hoistedGlobal = chartInstance->c3_is_active_c3_Ali_model;
  c3_b_u = c3_hoistedGlobal;
  c3_c_y = NULL;
  sf_mex_assign(&c3_c_y, sf_mex_create("y", &c3_b_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c3_y, 1, c3_c_y);
  sf_mex_assign(&c3_st, c3_y, false);
  return c3_st;
}

static void set_sim_state_c3_Ali_model(SFc3_Ali_modelInstanceStruct
  *chartInstance, const mxArray *c3_st)
{
  const mxArray *c3_u;
  real_T c3_dv0[8];
  int32_T c3_i1;
  chartInstance->c3_doneDoubleBufferReInit = true;
  c3_u = sf_mex_dup(c3_st);
  c3_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c3_u, 0)), "qdot",
                      c3_dv0);
  for (c3_i1 = 0; c3_i1 < 8; c3_i1++) {
    (*chartInstance->c3_qdot)[c3_i1] = c3_dv0[c3_i1];
  }

  chartInstance->c3_is_active_c3_Ali_model = c3_f_emlrt_marshallIn(chartInstance,
    sf_mex_dup(sf_mex_getcell(c3_u, 1)), "is_active_c3_Ali_model");
  sf_mex_destroy(&c3_u);
  c3_update_debugger_state_c3_Ali_model(chartInstance);
  sf_mex_destroy(&c3_st);
}

static void finalize_c3_Ali_model(SFc3_Ali_modelInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c3_Ali_model(SFc3_Ali_modelInstanceStruct *chartInstance)
{
  int32_T c3_i2;
  int32_T c3_i3;
  real_T c3_b_q[8];
  int32_T c3_i4;
  real_T c3_b_tau[4];
  uint32_T c3_debug_family_var_map[14];
  real_T c3_tau_l1;
  real_T c3_tau_r1;
  real_T c3_tau_l2;
  real_T c3_tau_r2;
  real_T c3_theta_1;
  real_T c3_theta_2;
  real_T c3_theta_dot1;
  real_T c3_theta_dot2;
  real_T c3_qdoubledot[4];
  real_T c3_nargin = 2.0;
  real_T c3_nargout = 1.0;
  real_T c3_b_qdot[8];
  int32_T c3_i5;
  real_T c3_dv1[4];
  int32_T c3_i6;
  int32_T c3_i7;
  int32_T c3_i8;
  int32_T c3_i9;
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 2U, chartInstance->c3_sfEvent);
  for (c3_i2 = 0; c3_i2 < 8; c3_i2++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c3_q)[c3_i2], 0U);
  }

  chartInstance->c3_sfEvent = CALL_EVENT;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 2U, chartInstance->c3_sfEvent);
  for (c3_i3 = 0; c3_i3 < 8; c3_i3++) {
    c3_b_q[c3_i3] = (*chartInstance->c3_q)[c3_i3];
  }

  for (c3_i4 = 0; c3_i4 < 4; c3_i4++) {
    c3_b_tau[c3_i4] = (*chartInstance->c3_tau)[c3_i4];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 14U, 14U, c3_debug_family_names,
    c3_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_tau_l1, 0U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_tau_r1, 1U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_tau_l2, 2U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_tau_r2, 3U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_theta_1, 4U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_theta_2, 5U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_theta_dot1, 6U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_theta_dot2, 7U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_qdoubledot, 8U, c3_b_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargin, 9U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargout, 10U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c3_b_q, 11U, c3_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c3_b_tau, 12U, c3_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_b_qdot, 13U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 4);
  for (c3_i5 = 0; c3_i5 < 8; c3_i5++) {
    c3_b_qdot[c3_i5] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 10);
  c3_tau_l1 = c3_b_tau[0];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 11);
  c3_tau_r1 = c3_b_tau[1];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 12);
  c3_tau_l2 = c3_b_tau[2];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 13);
  c3_tau_r2 = c3_b_tau[3];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 17);
  c3_theta_1 = c3_b_q[2];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 18);
  c3_theta_2 = c3_b_q[3];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 19);
  c3_theta_dot1 = c3_b_q[6];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 20);
  c3_theta_dot2 = c3_b_q[7];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 24);
  c3_Qdotdot(chartInstance, c3_tau_l1, c3_tau_l2, c3_tau_r1, c3_tau_r2,
             c3_theta_1, c3_theta_2, c3_theta_dot1, c3_theta_dot2, c3_dv1);
  for (c3_i6 = 0; c3_i6 < 4; c3_i6++) {
    c3_qdoubledot[c3_i6] = c3_dv1[c3_i6];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 28);
  c3_b_qdot[0] = c3_b_q[4];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 29);
  c3_b_qdot[1] = c3_b_q[5];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 30);
  c3_b_qdot[2] = c3_b_q[6];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 31);
  c3_b_qdot[3] = c3_b_q[7];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 32);
  c3_b_qdot[4] = c3_qdoubledot[0];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 33);
  c3_b_qdot[5] = c3_qdoubledot[1];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 34);
  c3_b_qdot[6] = c3_qdoubledot[2];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 35);
  c3_b_qdot[7] = c3_qdoubledot[3];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, -35);
  _SFD_SYMBOL_SCOPE_POP();
  for (c3_i7 = 0; c3_i7 < 8; c3_i7++) {
    (*chartInstance->c3_qdot)[c3_i7] = c3_b_qdot[c3_i7];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 2U, chartInstance->c3_sfEvent);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_Ali_modelMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c3_i8 = 0; c3_i8 < 8; c3_i8++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c3_qdot)[c3_i8], 1U);
  }

  for (c3_i9 = 0; c3_i9 < 4; c3_i9++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c3_tau)[c3_i9], 2U);
  }
}

static void mdl_start_c3_Ali_model(SFc3_Ali_modelInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void initSimStructsc3_Ali_model(SFc3_Ali_modelInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void c3_Qdotdot(SFc3_Ali_modelInstanceStruct *chartInstance, real_T
  c3_tau_l1, real_T c3_tau_l2, real_T c3_tau_r1, real_T c3_tau_r2, real_T
  c3_theta_1, real_T c3_theta_2, real_T c3_theta_dot1, real_T c3_theta_dot2,
  real_T c3_qdotdot[4])
{
  uint32_T c3_debug_family_var_map[83];
  real_T c3_t2;
  real_T c3_t3;
  real_T c3_t4;
  real_T c3_t5;
  real_T c3_t6;
  real_T c3_t7;
  real_T c3_t8;
  real_T c3_t9;
  real_T c3_t10;
  real_T c3_t11;
  real_T c3_t12;
  real_T c3_t13;
  real_T c3_t14;
  real_T c3_t15;
  real_T c3_t16;
  real_T c3_t17;
  real_T c3_t18;
  real_T c3_t19;
  real_T c3_t25;
  real_T c3_t26;
  real_T c3_t27;
  real_T c3_t28;
  real_T c3_t20;
  real_T c3_t21;
  real_T c3_t22;
  real_T c3_t23;
  real_T c3_t24;
  real_T c3_t29;
  real_T c3_t30;
  real_T c3_t31;
  real_T c3_t32;
  real_T c3_t33;
  real_T c3_t61;
  real_T c3_t62;
  real_T c3_t34;
  real_T c3_t35;
  real_T c3_t36;
  real_T c3_t52;
  real_T c3_t37;
  real_T c3_t38;
  real_T c3_t39;
  real_T c3_t40;
  real_T c3_t41;
  real_T c3_t42;
  real_T c3_t43;
  real_T c3_t44;
  real_T c3_t56;
  real_T c3_t57;
  real_T c3_t45;
  real_T c3_t46;
  real_T c3_t47;
  real_T c3_t48;
  real_T c3_t49;
  real_T c3_t50;
  real_T c3_t63;
  real_T c3_t64;
  real_T c3_t51;
  real_T c3_t53;
  real_T c3_t54;
  real_T c3_t55;
  real_T c3_t58;
  real_T c3_t59;
  real_T c3_t60;
  real_T c3_t65;
  real_T c3_t66;
  real_T c3_t67;
  real_T c3_t68;
  real_T c3_t69;
  real_T c3_t70;
  real_T c3_t71;
  real_T c3_t72;
  real_T c3_t73;
  real_T c3_nargin = 8.0;
  real_T c3_nargout = 1.0;
  real_T c3_x;
  real_T c3_b_x;
  real_T c3_c_x;
  real_T c3_d_x;
  real_T c3_e_x;
  real_T c3_f_x;
  real_T c3_g_x;
  real_T c3_h_x;
  real_T c3_i_x;
  real_T c3_j_x;
  real_T c3_k_x;
  real_T c3_l_x;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 83U, 83U, c3_b_debug_family_names,
    c3_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t2, 0U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t3, 1U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t4, 2U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t5, 3U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t6, 4U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t7, 5U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t8, 6U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t9, 7U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t10, 8U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t11, 9U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t12, 10U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t13, 11U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t14, 12U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t15, 13U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t16, 14U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t17, 15U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t18, 16U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t19, 17U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t25, 18U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t26, 19U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t27, 20U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t28, 21U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t20, 22U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t21, 23U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t22, 24U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t23, 25U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t24, 26U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t29, 27U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t30, 28U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t31, 29U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t32, 30U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t33, 31U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t61, 32U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t62, 33U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t34, 34U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t35, 35U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t36, 36U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t52, 37U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t37, 38U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t38, 39U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t39, 40U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t40, 41U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t41, 42U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t42, 43U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t43, 44U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t44, 45U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t56, 46U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t57, 47U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t45, 48U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t46, 49U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t47, 50U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t48, 51U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t49, 52U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t50, 53U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t63, 54U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t64, 55U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t51, 56U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t53, 57U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t54, 58U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t55, 59U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t58, 60U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t59, 61U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t60, 62U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t65, 63U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t66, 64U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t67, 65U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t68, 66U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t69, 67U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t70, 68U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t71, 69U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t72, 70U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t73, 71U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargin, 72U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargout, 73U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_tau_l1, 74U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_tau_l2, 75U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_tau_r1, 76U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_tau_r2, 77U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_theta_1, 78U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_theta_2, 79U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_theta_dot1, 80U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_theta_dot2, 81U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_qdotdot, 82U, c3_b_sf_marshallOut,
    c3_c_sf_marshallIn);
  CV_SCRIPT_FCN(0, 0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 8);
  c3_t2 = c3_theta_1 - c3_theta_2;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 9);
  c3_x = c3_theta_1;
  c3_t3 = c3_x;
  c3_b_x = c3_t3;
  c3_t3 = c3_b_x;
  c3_t3 = muDoubleScalarCos(c3_t3);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 10);
  c3_c_x = c3_theta_2;
  c3_t4 = c3_c_x;
  c3_d_x = c3_t4;
  c3_t4 = c3_d_x;
  c3_t4 = muDoubleScalarSin(c3_t4);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 11);
  c3_e_x = c3_theta_1;
  c3_t5 = c3_e_x;
  c3_f_x = c3_t5;
  c3_t5 = c3_f_x;
  c3_t5 = muDoubleScalarSin(c3_t5);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 12);
  c3_t6 = c3_power(chartInstance, c3_t3);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 13);
  c3_g_x = c3_theta_2;
  c3_t7 = c3_g_x;
  c3_h_x = c3_t7;
  c3_t7 = c3_h_x;
  c3_t7 = muDoubleScalarCos(c3_t7);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 14);
  c3_i_x = c3_t2;
  c3_t8 = c3_i_x;
  c3_j_x = c3_t8;
  c3_t8 = c3_j_x;
  c3_t8 = muDoubleScalarCos(c3_t8);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 15);
  c3_t9 = c3_power(chartInstance, c3_t4);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 16);
  c3_t10 = c3_power(chartInstance, c3_t7);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 17);
  c3_t11 = c3_power(chartInstance, c3_t5);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 18);
  c3_t12 = c3_power(chartInstance, c3_t8);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 19);
  c3_t13 = c3_power(chartInstance, c3_theta_dot1);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 20);
  c3_t14 = c3_t6 * 1.819875895536978E+43;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 21);
  c3_t15 = c3_t10 * 1.819875895536978E+43;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 22);
  c3_t16 = c3_t11 * 1.819875895536978E+43;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 23);
  c3_t17 = c3_t9 * 1.819875895536978E+43;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 24);
  c3_t18 = c3_t12 * 2.8573578913191159E+43;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 25);
  c3_t19 = c3_t3 * c3_t4 * c3_t5 * c3_t7 * 1.4286789456595579E+43;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 26);
  c3_t25 = c3_t6 * c3_t9 * 7.1433947282977922E+42;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 27);
  c3_t26 = c3_t10 * c3_t11 * 7.1433947282977922E+42;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 28);
  c3_t27 = c3_t3 * c3_t7 * c3_t8 * 2.8573578913191169E+43;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 29);
  c3_t28 = c3_t4 * c3_t5 * c3_t8 * 2.8573578913191169E+43;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 30);
  c3_t20 = (((((((((c3_t14 + c3_t15) + c3_t16) + c3_t17) + c3_t18) + c3_t19) -
               c3_t25) - c3_t26) - c3_t27) - c3_t28) - 4.636378642267367E+43;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 31);
  c3_t21 = c3_rdivide(chartInstance, 1.0, c3_t20);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 32);
  c3_t22 = c3_power(chartInstance, c3_theta_dot2);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 33);
  c3_k_x = c3_t2;
  c3_t23 = c3_k_x;
  c3_l_x = c3_t23;
  c3_t23 = c3_l_x;
  c3_t23 = muDoubleScalarSin(c3_t23);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 34);
  c3_t24 = c3_t23 * 0.1258;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 35);
  c3_t29 = c3_t12 * 1.4629672403553879E+37;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 36);
  c3_t30 = c3_t5 * c3_tau_l1;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 37);
  c3_t31 = c3_t4 * c3_tau_l2;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 38);
  c3_t32 = c3_t5 * c3_tau_r1;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 39);
  c3_t33 = c3_t4 * c3_tau_r2;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 40);
  c3_t61 = c3_t5 * c3_t13 * 0.07548;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 41);
  c3_t62 = c3_t4 * c3_t22 * 0.07548;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 42);
  c3_t34 = ((((c3_t30 + c3_t31) + c3_t32) + c3_t33) - c3_t61) - c3_t62;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 43);
  c3_t35 = c3_tau_r2 * c3_rdivide(chartInstance, 17.0, 200.0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 44);
  c3_t36 = c3_t13 * c3_t23 * 0.009495384;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 45);
  c3_t52 = c3_tau_l2 * c3_rdivide(chartInstance, 17.0, 200.0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 46);
  c3_t37 = (c3_t35 + c3_t36) - c3_t52;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 47);
  c3_t38 = c3_t3 * c3_t5 * 3.486257080019396E+15;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 48);
  c3_t39 = c3_t4 * c3_t7 * 3.486257080019396E+15;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 49);
  c3_t40 = ((c3_t38 + c3_t39) - c3_t3 * c3_t4 * c3_t8 * 2.736858102024945E+15) -
    c3_t5 * c3_t7 * c3_t8 * 2.736858102024945E+15;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 50);
  c3_t41 = c3_t3 * c3_tau_l1;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 51);
  c3_t42 = c3_t7 * c3_tau_l2;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 52);
  c3_t43 = c3_t3 * c3_tau_r1;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 53);
  c3_t44 = c3_t7 * c3_tau_r2;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 54);
  c3_t56 = c3_t3 * c3_t13 * 0.07548;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 55);
  c3_t57 = c3_t7 * c3_t22 * 0.07548;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 56);
  c3_t45 = ((((c3_t41 + c3_t42) + c3_t43) + c3_t44) - c3_t56) - c3_t57;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 57);
  c3_t46 = c3_tau_r1 * c3_rdivide(chartInstance, 17.0, 200.0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 58);
  c3_t47 = c3_t24 - c3_rdivide(chartInstance, 17.0, 200.0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 59);
  c3_t48 = c3_t47 * c3_tau_l2;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 60);
  c3_t49 = c3_t24 + c3_rdivide(chartInstance, 17.0, 200.0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 61);
  c3_t50 = c3_t49 * c3_tau_r2;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 62);
  c3_t63 = c3_tau_l1 * c3_rdivide(chartInstance, 17.0, 200.0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 63);
  c3_t64 = c3_t22 * c3_t23 * 0.009495384;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 64);
  c3_t51 = (((c3_t46 + c3_t48) + c3_t50) - c3_t63) - c3_t64;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 65);
  c3_t53 = c3_t5 * 6.8090958594128835E+21;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 66);
  c3_t54 = c3_t3 * c3_t4 * c3_t7 * 2.6727129902587361E+21;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 67);
  c3_t55 = ((c3_t53 + c3_t54) - c3_t4 * c3_t8 * 5.3454259805174712E+21) - c3_t5 *
    c3_t10 * 2.6727129902587361E+21;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 68);
  c3_t58 = c3_t3 * 6.8090958594128835E+21;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 69);
  c3_t59 = c3_t4 * c3_t5 * c3_t7 * 2.6727129902587361E+21;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 70);
  c3_t60 = ((c3_t58 + c3_t59) - c3_t3 * c3_t9 * 2.6727129902587361E+21) - c3_t7 *
    c3_t8 * 5.3454259805174712E+21;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 71);
  c3_t65 = c3_t4 * 6.8090958594128835E+21;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 72);
  c3_t66 = c3_t3 * c3_t5 * c3_t7 * 2.6727129902587361E+21;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 73);
  c3_t67 = ((c3_t65 + c3_t66) - c3_t4 * c3_t6 * 2.6727129902587361E+21) - c3_t5 *
    c3_t8 * 5.3454259805174712E+21;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 74);
  c3_t68 = c3_t3 * c3_t7 * 8.909043300862452E+20;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 75);
  c3_t69 = c3_t4 * c3_t5 * 8.909043300862452E+20;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 76);
  c3_t70 = (c3_t8 * -1.78180866017249E+21 + c3_t68) + c3_t69;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 77);
  c3_t71 = c3_t7 * 6.8090958594128835E+21;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 78);
  c3_t72 = c3_t3 * c3_t4 * c3_t5 * 2.6727129902587361E+21;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 79);
  c3_t73 = ((c3_t71 + c3_t72) - c3_t3 * c3_t8 * 5.3454259805174712E+21) - c3_t7 *
    c3_t11 * 2.6727129902587361E+21;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 80);
  c3_qdotdot[0] = ((c3_t21 * c3_t45 * ((((c3_t6 * 9.3177645851493275E+36 +
    c3_t10 * 9.3177645851493275E+36) + c3_t29) - c3_t3 * c3_t7 * c3_t8 *
    1.4629672403553879E+37) - 2.373825864840892E+37) * 1.627604166666667E+6 +
                    c3_t21 * c3_t34 * c3_t40 * 4.3501187992492442E+27) + c3_t21 *
                   c3_t37 * c3_t67 * 3.5409552070200521E+22) + c3_t21 * c3_t51 *
    c3_t55 * 3.5409552070200521E+22;
  c3_qdotdot[1] = ((c3_t21 * c3_t34 * ((((c3_t9 * 9.3177645851493275E+36 +
    c3_t11 * 9.3177645851493275E+36) + c3_t29) - c3_t4 * c3_t5 * c3_t8 *
    1.4629672403553879E+37) - 2.373825864840892E+37) * 1.627604166666667E+6 +
                    c3_t21 * c3_t40 * c3_t45 * 4.3501187992492442E+27) - c3_t21 *
                   c3_t37 * c3_t73 * 3.5409552070200521E+22) - c3_t21 * c3_t51 *
    c3_t60 * 3.5409552070200521E+22;
  c3_qdotdot[2] = ((c3_t21 * c3_t51 * ((c3_t9 * 6.681782475646839E+20 + c3_t10 *
    6.681782475646839E+20) - 1.7022739648532211E+21) * 2.251799813685248E+24 -
                    c3_t21 * c3_t34 * c3_t60 * 3.5409552070200521E+22) + c3_t21 *
                   c3_t45 * c3_t55 * 3.5409552070200521E+22) - c3_t21 * c3_t37 *
    c3_t70 * 1.688849860263936E+24;
  c3_qdotdot[3] = ((c3_t21 * c3_t37 * ((c3_t6 * 6.681782475646839E+20 + c3_t11 *
    6.681782475646839E+20) - 1.7022739648532211E+21) * 2.251799813685248E+24 -
                    c3_t21 * c3_t34 * c3_t73 * 3.5409552070200521E+22) + c3_t21 *
                   c3_t45 * c3_t67 * 3.5409552070200521E+22) - c3_t21 * c3_t51 *
    c3_t70 * 1.688849860263936E+24;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, -80);
  _SFD_SYMBOL_SCOPE_POP();
}

static void init_script_number_translation(uint32_T c3_machineNumber, uint32_T
  c3_chartNumber, uint32_T c3_instanceNumber)
{
  (void)c3_machineNumber;
  _SFD_SCRIPT_TRANSLATION(c3_chartNumber, c3_instanceNumber, 0U,
    sf_debug_get_script_id(
    "C:\\Users\\nlashkar\\Google Drive\\F15\\Multi_robots_control\\4. Computed Torque Control\\Ali\\Qdotdot.m"));
}

static const mxArray *c3_sf_marshallOut(void *chartInstanceVoid, void *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i10;
  real_T c3_b_inData[8];
  int32_T c3_i11;
  real_T c3_u[8];
  const mxArray *c3_y = NULL;
  SFc3_Ali_modelInstanceStruct *chartInstance;
  chartInstance = (SFc3_Ali_modelInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  for (c3_i10 = 0; c3_i10 < 8; c3_i10++) {
    c3_b_inData[c3_i10] = (*(real_T (*)[8])c3_inData)[c3_i10];
  }

  for (c3_i11 = 0; c3_i11 < 8; c3_i11++) {
    c3_u[c3_i11] = c3_b_inData[c3_i11];
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 1, 8), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static void c3_emlrt_marshallIn(SFc3_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c3_b_qdot, const char_T *c3_identifier, real_T c3_y[8])
{
  emlrtMsgIdentifier c3_thisId;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_b_qdot), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_b_qdot);
}

static void c3_b_emlrt_marshallIn(SFc3_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[8])
{
  real_T c3_dv2[8];
  int32_T c3_i12;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv2, 1, 0, 0U, 1, 0U, 1, 8);
  for (c3_i12 = 0; c3_i12 < 8; c3_i12++) {
    c3_y[c3_i12] = c3_dv2[c3_i12];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_b_qdot;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[8];
  int32_T c3_i13;
  SFc3_Ali_modelInstanceStruct *chartInstance;
  chartInstance = (SFc3_Ali_modelInstanceStruct *)chartInstanceVoid;
  c3_b_qdot = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_b_qdot), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_b_qdot);
  for (c3_i13 = 0; c3_i13 < 8; c3_i13++) {
    (*(real_T (*)[8])c3_outData)[c3_i13] = c3_y[c3_i13];
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_b_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i14;
  real_T c3_b_inData[4];
  int32_T c3_i15;
  real_T c3_u[4];
  const mxArray *c3_y = NULL;
  SFc3_Ali_modelInstanceStruct *chartInstance;
  chartInstance = (SFc3_Ali_modelInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  for (c3_i14 = 0; c3_i14 < 4; c3_i14++) {
    c3_b_inData[c3_i14] = (*(real_T (*)[4])c3_inData)[c3_i14];
  }

  for (c3_i15 = 0; c3_i15 < 4; c3_i15++) {
    c3_u[c3_i15] = c3_b_inData[c3_i15];
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 1, 4), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static const mxArray *c3_c_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  real_T c3_u;
  const mxArray *c3_y = NULL;
  SFc3_Ali_modelInstanceStruct *chartInstance;
  chartInstance = (SFc3_Ali_modelInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_u = *(real_T *)c3_inData;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", &c3_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static real_T c3_c_emlrt_marshallIn(SFc3_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId)
{
  real_T c3_y;
  real_T c3_d0;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), &c3_d0, 1, 0, 0U, 0, 0U, 0);
  c3_y = c3_d0;
  sf_mex_destroy(&c3_u);
  return c3_y;
}

static void c3_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_nargout;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y;
  SFc3_Ali_modelInstanceStruct *chartInstance;
  chartInstance = (SFc3_Ali_modelInstanceStruct *)chartInstanceVoid;
  c3_nargout = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_y = c3_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_nargout), &c3_thisId);
  sf_mex_destroy(&c3_nargout);
  *(real_T *)c3_outData = c3_y;
  sf_mex_destroy(&c3_mxArrayInData);
}

static void c3_d_emlrt_marshallIn(SFc3_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[4])
{
  real_T c3_dv3[4];
  int32_T c3_i16;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv3, 1, 0, 0U, 1, 0U, 1, 4);
  for (c3_i16 = 0; c3_i16 < 4; c3_i16++) {
    c3_y[c3_i16] = c3_dv3[c3_i16];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_qdoubledot;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[4];
  int32_T c3_i17;
  SFc3_Ali_modelInstanceStruct *chartInstance;
  chartInstance = (SFc3_Ali_modelInstanceStruct *)chartInstanceVoid;
  c3_qdoubledot = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_qdoubledot), &c3_thisId,
                        c3_y);
  sf_mex_destroy(&c3_qdoubledot);
  for (c3_i17 = 0; c3_i17 < 4; c3_i17++) {
    (*(real_T (*)[4])c3_outData)[c3_i17] = c3_y[c3_i17];
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

const mxArray *sf_c3_Ali_model_get_eml_resolved_functions_info(void)
{
  const mxArray *c3_nameCaptureInfo = NULL;
  c3_nameCaptureInfo = NULL;
  sf_mex_assign(&c3_nameCaptureInfo, sf_mex_createstruct("structure", 2, 20, 1),
                false);
  c3_info_helper(&c3_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c3_nameCaptureInfo);
  return c3_nameCaptureInfo;
}

static void c3_info_helper(const mxArray **c3_info)
{
  const mxArray *c3_rhs0 = NULL;
  const mxArray *c3_lhs0 = NULL;
  const mxArray *c3_rhs1 = NULL;
  const mxArray *c3_lhs1 = NULL;
  const mxArray *c3_rhs2 = NULL;
  const mxArray *c3_lhs2 = NULL;
  const mxArray *c3_rhs3 = NULL;
  const mxArray *c3_lhs3 = NULL;
  const mxArray *c3_rhs4 = NULL;
  const mxArray *c3_lhs4 = NULL;
  const mxArray *c3_rhs5 = NULL;
  const mxArray *c3_lhs5 = NULL;
  const mxArray *c3_rhs6 = NULL;
  const mxArray *c3_lhs6 = NULL;
  const mxArray *c3_rhs7 = NULL;
  const mxArray *c3_lhs7 = NULL;
  const mxArray *c3_rhs8 = NULL;
  const mxArray *c3_lhs8 = NULL;
  const mxArray *c3_rhs9 = NULL;
  const mxArray *c3_lhs9 = NULL;
  const mxArray *c3_rhs10 = NULL;
  const mxArray *c3_lhs10 = NULL;
  const mxArray *c3_rhs11 = NULL;
  const mxArray *c3_lhs11 = NULL;
  const mxArray *c3_rhs12 = NULL;
  const mxArray *c3_lhs12 = NULL;
  const mxArray *c3_rhs13 = NULL;
  const mxArray *c3_lhs13 = NULL;
  const mxArray *c3_rhs14 = NULL;
  const mxArray *c3_lhs14 = NULL;
  const mxArray *c3_rhs15 = NULL;
  const mxArray *c3_lhs15 = NULL;
  const mxArray *c3_rhs16 = NULL;
  const mxArray *c3_lhs16 = NULL;
  const mxArray *c3_rhs17 = NULL;
  const mxArray *c3_lhs17 = NULL;
  const mxArray *c3_rhs18 = NULL;
  const mxArray *c3_lhs18 = NULL;
  const mxArray *c3_rhs19 = NULL;
  const mxArray *c3_lhs19 = NULL;
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("Qdotdot"), "name", "name", 0);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[E]C:/Users/nlashkar/Google Drive/F15/Multi_robots_control/4. Computed Torque Control/Ali/Qdotdot.m"),
                  "resolved", "resolved", 0);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1445640296U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c3_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[E]C:/Users/nlashkar/Google Drive/F15/Multi_robots_control/4. Computed Torque Control/Ali/Qdotdot.m"),
                  "context", "context", 1);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("cos"), "name", "name", 1);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 1);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1395346496U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c3_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 2);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 2);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 2);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286840322U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c3_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[E]C:/Users/nlashkar/Google Drive/F15/Multi_robots_control/4. Computed Torque Control/Ali/Qdotdot.m"),
                  "context", "context", 3);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("sin"), "name", "name", 3);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 3);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1395346504U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c3_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 4);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286840336U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c3_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[E]C:/Users/nlashkar/Google Drive/F15/Multi_robots_control/4. Computed Torque Control/Ali/Qdotdot.m"),
                  "context", "context", 5);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("power"), "name", "name", 5);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "resolved",
                  "resolved", 5);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1395346506U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c3_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 6);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1395949856U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c3_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 7);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 7);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1376002288U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c3_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 8);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1410829370U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c3_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 9);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 9);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1376002288U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c3_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "context", "context", 10);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.scalexpAlloc"),
                  "name", "name", 10);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1410829370U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c3_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 11);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("floor"), "name", "name", 11);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "resolved",
                  "resolved", 11);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363731854U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c3_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 12);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 12);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1395949856U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c3_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 13);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 13);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286840326U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c3_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 14);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 14);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 14);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1376002288U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c3_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[E]C:/Users/nlashkar/Google Drive/F15/Multi_robots_control/4. Computed Torque Control/Ali/Qdotdot.m"),
                  "context", "context", 15);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("rdivide"), "name", "name", 15);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 15);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363731880U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c3_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 16);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 16);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1395949856U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c3_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 17);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 17);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 17);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286840396U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c3_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 18);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_div"), "name", "name", 18);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 18);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1386445552U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c3_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 19);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 19);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 19);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1410829370U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c3_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs19), "lhs", "lhs",
                  19);
  sf_mex_destroy(&c3_rhs0);
  sf_mex_destroy(&c3_lhs0);
  sf_mex_destroy(&c3_rhs1);
  sf_mex_destroy(&c3_lhs1);
  sf_mex_destroy(&c3_rhs2);
  sf_mex_destroy(&c3_lhs2);
  sf_mex_destroy(&c3_rhs3);
  sf_mex_destroy(&c3_lhs3);
  sf_mex_destroy(&c3_rhs4);
  sf_mex_destroy(&c3_lhs4);
  sf_mex_destroy(&c3_rhs5);
  sf_mex_destroy(&c3_lhs5);
  sf_mex_destroy(&c3_rhs6);
  sf_mex_destroy(&c3_lhs6);
  sf_mex_destroy(&c3_rhs7);
  sf_mex_destroy(&c3_lhs7);
  sf_mex_destroy(&c3_rhs8);
  sf_mex_destroy(&c3_lhs8);
  sf_mex_destroy(&c3_rhs9);
  sf_mex_destroy(&c3_lhs9);
  sf_mex_destroy(&c3_rhs10);
  sf_mex_destroy(&c3_lhs10);
  sf_mex_destroy(&c3_rhs11);
  sf_mex_destroy(&c3_lhs11);
  sf_mex_destroy(&c3_rhs12);
  sf_mex_destroy(&c3_lhs12);
  sf_mex_destroy(&c3_rhs13);
  sf_mex_destroy(&c3_lhs13);
  sf_mex_destroy(&c3_rhs14);
  sf_mex_destroy(&c3_lhs14);
  sf_mex_destroy(&c3_rhs15);
  sf_mex_destroy(&c3_lhs15);
  sf_mex_destroy(&c3_rhs16);
  sf_mex_destroy(&c3_lhs16);
  sf_mex_destroy(&c3_rhs17);
  sf_mex_destroy(&c3_lhs17);
  sf_mex_destroy(&c3_rhs18);
  sf_mex_destroy(&c3_lhs18);
  sf_mex_destroy(&c3_rhs19);
  sf_mex_destroy(&c3_lhs19);
}

static const mxArray *c3_emlrt_marshallOut(const char * c3_u)
{
  const mxArray *c3_y = NULL;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c3_u)), false);
  return c3_y;
}

static const mxArray *c3_b_emlrt_marshallOut(const uint32_T c3_u)
{
  const mxArray *c3_y = NULL;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", &c3_u, 7, 0U, 0U, 0U, 0), false);
  return c3_y;
}

static real_T c3_power(SFc3_Ali_modelInstanceStruct *chartInstance, real_T c3_a)
{
  real_T c3_b_a;
  real_T c3_ak;
  real_T c3_c_a;
  c3_b_a = c3_a;
  c3_eml_scalar_eg(chartInstance);
  c3_ak = c3_b_a;
  c3_c_a = c3_ak;
  c3_eml_scalar_eg(chartInstance);
  return c3_c_a * c3_c_a;
}

static void c3_eml_scalar_eg(SFc3_Ali_modelInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static real_T c3_rdivide(SFc3_Ali_modelInstanceStruct *chartInstance, real_T
  c3_x, real_T c3_y)
{
  real_T c3_b_x;
  real_T c3_b_y;
  real_T c3_c_x;
  real_T c3_c_y;
  (void)chartInstance;
  c3_b_x = c3_x;
  c3_b_y = c3_y;
  c3_c_x = c3_b_x;
  c3_c_y = c3_b_y;
  return c3_c_x / c3_c_y;
}

static const mxArray *c3_d_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_u;
  const mxArray *c3_y = NULL;
  SFc3_Ali_modelInstanceStruct *chartInstance;
  chartInstance = (SFc3_Ali_modelInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_u = *(int32_T *)c3_inData;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", &c3_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static int32_T c3_e_emlrt_marshallIn(SFc3_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId)
{
  int32_T c3_y;
  int32_T c3_i18;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), &c3_i18, 1, 6, 0U, 0, 0U, 0);
  c3_y = c3_i18;
  sf_mex_destroy(&c3_u);
  return c3_y;
}

static void c3_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_b_sfEvent;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  int32_T c3_y;
  SFc3_Ali_modelInstanceStruct *chartInstance;
  chartInstance = (SFc3_Ali_modelInstanceStruct *)chartInstanceVoid;
  c3_b_sfEvent = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_y = c3_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_b_sfEvent),
    &c3_thisId);
  sf_mex_destroy(&c3_b_sfEvent);
  *(int32_T *)c3_outData = c3_y;
  sf_mex_destroy(&c3_mxArrayInData);
}

static uint8_T c3_f_emlrt_marshallIn(SFc3_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c3_b_is_active_c3_Ali_model, const char_T *c3_identifier)
{
  uint8_T c3_y;
  emlrtMsgIdentifier c3_thisId;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_y = c3_g_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c3_b_is_active_c3_Ali_model), &c3_thisId);
  sf_mex_destroy(&c3_b_is_active_c3_Ali_model);
  return c3_y;
}

static uint8_T c3_g_emlrt_marshallIn(SFc3_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId)
{
  uint8_T c3_y;
  uint8_T c3_u0;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), &c3_u0, 1, 3, 0U, 0, 0U, 0);
  c3_y = c3_u0;
  sf_mex_destroy(&c3_u);
  return c3_y;
}

static void init_dsm_address_info(SFc3_Ali_modelInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void init_simulink_io_address(SFc3_Ali_modelInstanceStruct *chartInstance)
{
  chartInstance->c3_q = (real_T (*)[8])ssGetInputPortSignal_wrapper
    (chartInstance->S, 0);
  chartInstance->c3_qdot = (real_T (*)[8])ssGetOutputPortSignal_wrapper
    (chartInstance->S, 1);
  chartInstance->c3_tau = (real_T (*)[4])ssGetInputPortSignal_wrapper
    (chartInstance->S, 1);
}

/* SFunction Glue Code */
#ifdef utFree
#undef utFree
#endif

#ifdef utMalloc
#undef utMalloc
#endif

#ifdef __cplusplus

extern "C" void *utMalloc(size_t size);
extern "C" void utFree(void*);

#else

extern void *utMalloc(size_t size);
extern void utFree(void*);

#endif

void sf_c3_Ali_model_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3781048134U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1769498609U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3329649307U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(604600189U);
}

mxArray* sf_c3_Ali_model_get_post_codegen_info(void);
mxArray *sf_c3_Ali_model_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals", "postCodegenInfo" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1, 1, sizeof
    (autoinheritanceFields)/sizeof(autoinheritanceFields[0]),
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("rw0gTy6bB2ja2XHY3DORmF");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,2,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(8);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(4);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(8);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  {
    mxArray* mxPostCodegenInfo = sf_c3_Ali_model_get_post_codegen_info();
    mxSetField(mxAutoinheritanceInfo,0,"postCodegenInfo",mxPostCodegenInfo);
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c3_Ali_model_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c3_Ali_model_jit_fallback_info(void)
{
  const char *infoFields[] = { "fallbackType", "fallbackReason",
    "incompatibleSymbol", };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 3, infoFields);
  mxArray *fallbackReason = mxCreateString("feature_off");
  mxArray *incompatibleSymbol = mxCreateString("");
  mxArray *fallbackType = mxCreateString("early");
  mxSetField(mxInfo, 0, infoFields[0], fallbackType);
  mxSetField(mxInfo, 0, infoFields[1], fallbackReason);
  mxSetField(mxInfo, 0, infoFields[2], incompatibleSymbol);
  return mxInfo;
}

mxArray *sf_c3_Ali_model_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

mxArray* sf_c3_Ali_model_get_post_codegen_info(void)
{
  const char* fieldNames[] = { "exportedFunctionsUsedByThisChart",
    "exportedFunctionsChecksum" };

  mwSize dims[2] = { 1, 1 };

  mxArray* mxPostCodegenInfo = mxCreateStructArray(2, dims, sizeof(fieldNames)/
    sizeof(fieldNames[0]), fieldNames);

  {
    mxArray* mxExportedFunctionsChecksum = mxCreateString("");
    mwSize exp_dims[2] = { 0, 1 };

    mxArray* mxExportedFunctionsUsedByThisChart = mxCreateCellArray(2, exp_dims);
    mxSetField(mxPostCodegenInfo, 0, "exportedFunctionsUsedByThisChart",
               mxExportedFunctionsUsedByThisChart);
    mxSetField(mxPostCodegenInfo, 0, "exportedFunctionsChecksum",
               mxExportedFunctionsChecksum);
  }

  return mxPostCodegenInfo;
}

static const mxArray *sf_get_sim_state_info_c3_Ali_model(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"qdot\",},{M[8],M[0],T\"is_active_c3_Ali_model\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c3_Ali_model_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc3_Ali_modelInstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc3_Ali_modelInstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _Ali_modelMachineNumber_,
           3,
           1,
           1,
           0,
           3,
           0,
           0,
           0,
           0,
           1,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           (void *)S);

        /* Each instance must initialize its own list of scripts */
        init_script_number_translation(_Ali_modelMachineNumber_,
          chartInstance->chartNumber,chartInstance->instanceNumber);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_Ali_modelMachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _Ali_modelMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,1,1,0,"q");
          _SFD_SET_DATA_PROPS(1,2,0,1,"qdot");
          _SFD_SET_DATA_PROPS(2,1,1,0,"tau");
          _SFD_STATE_INFO(0,0,2);
          _SFD_CH_SUBSTATE_COUNT(0);
          _SFD_CH_SUBSTATE_DECOMP(0);
        }

        _SFD_CV_INIT_CHART(0,0,0,0);

        {
          _SFD_CV_INIT_STATE(0,0,0,0,0,0,NULL,NULL);
        }

        _SFD_CV_INIT_TRANS(0,0,NULL,NULL,0,NULL);

        /* Initialization of MATLAB Function Model Coverage */
        _SFD_CV_INIT_EML(0,1,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,453);
        _SFD_CV_INIT_SCRIPT(0,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_SCRIPT_FCN(0,0,"Qdotdot",0,-1,3431);

        {
          unsigned int dimVector[1];
          dimVector[0]= 8;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 8;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)
            c3_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c3_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_VALUE_PTR(0U, *chartInstance->c3_q);
        _SFD_SET_DATA_VALUE_PTR(1U, *chartInstance->c3_qdot);
        _SFD_SET_DATA_VALUE_PTR(2U, *chartInstance->c3_tau);
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _Ali_modelMachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "L8yfMi2qPpifF1yCkuVosE";
}

static void sf_opaque_initialize_c3_Ali_model(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc3_Ali_modelInstanceStruct*) chartInstanceVar
    )->S,0);
  initialize_params_c3_Ali_model((SFc3_Ali_modelInstanceStruct*)
    chartInstanceVar);
  initialize_c3_Ali_model((SFc3_Ali_modelInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c3_Ali_model(void *chartInstanceVar)
{
  enable_c3_Ali_model((SFc3_Ali_modelInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c3_Ali_model(void *chartInstanceVar)
{
  disable_c3_Ali_model((SFc3_Ali_modelInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c3_Ali_model(void *chartInstanceVar)
{
  sf_gateway_c3_Ali_model((SFc3_Ali_modelInstanceStruct*) chartInstanceVar);
}

static const mxArray* sf_opaque_get_sim_state_c3_Ali_model(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  return get_sim_state_c3_Ali_model((SFc3_Ali_modelInstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
}

static void sf_opaque_set_sim_state_c3_Ali_model(SimStruct* S, const mxArray *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  set_sim_state_c3_Ali_model((SFc3_Ali_modelInstanceStruct*)
    chartInfo->chartInstance, st);
}

static void sf_opaque_terminate_c3_Ali_model(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc3_Ali_modelInstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_Ali_model_optimization_info();
    }

    finalize_c3_Ali_model((SFc3_Ali_modelInstanceStruct*) chartInstanceVar);
    utFree(chartInstanceVar);
    if (crtInfo != NULL) {
      utFree(crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc3_Ali_model((SFc3_Ali_modelInstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c3_Ali_model(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    initialize_params_c3_Ali_model((SFc3_Ali_modelInstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c3_Ali_model(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_Ali_model_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,3);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,3,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,3,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,3);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,3,2);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,3,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 2; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,3);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(2991611473U));
  ssSetChecksum1(S,(300441850U));
  ssSetChecksum2(S,(3142670611U));
  ssSetChecksum3(S,(2902339027U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c3_Ali_model(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c3_Ali_model(SimStruct *S)
{
  SFc3_Ali_modelInstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc3_Ali_modelInstanceStruct *)utMalloc(sizeof
    (SFc3_Ali_modelInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc3_Ali_modelInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c3_Ali_model;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c3_Ali_model;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c3_Ali_model;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c3_Ali_model;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c3_Ali_model;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c3_Ali_model;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c3_Ali_model;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c3_Ali_model;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c3_Ali_model;
  chartInstance->chartInfo.mdlStart = mdlStart_c3_Ali_model;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c3_Ali_model;
  chartInstance->chartInfo.extModeExec = NULL;
  chartInstance->chartInfo.restoreLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.restoreBeforeLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.storeCurrentConfiguration = NULL;
  chartInstance->chartInfo.callAtomicSubchartUserFcn = NULL;
  chartInstance->chartInfo.callAtomicSubchartAutoFcn = NULL;
  chartInstance->chartInfo.debugInstance = sfGlobalDebugInstanceStruct;
  chartInstance->S = S;
  crtInfo->checksum = SF_RUNTIME_INFO_CHECKSUM;
  crtInfo->instanceInfo = (&(chartInstance->chartInfo));
  crtInfo->isJITEnabled = false;
  crtInfo->compiledInfo = NULL;
  ssSetUserData(S,(void *)(crtInfo));  /* register the chart instance with simstruct */
  init_dsm_address_info(chartInstance);
  init_simulink_io_address(chartInstance);
  if (!sim_mode_is_rtw_gen(S)) {
  }

  sf_opaque_init_subchart_simstructs(chartInstance->chartInfo.chartInstance);
  chart_debug_initialization(S,1);
}

void c3_Ali_model_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c3_Ali_model(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c3_Ali_model(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c3_Ali_model(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c3_Ali_model_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
