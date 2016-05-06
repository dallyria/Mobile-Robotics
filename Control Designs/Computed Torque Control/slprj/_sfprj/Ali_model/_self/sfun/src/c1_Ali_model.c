/* Include files */

#include <stddef.h>
#include "blas.h"
#include "Ali_model_sfun.h"
#include "c1_Ali_model.h"
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
static const char * c1_debug_family_names[17] = { "theta_1", "theta_2",
  "theta_dot1", "theta_dot2", "x_dotdot", "y_dotdot", "theta_dotdot1",
  "theta_dotdot2", "nargin", "nargout", "kp", "kd", "e", "e_dot", "q", "q_dot",
  "Tau" };

static const char * c1_b_debug_family_names[45] = { "t2", "t3", "t4", "t5", "t6",
  "t7", "t8", "t9", "t10", "t11", "t18", "t12", "t13", "t14", "t15", "t16",
  "t17", "t19", "t20", "t21", "t22", "t33", "t34", "t23", "t24", "t25", "t26",
  "t27", "t28", "t29", "t30", "t31", "t32", "t35", "nargin", "nargout",
  "theta_1", "theta_2", "theta_dot1", "theta_dot2", "theta_dotdot1",
  "theta_dotdot2", "x_dotdot", "y_dotdot", "tau_computed" };

static const char * c1_c_debug_family_names[18] = { "t2", "t3", "t4", "t5", "t6",
  "t12", "t7", "t8", "t9", "t10", "t11", "t13", "t14", "nargin", "nargout",
  "theta_1", "theta_2", "out1" };

/* Function Declarations */
static void initialize_c1_Ali_model(SFc1_Ali_modelInstanceStruct *chartInstance);
static void initialize_params_c1_Ali_model(SFc1_Ali_modelInstanceStruct
  *chartInstance);
static void enable_c1_Ali_model(SFc1_Ali_modelInstanceStruct *chartInstance);
static void disable_c1_Ali_model(SFc1_Ali_modelInstanceStruct *chartInstance);
static void c1_update_debugger_state_c1_Ali_model(SFc1_Ali_modelInstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c1_Ali_model(SFc1_Ali_modelInstanceStruct
  *chartInstance);
static void set_sim_state_c1_Ali_model(SFc1_Ali_modelInstanceStruct
  *chartInstance, const mxArray *c1_st);
static void finalize_c1_Ali_model(SFc1_Ali_modelInstanceStruct *chartInstance);
static void sf_gateway_c1_Ali_model(SFc1_Ali_modelInstanceStruct *chartInstance);
static void mdl_start_c1_Ali_model(SFc1_Ali_modelInstanceStruct *chartInstance);
static void c1_chartstep_c1_Ali_model(SFc1_Ali_modelInstanceStruct
  *chartInstance);
static void initSimStructsc1_Ali_model(SFc1_Ali_modelInstanceStruct
  *chartInstance);
static void c1_B_inverse(SFc1_Ali_modelInstanceStruct *chartInstance, real_T
  c1_theta_1, real_T c1_theta_2, real_T c1_out1[16]);
static void init_script_number_translation(uint32_T c1_machineNumber, uint32_T
  c1_chartNumber, uint32_T c1_instanceNumber);
static const mxArray *c1_sf_marshallOut(void *chartInstanceVoid, void *c1_inData);
static void c1_emlrt_marshallIn(SFc1_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c1_b_Tau, const char_T *c1_identifier, real_T c1_y[4]);
static void c1_b_emlrt_marshallIn(SFc1_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[4]);
static void c1_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_b_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static const mxArray *c1_c_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static const mxArray *c1_d_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static real_T c1_c_emlrt_marshallIn(SFc1_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static void c1_d_emlrt_marshallIn(SFc1_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[16]);
static void c1_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static void c1_info_helper(const mxArray **c1_info);
static const mxArray *c1_emlrt_marshallOut(const char * c1_u);
static const mxArray *c1_b_emlrt_marshallOut(const uint32_T c1_u);
static real_T c1_power(SFc1_Ali_modelInstanceStruct *chartInstance, real_T c1_a);
static void c1_eml_scalar_eg(SFc1_Ali_modelInstanceStruct *chartInstance);
static real_T c1_rdivide(SFc1_Ali_modelInstanceStruct *chartInstance, real_T
  c1_x, real_T c1_y);
static void c1_b_eml_scalar_eg(SFc1_Ali_modelInstanceStruct *chartInstance);
static void c1_eml_xgemm(SFc1_Ali_modelInstanceStruct *chartInstance, real_T
  c1_A[16], real_T c1_B[4], real_T c1_C[4], real_T c1_b_C[4]);
static const mxArray *c1_e_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static int32_T c1_e_emlrt_marshallIn(SFc1_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static uint8_T c1_f_emlrt_marshallIn(SFc1_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c1_b_is_active_c1_Ali_model, const char_T *c1_identifier);
static uint8_T c1_g_emlrt_marshallIn(SFc1_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_b_eml_xgemm(SFc1_Ali_modelInstanceStruct *chartInstance, real_T
  c1_A[16], real_T c1_B[4], real_T c1_C[4]);
static void init_dsm_address_info(SFc1_Ali_modelInstanceStruct *chartInstance);
static void init_simulink_io_address(SFc1_Ali_modelInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c1_Ali_model(SFc1_Ali_modelInstanceStruct *chartInstance)
{
  chartInstance->c1_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c1_is_active_c1_Ali_model = 0U;
}

static void initialize_params_c1_Ali_model(SFc1_Ali_modelInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void enable_c1_Ali_model(SFc1_Ali_modelInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c1_Ali_model(SFc1_Ali_modelInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c1_update_debugger_state_c1_Ali_model(SFc1_Ali_modelInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c1_Ali_model(SFc1_Ali_modelInstanceStruct
  *chartInstance)
{
  const mxArray *c1_st;
  const mxArray *c1_y = NULL;
  int32_T c1_i0;
  real_T c1_u[4];
  const mxArray *c1_b_y = NULL;
  uint8_T c1_hoistedGlobal;
  uint8_T c1_b_u;
  const mxArray *c1_c_y = NULL;
  c1_st = NULL;
  c1_st = NULL;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_createcellmatrix(2, 1), false);
  for (c1_i0 = 0; c1_i0 < 4; c1_i0++) {
    c1_u[c1_i0] = (*chartInstance->c1_Tau)[c1_i0];
  }

  c1_b_y = NULL;
  sf_mex_assign(&c1_b_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 1, 4), false);
  sf_mex_setcell(c1_y, 0, c1_b_y);
  c1_hoistedGlobal = chartInstance->c1_is_active_c1_Ali_model;
  c1_b_u = c1_hoistedGlobal;
  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", &c1_b_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c1_y, 1, c1_c_y);
  sf_mex_assign(&c1_st, c1_y, false);
  return c1_st;
}

static void set_sim_state_c1_Ali_model(SFc1_Ali_modelInstanceStruct
  *chartInstance, const mxArray *c1_st)
{
  const mxArray *c1_u;
  real_T c1_dv0[4];
  int32_T c1_i1;
  chartInstance->c1_doneDoubleBufferReInit = true;
  c1_u = sf_mex_dup(c1_st);
  c1_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 0)), "Tau",
                      c1_dv0);
  for (c1_i1 = 0; c1_i1 < 4; c1_i1++) {
    (*chartInstance->c1_Tau)[c1_i1] = c1_dv0[c1_i1];
  }

  chartInstance->c1_is_active_c1_Ali_model = c1_f_emlrt_marshallIn(chartInstance,
    sf_mex_dup(sf_mex_getcell(c1_u, 1)), "is_active_c1_Ali_model");
  sf_mex_destroy(&c1_u);
  c1_update_debugger_state_c1_Ali_model(chartInstance);
  sf_mex_destroy(&c1_st);
}

static void finalize_c1_Ali_model(SFc1_Ali_modelInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c1_Ali_model(SFc1_Ali_modelInstanceStruct *chartInstance)
{
  int32_T c1_i2;
  int32_T c1_i3;
  int32_T c1_i4;
  int32_T c1_i5;
  int32_T c1_i6;
  int32_T c1_i7;
  int32_T c1_i8;
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
  for (c1_i2 = 0; c1_i2 < 16; c1_i2++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_kp)[c1_i2], 0U);
  }

  for (c1_i3 = 0; c1_i3 < 16; c1_i3++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_kd)[c1_i3], 1U);
  }

  for (c1_i4 = 0; c1_i4 < 4; c1_i4++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_e)[c1_i4], 2U);
  }

  for (c1_i5 = 0; c1_i5 < 4; c1_i5++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_e_dot)[c1_i5], 3U);
  }

  for (c1_i6 = 0; c1_i6 < 8; c1_i6++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_q)[c1_i6], 4U);
  }

  chartInstance->c1_sfEvent = CALL_EVENT;
  c1_chartstep_c1_Ali_model(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_Ali_modelMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c1_i7 = 0; c1_i7 < 4; c1_i7++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_Tau)[c1_i7], 5U);
  }

  for (c1_i8 = 0; c1_i8 < 8; c1_i8++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_q_dot)[c1_i8], 6U);
  }
}

static void mdl_start_c1_Ali_model(SFc1_Ali_modelInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c1_chartstep_c1_Ali_model(SFc1_Ali_modelInstanceStruct
  *chartInstance)
{
  int32_T c1_i9;
  real_T c1_b_kp[16];
  int32_T c1_i10;
  real_T c1_b_kd[16];
  int32_T c1_i11;
  real_T c1_b_e[4];
  int32_T c1_i12;
  real_T c1_b_e_dot[4];
  int32_T c1_i13;
  real_T c1_b_q[8];
  int32_T c1_i14;
  real_T c1_b_q_dot[8];
  uint32_T c1_debug_family_var_map[17];
  real_T c1_theta_1;
  real_T c1_theta_2;
  real_T c1_theta_dot1;
  real_T c1_theta_dot2;
  real_T c1_x_dotdot;
  real_T c1_y_dotdot;
  real_T c1_theta_dotdot1;
  real_T c1_theta_dotdot2;
  real_T c1_nargin = 6.0;
  real_T c1_nargout = 1.0;
  real_T c1_b_Tau[4];
  real_T c1_b_theta_1;
  real_T c1_b_theta_2;
  real_T c1_b_theta_dot1;
  real_T c1_b_theta_dot2;
  real_T c1_b_theta_dotdot1;
  real_T c1_b_theta_dotdot2;
  real_T c1_b_x_dotdot;
  real_T c1_b_y_dotdot;
  uint32_T c1_b_debug_family_var_map[45];
  real_T c1_t2;
  real_T c1_t3;
  real_T c1_t4;
  real_T c1_t5;
  real_T c1_t6;
  real_T c1_t7;
  real_T c1_t8;
  real_T c1_t9;
  real_T c1_t10;
  real_T c1_t11;
  real_T c1_t18;
  real_T c1_t12;
  real_T c1_t13;
  real_T c1_t14;
  real_T c1_t15;
  real_T c1_t16;
  real_T c1_t17;
  real_T c1_t19;
  real_T c1_t20;
  real_T c1_t21;
  real_T c1_t22;
  real_T c1_t33;
  real_T c1_t34;
  real_T c1_t23;
  real_T c1_t24;
  real_T c1_t25;
  real_T c1_t26;
  real_T c1_t27;
  real_T c1_t28;
  real_T c1_t29;
  real_T c1_t30;
  real_T c1_t31;
  real_T c1_t32;
  real_T c1_t35;
  real_T c1_b_nargin = 8.0;
  real_T c1_b_nargout = 1.0;
  real_T c1_x;
  real_T c1_b_x;
  real_T c1_c_x;
  real_T c1_d_x;
  real_T c1_e_x;
  real_T c1_f_x;
  real_T c1_g_x;
  real_T c1_h_x;
  real_T c1_i_x;
  real_T c1_j_x;
  real_T c1_k_x;
  real_T c1_l_x;
  real_T c1_d0;
  real_T c1_d1;
  real_T c1_d2;
  real_T c1_d3;
  real_T c1_d4;
  real_T c1_d5;
  real_T c1_d6;
  real_T c1_d7;
  real_T c1_d8;
  int32_T c1_i15;
  real_T c1_a[16];
  int32_T c1_i16;
  real_T c1_b[4];
  int32_T c1_i17;
  real_T c1_y[4];
  int32_T c1_i18;
  real_T c1_b_a[16];
  int32_T c1_i19;
  real_T c1_b_b[4];
  int32_T c1_i20;
  int32_T c1_i21;
  int32_T c1_i22;
  real_T c1_b_y[4];
  int32_T c1_i23;
  real_T c1_c_a[16];
  int32_T c1_i24;
  real_T c1_c_b[4];
  int32_T c1_i25;
  int32_T c1_i26;
  int32_T c1_i27;
  real_T c1_d_a[16];
  int32_T c1_i28;
  real_T c1_c_y[4];
  int32_T c1_i29;
  int32_T c1_i30;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
  for (c1_i9 = 0; c1_i9 < 16; c1_i9++) {
    c1_b_kp[c1_i9] = (*chartInstance->c1_kp)[c1_i9];
  }

  for (c1_i10 = 0; c1_i10 < 16; c1_i10++) {
    c1_b_kd[c1_i10] = (*chartInstance->c1_kd)[c1_i10];
  }

  for (c1_i11 = 0; c1_i11 < 4; c1_i11++) {
    c1_b_e[c1_i11] = (*chartInstance->c1_e)[c1_i11];
  }

  for (c1_i12 = 0; c1_i12 < 4; c1_i12++) {
    c1_b_e_dot[c1_i12] = (*chartInstance->c1_e_dot)[c1_i12];
  }

  for (c1_i13 = 0; c1_i13 < 8; c1_i13++) {
    c1_b_q[c1_i13] = (*chartInstance->c1_q)[c1_i13];
  }

  for (c1_i14 = 0; c1_i14 < 8; c1_i14++) {
    c1_b_q_dot[c1_i14] = (*chartInstance->c1_q_dot)[c1_i14];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 17U, 17U, c1_debug_family_names,
    c1_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_theta_1, 0U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_theta_2, 1U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_theta_dot1, 2U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_theta_dot2, 3U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_x_dotdot, 4U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_y_dotdot, 5U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_theta_dotdot1, 6U,
    c1_d_sf_marshallOut, c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_theta_dotdot2, 7U,
    c1_d_sf_marshallOut, c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargin, 8U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargout, 9U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_kp, 10U, c1_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_kd, 11U, c1_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_e, 12U, c1_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_e_dot, 13U, c1_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_q, 14U, c1_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_q_dot, 15U, c1_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_b_Tau, 16U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 3);
  c1_theta_1 = c1_b_q[2];
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 4);
  c1_theta_2 = c1_b_q[3];
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 6);
  c1_theta_dot1 = c1_b_q[6];
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 7);
  c1_theta_dot2 = c1_b_q[7];
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 9);
  c1_x_dotdot = c1_b_q_dot[4];
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 10);
  c1_y_dotdot = c1_b_q_dot[5];
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 11);
  c1_theta_dotdot1 = c1_b_q_dot[6];
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 12);
  c1_theta_dotdot2 = c1_b_q_dot[7];
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 26);
  c1_b_theta_1 = c1_theta_1;
  c1_b_theta_2 = c1_theta_2;
  c1_b_theta_dot1 = c1_theta_dot1;
  c1_b_theta_dot2 = c1_theta_dot2;
  c1_b_theta_dotdot1 = c1_theta_dotdot1;
  c1_b_theta_dotdot2 = c1_theta_dotdot2;
  c1_b_x_dotdot = c1_x_dotdot;
  c1_b_y_dotdot = c1_y_dotdot;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 45U, 45U, c1_b_debug_family_names,
    c1_b_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t2, 0U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t3, 1U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t4, 2U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t5, 3U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t6, 4U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t7, 5U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t8, 6U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t9, 7U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t10, 8U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t11, 9U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t18, 10U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t12, 11U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t13, 12U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t14, 13U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t15, 14U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t16, 15U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t17, 16U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t19, 17U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t20, 18U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t21, 19U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t22, 20U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t33, 21U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t34, 22U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t23, 23U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t24, 24U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t25, 25U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t26, 26U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t27, 27U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t28, 28U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t29, 29U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t30, 30U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t31, 31U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t32, 32U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t35, 33U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_nargin, 34U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_nargout, 35U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_theta_1, 36U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_theta_2, 37U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_theta_dot1, 38U,
    c1_d_sf_marshallOut, c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_theta_dot2, 39U,
    c1_d_sf_marshallOut, c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_theta_dotdot1, 40U,
    c1_d_sf_marshallOut, c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_theta_dotdot2, 41U,
    c1_d_sf_marshallOut, c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_x_dotdot, 42U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_y_dotdot, 43U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_b_Tau, 44U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  CV_SCRIPT_FCN(0, 0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 8);
  c1_t2 = c1_b_theta_1 - c1_b_theta_2;
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 9);
  c1_x = c1_t2;
  c1_t3 = c1_x;
  c1_b_x = c1_t3;
  c1_t3 = c1_b_x;
  c1_t3 = muDoubleScalarSin(c1_t3);
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 10);
  c1_c_x = c1_t2;
  c1_t4 = c1_c_x;
  c1_d_x = c1_t4;
  c1_t4 = c1_d_x;
  c1_t4 = muDoubleScalarCos(c1_t4);
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 11);
  c1_e_x = c1_b_theta_1;
  c1_t5 = c1_e_x;
  c1_f_x = c1_t5;
  c1_t5 = c1_f_x;
  c1_t5 = muDoubleScalarCos(c1_t5);
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 12);
  c1_g_x = c1_b_theta_2;
  c1_t6 = c1_g_x;
  c1_h_x = c1_t6;
  c1_t6 = c1_h_x;
  c1_t6 = muDoubleScalarSin(c1_t6);
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 13);
  c1_i_x = c1_b_theta_2;
  c1_t7 = c1_i_x;
  c1_j_x = c1_t7;
  c1_t7 = c1_j_x;
  c1_t7 = muDoubleScalarCos(c1_t7);
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 14);
  c1_k_x = c1_b_theta_1;
  c1_t8 = c1_k_x;
  c1_l_x = c1_t8;
  c1_t8 = c1_l_x;
  c1_t8 = muDoubleScalarSin(c1_t8);
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 15);
  c1_t9 = c1_power(chartInstance, c1_b_theta_dot1);
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 16);
  c1_t10 = c1_power(chartInstance, c1_b_theta_dot2);
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 17);
  c1_t11 = c1_t5 * c1_t6;
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 18);
  c1_t18 = c1_t7 * c1_t8;
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 19);
  c1_t12 = c1_t11 - c1_t18;
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 20);
  c1_t13 = c1_rdivide(chartInstance, 1.0, c1_t12);
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 21);
  c1_t14 = c1_b_theta_dotdot2 * 0.071149317647058827;
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 22);
  c1_t15 = c1_t5 * c1_b_y_dotdot * c1_rdivide(chartInstance, 111.0, 250.0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 23);
  c1_t16 = c1_t6 * c1_b_x_dotdot * c1_rdivide(chartInstance, 111.0, 250.0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 24);
  c1_t17 = c1_t4 * c1_b_theta_dotdot1 * 0.055855199999999987;
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 25);
  c1_t19 = c1_t7 * 25.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 26);
  c1_t20 = c1_b_y_dotdot * c1_rdivide(chartInstance, 6.0, 5.0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 27);
  c1_t21 = c1_t8 * c1_t9 * 0.07548;
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 28);
  c1_t22 = c1_t6 * c1_t10 * 0.07548;
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 29);
  c1_t33 = c1_t5 * c1_b_theta_dotdot1 * 0.07548;
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 30);
  c1_t34 = c1_t7 * c1_b_theta_dotdot2 * 0.07548;
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 31);
  c1_t23 = (((c1_t20 + c1_t21) + c1_t22) - c1_t33) - c1_t34;
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 32);
  c1_t24 = c1_t6 * 25.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 33);
  c1_t25 = c1_b_x_dotdot * c1_rdivide(chartInstance, 6.0, 5.0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 34);
  c1_t26 = c1_t5 * c1_t9 * 0.07548;
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 35);
  c1_t27 = c1_t7 * c1_t10 * 0.07548;
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 36);
  c1_t28 = c1_t8 * c1_b_theta_dotdot1 * 0.07548;
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 37);
  c1_t29 = c1_t6 * c1_b_theta_dotdot2 * 0.07548;
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 38);
  c1_t30 = (((c1_t25 + c1_t26) + c1_t27) + c1_t28) + c1_t29;
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 39);
  c1_t31 = c1_t3 * c1_t9 * 0.0558552;
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 40);
  c1_t32 = c1_t7 * c1_b_y_dotdot * c1_rdivide(chartInstance, 111.0, 250.0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 41);
  c1_t35 = c1_t5 * c1_t13 * c1_t23 * c1_rdivide(chartInstance, 1.0, 2.0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 42);
  c1_d0 = c1_rdivide(chartInstance, 111.0, 250.0);
  c1_d1 = c1_rdivide(chartInstance, 111.0, 250.0);
  c1_d2 = c1_rdivide(chartInstance, 1.0, 50.0);
  c1_d3 = c1_rdivide(chartInstance, 1.0, 50.0);
  c1_d4 = c1_rdivide(chartInstance, 111.0, 250.0);
  c1_d5 = c1_rdivide(chartInstance, 1.0, 50.0);
  c1_d6 = c1_rdivide(chartInstance, 1.0, 50.0);
  c1_d7 = c1_rdivide(chartInstance, 1.0, 2.0);
  c1_d8 = c1_rdivide(chartInstance, 1.0, 2.0);
  c1_b_Tau[0] = ((((((((((c1_t14 + c1_t15) + c1_t16) + c1_t17) -
                       c1_b_theta_dotdot1 * 0.071149317647058827) - c1_t3 *
                      c1_t9 * 0.0558552) - c1_t3 * c1_t10 * 0.0558552) - c1_t4 *
                    c1_b_theta_dotdot2 * 0.055855199999999987) - c1_t8 *
                   c1_b_x_dotdot * c1_d0) - c1_t7 * c1_b_y_dotdot * c1_d1) -
                 c1_t13 * c1_t23 * (c1_t19 - c1_t3 * c1_t5 * 37.0) * c1_d2) +
    c1_t13 * c1_t30 * (c1_t24 - c1_t3 * c1_t8 * 37.0) * c1_d3;
  c1_b_Tau[1] = ((((((((((-c1_t14 - c1_t15) - c1_t16) - c1_t17) + c1_t31) +
                      c1_t32) + c1_b_theta_dotdot1 * 0.071149317647058827) +
                    c1_t3 * c1_t10 * 0.0558552) + c1_t4 * c1_b_theta_dotdot2 *
                   0.055855199999999987) + c1_t8 * c1_b_x_dotdot * c1_d4) -
                 c1_t13 * c1_t23 * (c1_t19 + c1_t3 * c1_t5 * 37.0) * c1_d5) +
    c1_t13 * c1_t30 * (c1_t24 + c1_t3 * c1_t8 * 37.0) * c1_d6;
  c1_b_Tau[2] = (((((-c1_t14 - c1_t16) - c1_t17) + c1_t31) + c1_t32) + c1_t35) -
    c1_t8 * c1_t13 * c1_t30 * c1_d7;
  c1_b_Tau[3] = (((((c1_t14 + c1_t16) + c1_t17) - c1_t31) - c1_t32) + c1_t35) -
    c1_t8 * c1_t13 * c1_t30 * c1_d8;
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, -42);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 27);
  for (c1_i15 = 0; c1_i15 < 16; c1_i15++) {
    c1_a[c1_i15] = c1_b_kp[c1_i15];
  }

  for (c1_i16 = 0; c1_i16 < 4; c1_i16++) {
    c1_b[c1_i16] = c1_b_e[c1_i16];
  }

  c1_b_eml_scalar_eg(chartInstance);
  c1_b_eml_scalar_eg(chartInstance);
  for (c1_i17 = 0; c1_i17 < 4; c1_i17++) {
    c1_y[c1_i17] = 0.0;
  }

  for (c1_i18 = 0; c1_i18 < 16; c1_i18++) {
    c1_b_a[c1_i18] = c1_a[c1_i18];
  }

  for (c1_i19 = 0; c1_i19 < 4; c1_i19++) {
    c1_b_b[c1_i19] = c1_b[c1_i19];
  }

  c1_b_eml_xgemm(chartInstance, c1_b_a, c1_b_b, c1_y);
  for (c1_i20 = 0; c1_i20 < 16; c1_i20++) {
    c1_a[c1_i20] = c1_b_kd[c1_i20];
  }

  for (c1_i21 = 0; c1_i21 < 4; c1_i21++) {
    c1_b[c1_i21] = c1_b_e_dot[c1_i21];
  }

  c1_b_eml_scalar_eg(chartInstance);
  c1_b_eml_scalar_eg(chartInstance);
  for (c1_i22 = 0; c1_i22 < 4; c1_i22++) {
    c1_b_y[c1_i22] = 0.0;
  }

  for (c1_i23 = 0; c1_i23 < 16; c1_i23++) {
    c1_c_a[c1_i23] = c1_a[c1_i23];
  }

  for (c1_i24 = 0; c1_i24 < 4; c1_i24++) {
    c1_c_b[c1_i24] = c1_b[c1_i24];
  }

  c1_b_eml_xgemm(chartInstance, c1_c_a, c1_c_b, c1_b_y);
  c1_B_inverse(chartInstance, c1_theta_1, c1_theta_2, c1_a);
  for (c1_i25 = 0; c1_i25 < 4; c1_i25++) {
    c1_y[c1_i25] += c1_b_y[c1_i25];
  }

  c1_b_eml_scalar_eg(chartInstance);
  c1_b_eml_scalar_eg(chartInstance);
  for (c1_i26 = 0; c1_i26 < 4; c1_i26++) {
    c1_b_y[c1_i26] = 0.0;
  }

  for (c1_i27 = 0; c1_i27 < 16; c1_i27++) {
    c1_d_a[c1_i27] = c1_a[c1_i27];
  }

  for (c1_i28 = 0; c1_i28 < 4; c1_i28++) {
    c1_c_y[c1_i28] = c1_y[c1_i28];
  }

  c1_b_eml_xgemm(chartInstance, c1_d_a, c1_c_y, c1_b_y);
  for (c1_i29 = 0; c1_i29 < 4; c1_i29++) {
    c1_b_Tau[c1_i29] += c1_b_y[c1_i29];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, -27);
  _SFD_SYMBOL_SCOPE_POP();
  for (c1_i30 = 0; c1_i30 < 4; c1_i30++) {
    (*chartInstance->c1_Tau)[c1_i30] = c1_b_Tau[c1_i30];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
}

static void initSimStructsc1_Ali_model(SFc1_Ali_modelInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void c1_B_inverse(SFc1_Ali_modelInstanceStruct *chartInstance, real_T
  c1_theta_1, real_T c1_theta_2, real_T c1_out1[16])
{
  uint32_T c1_debug_family_var_map[18];
  real_T c1_t2;
  real_T c1_t3;
  real_T c1_t4;
  real_T c1_t5;
  real_T c1_t6;
  real_T c1_t12;
  real_T c1_t7;
  real_T c1_t8;
  real_T c1_t9;
  real_T c1_t10;
  real_T c1_t11;
  real_T c1_t13;
  real_T c1_t14;
  real_T c1_nargin = 2.0;
  real_T c1_nargout = 1.0;
  real_T c1_x;
  real_T c1_b_x;
  real_T c1_c_x;
  real_T c1_d_x;
  real_T c1_e_x;
  real_T c1_f_x;
  real_T c1_g_x;
  real_T c1_h_x;
  real_T c1_i_x;
  real_T c1_j_x;
  real_T c1_d9;
  real_T c1_d10;
  real_T c1_d11;
  real_T c1_d12;
  real_T c1_d13;
  real_T c1_d14;
  real_T c1_d15;
  real_T c1_d16;
  real_T c1_d17;
  real_T c1_d18;
  real_T c1_d19;
  real_T c1_d20;
  real_T c1_k_x[16];
  int32_T c1_k;
  int32_T c1_b_k;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 18U, 18U, c1_c_debug_family_names,
    c1_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t2, 0U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t3, 1U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t4, 2U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t5, 3U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t6, 4U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t12, 5U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t7, 6U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t8, 7U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t9, 8U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t10, 9U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t11, 10U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t13, 11U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_t14, 12U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargin, 13U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargout, 14U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_theta_1, 15U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_theta_2, 16U, c1_d_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_out1, 17U, c1_c_sf_marshallOut,
    c1_c_sf_marshallIn);
  CV_SCRIPT_FCN(1, 0);
  _SFD_SCRIPT_CALL(1U, chartInstance->c1_sfEvent, 8);
  c1_x = c1_theta_2;
  c1_t2 = c1_x;
  c1_b_x = c1_t2;
  c1_t2 = c1_b_x;
  c1_t2 = muDoubleScalarSin(c1_t2);
  _SFD_SCRIPT_CALL(1U, chartInstance->c1_sfEvent, 9);
  c1_c_x = c1_theta_1;
  c1_t3 = c1_c_x;
  c1_d_x = c1_t3;
  c1_t3 = c1_d_x;
  c1_t3 = muDoubleScalarSin(c1_t3);
  _SFD_SCRIPT_CALL(1U, chartInstance->c1_sfEvent, 10);
  c1_e_x = c1_theta_1;
  c1_t4 = c1_e_x;
  c1_f_x = c1_t4;
  c1_t4 = c1_f_x;
  c1_t4 = muDoubleScalarCos(c1_t4);
  _SFD_SCRIPT_CALL(1U, chartInstance->c1_sfEvent, 11);
  c1_t5 = c1_t2 * c1_t4;
  _SFD_SCRIPT_CALL(1U, chartInstance->c1_sfEvent, 12);
  c1_g_x = c1_theta_2;
  c1_t6 = c1_g_x;
  c1_h_x = c1_t6;
  c1_t6 = c1_h_x;
  c1_t6 = muDoubleScalarCos(c1_t6);
  _SFD_SCRIPT_CALL(1U, chartInstance->c1_sfEvent, 13);
  c1_t12 = c1_t3 * c1_t6;
  _SFD_SCRIPT_CALL(1U, chartInstance->c1_sfEvent, 14);
  c1_t7 = c1_t5 - c1_t12;
  _SFD_SCRIPT_CALL(1U, chartInstance->c1_sfEvent, 15);
  c1_t8 = c1_rdivide(chartInstance, 1.0, c1_t7);
  _SFD_SCRIPT_CALL(1U, chartInstance->c1_sfEvent, 16);
  c1_t9 = c1_theta_1 - c1_theta_2;
  _SFD_SCRIPT_CALL(1U, chartInstance->c1_sfEvent, 17);
  c1_i_x = c1_t9;
  c1_t10 = c1_i_x;
  c1_j_x = c1_t10;
  c1_t10 = c1_j_x;
  c1_t10 = muDoubleScalarSin(c1_t10);
  _SFD_SCRIPT_CALL(1U, chartInstance->c1_sfEvent, 18);
  c1_t11 = c1_t2 * 25.0;
  _SFD_SCRIPT_CALL(1U, chartInstance->c1_sfEvent, 19);
  c1_t13 = c1_t6 * 25.0;
  _SFD_SCRIPT_CALL(1U, chartInstance->c1_sfEvent, 20);
  c1_t14 = c1_t4 * c1_t8 * c1_rdivide(chartInstance, 1.0, 2.0);
  _SFD_SCRIPT_CALL(1U, chartInstance->c1_sfEvent, 21);
  c1_d9 = c1_rdivide(chartInstance, 1.0, 50.0);
  c1_d10 = c1_rdivide(chartInstance, 1.0, 50.0);
  c1_d11 = c1_rdivide(chartInstance, -1.0, 2.0);
  c1_d12 = c1_rdivide(chartInstance, -1.0, 2.0);
  c1_d13 = c1_rdivide(chartInstance, -1.0, 50.0);
  c1_d14 = c1_rdivide(chartInstance, -1.0, 50.0);
  c1_d15 = c1_rdivide(chartInstance, -100.0, 17.0);
  c1_d16 = c1_rdivide(chartInstance, 100.0, 17.0);
  c1_d17 = c1_rdivide(chartInstance, 100.0, 17.0);
  c1_d18 = c1_rdivide(chartInstance, -100.0, 17.0);
  c1_d19 = c1_rdivide(chartInstance, -100.0, 17.0);
  c1_d20 = c1_rdivide(chartInstance, 100.0, 17.0);
  c1_k_x[0] = c1_t8 * (c1_t11 - c1_t3 * c1_t10 * 37.0) * c1_d9;
  c1_k_x[1] = c1_t8 * (c1_t11 + c1_t3 * c1_t10 * 37.0) * c1_d10;
  c1_k_x[2] = c1_t3 * c1_t8 * c1_d11;
  c1_k_x[3] = c1_t3 * c1_t8 * c1_d12;
  c1_k_x[4] = c1_t8 * (c1_t13 - c1_t4 * c1_t10 * 37.0) * c1_d13;
  c1_k_x[5] = c1_t8 * (c1_t13 + c1_t4 * c1_t10 * 37.0) * c1_d14;
  c1_k_x[6] = c1_t14;
  c1_k_x[7] = c1_t14;
  c1_k_x[8] = c1_d15;
  c1_k_x[9] = c1_d16;
  c1_k_x[10] = 0.0;
  c1_k_x[11] = 0.0;
  c1_k_x[12] = c1_d17;
  c1_k_x[13] = c1_d18;
  c1_k_x[14] = c1_d19;
  c1_k_x[15] = c1_d20;
  for (c1_k = 1; c1_k < 17; c1_k++) {
    c1_b_k = c1_k;
    c1_out1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c1_b_k), 1, 16, 1, 0) - 1] = c1_k_x[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 16, 1, 0) - 1];
  }

  _SFD_SCRIPT_CALL(1U, chartInstance->c1_sfEvent, -21);
  _SFD_SYMBOL_SCOPE_POP();
}

static void init_script_number_translation(uint32_T c1_machineNumber, uint32_T
  c1_chartNumber, uint32_T c1_instanceNumber)
{
  (void)c1_machineNumber;
  _SFD_SCRIPT_TRANSLATION(c1_chartNumber, c1_instanceNumber, 0U,
    sf_debug_get_script_id(
    "C:\\Users\\nlashkar\\Google Drive\\F15\\Multi_robots_control\\4. Computed Torque Control\\Ali\\Tau_computed.m"));
  _SFD_SCRIPT_TRANSLATION(c1_chartNumber, c1_instanceNumber, 1U,
    sf_debug_get_script_id(
    "C:\\Users\\nlashkar\\Google Drive\\F15\\Multi_robots_control\\4. Computed Torque Control\\Ali\\B_inverse.m"));
}

static const mxArray *c1_sf_marshallOut(void *chartInstanceVoid, void *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i31;
  real_T c1_b_inData[4];
  int32_T c1_i32;
  real_T c1_u[4];
  const mxArray *c1_y = NULL;
  SFc1_Ali_modelInstanceStruct *chartInstance;
  chartInstance = (SFc1_Ali_modelInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  for (c1_i31 = 0; c1_i31 < 4; c1_i31++) {
    c1_b_inData[c1_i31] = (*(real_T (*)[4])c1_inData)[c1_i31];
  }

  for (c1_i32 = 0; c1_i32 < 4; c1_i32++) {
    c1_u[c1_i32] = c1_b_inData[c1_i32];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 1, 4), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static void c1_emlrt_marshallIn(SFc1_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c1_b_Tau, const char_T *c1_identifier, real_T c1_y[4])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_Tau), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_Tau);
}

static void c1_b_emlrt_marshallIn(SFc1_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[4])
{
  real_T c1_dv1[4];
  int32_T c1_i33;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv1, 1, 0, 0U, 1, 0U, 1, 4);
  for (c1_i33 = 0; c1_i33 < 4; c1_i33++) {
    c1_y[c1_i33] = c1_dv1[c1_i33];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_b_Tau;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y[4];
  int32_T c1_i34;
  SFc1_Ali_modelInstanceStruct *chartInstance;
  chartInstance = (SFc1_Ali_modelInstanceStruct *)chartInstanceVoid;
  c1_b_Tau = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_Tau), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_Tau);
  for (c1_i34 = 0; c1_i34 < 4; c1_i34++) {
    (*(real_T (*)[4])c1_outData)[c1_i34] = c1_y[c1_i34];
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_b_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i35;
  real_T c1_b_inData[8];
  int32_T c1_i36;
  real_T c1_u[8];
  const mxArray *c1_y = NULL;
  SFc1_Ali_modelInstanceStruct *chartInstance;
  chartInstance = (SFc1_Ali_modelInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  for (c1_i35 = 0; c1_i35 < 8; c1_i35++) {
    c1_b_inData[c1_i35] = (*(real_T (*)[8])c1_inData)[c1_i35];
  }

  for (c1_i36 = 0; c1_i36 < 8; c1_i36++) {
    c1_u[c1_i36] = c1_b_inData[c1_i36];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 1, 8), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static const mxArray *c1_c_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i37;
  int32_T c1_i38;
  int32_T c1_i39;
  real_T c1_b_inData[16];
  int32_T c1_i40;
  int32_T c1_i41;
  int32_T c1_i42;
  real_T c1_u[16];
  const mxArray *c1_y = NULL;
  SFc1_Ali_modelInstanceStruct *chartInstance;
  chartInstance = (SFc1_Ali_modelInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_i37 = 0;
  for (c1_i38 = 0; c1_i38 < 4; c1_i38++) {
    for (c1_i39 = 0; c1_i39 < 4; c1_i39++) {
      c1_b_inData[c1_i39 + c1_i37] = (*(real_T (*)[16])c1_inData)[c1_i39 +
        c1_i37];
    }

    c1_i37 += 4;
  }

  c1_i40 = 0;
  for (c1_i41 = 0; c1_i41 < 4; c1_i41++) {
    for (c1_i42 = 0; c1_i42 < 4; c1_i42++) {
      c1_u[c1_i42 + c1_i40] = c1_b_inData[c1_i42 + c1_i40];
    }

    c1_i40 += 4;
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 2, 4, 4), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static const mxArray *c1_d_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  real_T c1_u;
  const mxArray *c1_y = NULL;
  SFc1_Ali_modelInstanceStruct *chartInstance;
  chartInstance = (SFc1_Ali_modelInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_u = *(real_T *)c1_inData;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", &c1_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static real_T c1_c_emlrt_marshallIn(SFc1_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  real_T c1_y;
  real_T c1_d21;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_d21, 1, 0, 0U, 0, 0U, 0);
  c1_y = c1_d21;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_nargout;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y;
  SFc1_Ali_modelInstanceStruct *chartInstance;
  chartInstance = (SFc1_Ali_modelInstanceStruct *)chartInstanceVoid;
  c1_nargout = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_nargout), &c1_thisId);
  sf_mex_destroy(&c1_nargout);
  *(real_T *)c1_outData = c1_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

static void c1_d_emlrt_marshallIn(SFc1_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[16])
{
  real_T c1_dv2[16];
  int32_T c1_i43;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv2, 1, 0, 0U, 1, 0U, 2, 4, 4);
  for (c1_i43 = 0; c1_i43 < 16; c1_i43++) {
    c1_y[c1_i43] = c1_dv2[c1_i43];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_out1;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y[16];
  int32_T c1_i44;
  int32_T c1_i45;
  int32_T c1_i46;
  SFc1_Ali_modelInstanceStruct *chartInstance;
  chartInstance = (SFc1_Ali_modelInstanceStruct *)chartInstanceVoid;
  c1_out1 = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_out1), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_out1);
  c1_i44 = 0;
  for (c1_i45 = 0; c1_i45 < 4; c1_i45++) {
    for (c1_i46 = 0; c1_i46 < 4; c1_i46++) {
      (*(real_T (*)[16])c1_outData)[c1_i46 + c1_i44] = c1_y[c1_i46 + c1_i44];
    }

    c1_i44 += 4;
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

const mxArray *sf_c1_Ali_model_get_eml_resolved_functions_info(void)
{
  const mxArray *c1_nameCaptureInfo = NULL;
  c1_nameCaptureInfo = NULL;
  sf_mex_assign(&c1_nameCaptureInfo, sf_mex_createstruct("structure", 2, 60, 1),
                false);
  c1_info_helper(&c1_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c1_nameCaptureInfo);
  return c1_nameCaptureInfo;
}

static void c1_info_helper(const mxArray **c1_info)
{
  const mxArray *c1_rhs0 = NULL;
  const mxArray *c1_lhs0 = NULL;
  const mxArray *c1_rhs1 = NULL;
  const mxArray *c1_lhs1 = NULL;
  const mxArray *c1_rhs2 = NULL;
  const mxArray *c1_lhs2 = NULL;
  const mxArray *c1_rhs3 = NULL;
  const mxArray *c1_lhs3 = NULL;
  const mxArray *c1_rhs4 = NULL;
  const mxArray *c1_lhs4 = NULL;
  const mxArray *c1_rhs5 = NULL;
  const mxArray *c1_lhs5 = NULL;
  const mxArray *c1_rhs6 = NULL;
  const mxArray *c1_lhs6 = NULL;
  const mxArray *c1_rhs7 = NULL;
  const mxArray *c1_lhs7 = NULL;
  const mxArray *c1_rhs8 = NULL;
  const mxArray *c1_lhs8 = NULL;
  const mxArray *c1_rhs9 = NULL;
  const mxArray *c1_lhs9 = NULL;
  const mxArray *c1_rhs10 = NULL;
  const mxArray *c1_lhs10 = NULL;
  const mxArray *c1_rhs11 = NULL;
  const mxArray *c1_lhs11 = NULL;
  const mxArray *c1_rhs12 = NULL;
  const mxArray *c1_lhs12 = NULL;
  const mxArray *c1_rhs13 = NULL;
  const mxArray *c1_lhs13 = NULL;
  const mxArray *c1_rhs14 = NULL;
  const mxArray *c1_lhs14 = NULL;
  const mxArray *c1_rhs15 = NULL;
  const mxArray *c1_lhs15 = NULL;
  const mxArray *c1_rhs16 = NULL;
  const mxArray *c1_lhs16 = NULL;
  const mxArray *c1_rhs17 = NULL;
  const mxArray *c1_lhs17 = NULL;
  const mxArray *c1_rhs18 = NULL;
  const mxArray *c1_lhs18 = NULL;
  const mxArray *c1_rhs19 = NULL;
  const mxArray *c1_lhs19 = NULL;
  const mxArray *c1_rhs20 = NULL;
  const mxArray *c1_lhs20 = NULL;
  const mxArray *c1_rhs21 = NULL;
  const mxArray *c1_lhs21 = NULL;
  const mxArray *c1_rhs22 = NULL;
  const mxArray *c1_lhs22 = NULL;
  const mxArray *c1_rhs23 = NULL;
  const mxArray *c1_lhs23 = NULL;
  const mxArray *c1_rhs24 = NULL;
  const mxArray *c1_lhs24 = NULL;
  const mxArray *c1_rhs25 = NULL;
  const mxArray *c1_lhs25 = NULL;
  const mxArray *c1_rhs26 = NULL;
  const mxArray *c1_lhs26 = NULL;
  const mxArray *c1_rhs27 = NULL;
  const mxArray *c1_lhs27 = NULL;
  const mxArray *c1_rhs28 = NULL;
  const mxArray *c1_lhs28 = NULL;
  const mxArray *c1_rhs29 = NULL;
  const mxArray *c1_lhs29 = NULL;
  const mxArray *c1_rhs30 = NULL;
  const mxArray *c1_lhs30 = NULL;
  const mxArray *c1_rhs31 = NULL;
  const mxArray *c1_lhs31 = NULL;
  const mxArray *c1_rhs32 = NULL;
  const mxArray *c1_lhs32 = NULL;
  const mxArray *c1_rhs33 = NULL;
  const mxArray *c1_lhs33 = NULL;
  const mxArray *c1_rhs34 = NULL;
  const mxArray *c1_lhs34 = NULL;
  const mxArray *c1_rhs35 = NULL;
  const mxArray *c1_lhs35 = NULL;
  const mxArray *c1_rhs36 = NULL;
  const mxArray *c1_lhs36 = NULL;
  const mxArray *c1_rhs37 = NULL;
  const mxArray *c1_lhs37 = NULL;
  const mxArray *c1_rhs38 = NULL;
  const mxArray *c1_lhs38 = NULL;
  const mxArray *c1_rhs39 = NULL;
  const mxArray *c1_lhs39 = NULL;
  const mxArray *c1_rhs40 = NULL;
  const mxArray *c1_lhs40 = NULL;
  const mxArray *c1_rhs41 = NULL;
  const mxArray *c1_lhs41 = NULL;
  const mxArray *c1_rhs42 = NULL;
  const mxArray *c1_lhs42 = NULL;
  const mxArray *c1_rhs43 = NULL;
  const mxArray *c1_lhs43 = NULL;
  const mxArray *c1_rhs44 = NULL;
  const mxArray *c1_lhs44 = NULL;
  const mxArray *c1_rhs45 = NULL;
  const mxArray *c1_lhs45 = NULL;
  const mxArray *c1_rhs46 = NULL;
  const mxArray *c1_lhs46 = NULL;
  const mxArray *c1_rhs47 = NULL;
  const mxArray *c1_lhs47 = NULL;
  const mxArray *c1_rhs48 = NULL;
  const mxArray *c1_lhs48 = NULL;
  const mxArray *c1_rhs49 = NULL;
  const mxArray *c1_lhs49 = NULL;
  const mxArray *c1_rhs50 = NULL;
  const mxArray *c1_lhs50 = NULL;
  const mxArray *c1_rhs51 = NULL;
  const mxArray *c1_lhs51 = NULL;
  const mxArray *c1_rhs52 = NULL;
  const mxArray *c1_lhs52 = NULL;
  const mxArray *c1_rhs53 = NULL;
  const mxArray *c1_lhs53 = NULL;
  const mxArray *c1_rhs54 = NULL;
  const mxArray *c1_lhs54 = NULL;
  const mxArray *c1_rhs55 = NULL;
  const mxArray *c1_lhs55 = NULL;
  const mxArray *c1_rhs56 = NULL;
  const mxArray *c1_lhs56 = NULL;
  const mxArray *c1_rhs57 = NULL;
  const mxArray *c1_lhs57 = NULL;
  const mxArray *c1_rhs58 = NULL;
  const mxArray *c1_lhs58 = NULL;
  const mxArray *c1_rhs59 = NULL;
  const mxArray *c1_lhs59 = NULL;
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("Tau_computed"), "name", "name",
                  0);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]C:/Users/nlashkar/Google Drive/F15/Multi_robots_control/4. Computed Torque Control/Ali/Tau_computed.m"),
                  "resolved", "resolved", 0);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1445640297U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c1_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]C:/Users/nlashkar/Google Drive/F15/Multi_robots_control/4. Computed Torque Control/Ali/Tau_computed.m"),
                  "context", "context", 1);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("sin"), "name", "name", 1);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 1);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395346504U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c1_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 2);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 2);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 2);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286840336U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c1_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]C:/Users/nlashkar/Google Drive/F15/Multi_robots_control/4. Computed Torque Control/Ali/Tau_computed.m"),
                  "context", "context", 3);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("cos"), "name", "name", 3);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 3);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395346496U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c1_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286840322U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c1_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]C:/Users/nlashkar/Google Drive/F15/Multi_robots_control/4. Computed Torque Control/Ali/Tau_computed.m"),
                  "context", "context", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("power"), "name", "name", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "resolved",
                  "resolved", 5);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395346506U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c1_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395949856U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c1_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1376002288U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c1_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410829370U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c1_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1376002288U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c1_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "context", "context", 10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalexpAlloc"),
                  "name", "name", 10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410829370U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c1_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("floor"), "name", "name", 11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "resolved",
                  "resolved", 11);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363731854U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c1_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395949856U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c1_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286840326U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c1_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 14);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1376002288U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c1_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]C:/Users/nlashkar/Google Drive/F15/Multi_robots_control/4. Computed Torque Control/Ali/Tau_computed.m"),
                  "context", "context", 15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("rdivide"), "name", "name", 15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 15);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363731880U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c1_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 16);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 16);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395949856U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c1_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 17);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 17);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 17);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286840396U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c1_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 18);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_div"), "name", "name", 18);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 18);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1386445552U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c1_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 19);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 19);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 19);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410829370U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c1_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 20);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("B_inverse"), "name", "name",
                  20);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]C:/Users/nlashkar/Google Drive/F15/Multi_robots_control/4. Computed Torque Control/Ali/B_inverse.m"),
                  "resolved", "resolved", 20);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1445640297U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c1_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]C:/Users/nlashkar/Google Drive/F15/Multi_robots_control/4. Computed Torque Control/Ali/B_inverse.m"),
                  "context", "context", 21);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("sin"), "name", "name", 21);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 21);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395346504U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c1_rhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]C:/Users/nlashkar/Google Drive/F15/Multi_robots_control/4. Computed Torque Control/Ali/B_inverse.m"),
                  "context", "context", 22);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("cos"), "name", "name", 22);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 22);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395346496U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c1_rhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]C:/Users/nlashkar/Google Drive/F15/Multi_robots_control/4. Computed Torque Control/Ali/B_inverse.m"),
                  "context", "context", 23);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("rdivide"), "name", "name", 23);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 23);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363731880U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c1_rhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]C:/Users/nlashkar/Google Drive/F15/Multi_robots_control/4. Computed Torque Control/Ali/B_inverse.m"),
                  "context", "context", 24);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("reshape"), "name", "name", 24);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/reshape.m"), "resolved",
                  "resolved", 24);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1378317582U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c1_rhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/reshape.m"), "context",
                  "context", 25);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 25);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 25);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323192178U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c1_rhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/reshape.m!reshape_varargin_to_size"),
                  "context", "context", 26);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 26);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 26);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323192178U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c1_rhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/reshape.m!reshape_varargin_to_size"),
                  "context", "context", 27);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_assert_valid_size_arg"),
                  "name", "name", 27);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "resolved", "resolved", 27);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1368204630U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c1_rhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 28);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 28);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 28);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395949856U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c1_rhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isintegral"),
                  "context", "context", 29);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isinf"), "name", "name", 29);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "resolved",
                  "resolved", 29);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363731856U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c1_rhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "context",
                  "context", 30);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 30);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 30);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395949856U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c1_rhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 31);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_is_integer_class"), "name",
                  "name", 31);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_integer_class.m"),
                  "resolved", "resolved", 31);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286840382U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c1_rhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 32);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmax"), "name", "name", 32);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 32);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362283482U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c1_rhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "context",
                  "context", 33);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 33);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 33);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1393352458U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c1_rhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 34);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmin"), "name", "name", 34);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 34);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 34);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362283482U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c1_rhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs34), "lhs", "lhs",
                  34);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "context",
                  "context", 35);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 35);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 35);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 35);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1393352458U), "fileTimeLo",
                  "fileTimeLo", 35);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 35);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 35);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 35);
  sf_mex_assign(&c1_rhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs35), "rhs", "rhs",
                  35);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs35), "lhs", "lhs",
                  35);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 36);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexIntRelop"),
                  "name", "name", 36);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 36);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"),
                  "resolved", "resolved", 36);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326749922U), "fileTimeLo",
                  "fileTimeLo", 36);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 36);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 36);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 36);
  sf_mex_assign(&c1_rhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs36), "rhs", "rhs",
                  36);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs36), "lhs", "lhs",
                  36);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!apply_float_relop"),
                  "context", "context", 37);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 37);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 37);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 37);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1393352458U), "fileTimeLo",
                  "fileTimeLo", 37);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 37);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 37);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 37);
  sf_mex_assign(&c1_rhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs37), "rhs", "rhs",
                  37);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs37), "lhs", "lhs",
                  37);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!float_class_contains_indexIntClass"),
                  "context", "context", 38);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 38);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 38);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 38);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326749596U), "fileTimeLo",
                  "fileTimeLo", 38);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 38);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 38);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 38);
  sf_mex_assign(&c1_rhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs38), "rhs", "rhs",
                  38);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs38), "lhs", "lhs",
                  38);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!is_signed_indexIntClass"),
                  "context", "context", 39);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmin"), "name", "name", 39);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 39);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 39);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362283482U), "fileTimeLo",
                  "fileTimeLo", 39);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 39);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 39);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 39);
  sf_mex_assign(&c1_rhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs39), "rhs", "rhs",
                  39);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs39), "lhs", "lhs",
                  39);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 40);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 40);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 40);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 40);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323192178U), "fileTimeLo",
                  "fileTimeLo", 40);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 40);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 40);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 40);
  sf_mex_assign(&c1_rhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs40), "rhs", "rhs",
                  40);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs40), "lhs", "lhs",
                  40);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 41);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmax"), "name", "name", 41);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 41);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 41);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362283482U), "fileTimeLo",
                  "fileTimeLo", 41);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 41);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 41);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 41);
  sf_mex_assign(&c1_rhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs41), "rhs", "rhs",
                  41);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs41), "lhs", "lhs",
                  41);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/reshape.m"), "context",
                  "context", 42);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 42);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 42);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 42);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1376002288U), "fileTimeLo",
                  "fileTimeLo", 42);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 42);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 42);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 42);
  sf_mex_assign(&c1_rhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs42), "rhs", "rhs",
                  42);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs42), "lhs", "lhs",
                  42);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/reshape.m"), "context",
                  "context", 43);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 43);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 43);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 43);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1397279022U), "fileTimeLo",
                  "fileTimeLo", 43);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 43);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 43);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 43);
  sf_mex_assign(&c1_rhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs43), "rhs", "rhs",
                  43);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs43), "lhs", "lhs",
                  43);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 44);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isfi"), "name", "name", 44);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 44);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/fixedpoint/isfi.m"), "resolved",
                  "resolved", 44);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1346531958U), "fileTimeLo",
                  "fileTimeLo", 44);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 44);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 44);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 44);
  sf_mex_assign(&c1_rhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs44), "rhs", "rhs",
                  44);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs44), "lhs", "lhs",
                  44);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/fixedpoint/isfi.m"), "context",
                  "context", 45);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isnumerictype"), "name",
                  "name", 45);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 45);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/fixedpoint/isnumerictype.m"), "resolved",
                  "resolved", 45);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1398897198U), "fileTimeLo",
                  "fileTimeLo", 45);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 45);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 45);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 45);
  sf_mex_assign(&c1_rhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs45), "rhs", "rhs",
                  45);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs45), "lhs", "lhs",
                  45);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 46);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmax"), "name", "name", 46);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 46);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 46);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362283482U), "fileTimeLo",
                  "fileTimeLo", 46);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 46);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 46);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 46);
  sf_mex_assign(&c1_rhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs46), "rhs", "rhs",
                  46);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs46), "lhs", "lhs",
                  46);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 47);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmin"), "name", "name", 47);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 47);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 47);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362283482U), "fileTimeLo",
                  "fileTimeLo", 47);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 47);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 47);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 47);
  sf_mex_assign(&c1_rhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs47), "rhs", "rhs",
                  47);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs47), "lhs", "lhs",
                  47);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 48);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 48);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 48);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 48);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1383898894U), "fileTimeLo",
                  "fileTimeLo", 48);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 48);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 48);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 48);
  sf_mex_assign(&c1_rhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs48), "rhs", "rhs",
                  48);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs48), "lhs", "lhs",
                  48);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 49);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 49);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 49);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 49);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395949856U), "fileTimeLo",
                  "fileTimeLo", 49);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 49);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 49);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 49);
  sf_mex_assign(&c1_rhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs49), "rhs", "rhs",
                  49);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs49), "lhs", "lhs",
                  49);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 50);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 50);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 50);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 50);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323192178U), "fileTimeLo",
                  "fileTimeLo", 50);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 50);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 50);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 50);
  sf_mex_assign(&c1_rhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs50), "rhs", "rhs",
                  50);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs50), "lhs", "lhs",
                  50);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 51);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 51);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 51);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 51);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1376002288U), "fileTimeLo",
                  "fileTimeLo", 51);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 51);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 51);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 51);
  sf_mex_assign(&c1_rhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs51), "rhs", "rhs",
                  51);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs51), "lhs", "lhs",
                  51);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 52);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  52);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 52);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 52);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1376002290U), "fileTimeLo",
                  "fileTimeLo", 52);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 52);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 52);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 52);
  sf_mex_assign(&c1_rhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs52), "rhs", "rhs",
                  52);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs52), "lhs", "lhs",
                  52);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 53);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 53);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 53);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 53);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410829372U), "fileTimeLo",
                  "fileTimeLo", 53);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 53);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 53);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 53);
  sf_mex_assign(&c1_rhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs53), "rhs", "rhs",
                  53);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs53), "lhs", "lhs",
                  53);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 54);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.xgemm"),
                  "name", "name", 54);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 54);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "resolved", "resolved", 54);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410829370U), "fileTimeLo",
                  "fileTimeLo", 54);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 54);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 54);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 54);
  sf_mex_assign(&c1_rhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs54), "rhs", "rhs",
                  54);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs54), "lhs", "lhs",
                  54);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 55);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 55);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 55);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 55);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410829370U), "fileTimeLo",
                  "fileTimeLo", 55);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 55);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 55);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 55);
  sf_mex_assign(&c1_rhs55, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs55, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs55), "rhs", "rhs",
                  55);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs55), "lhs", "lhs",
                  55);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 56);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 56);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 56);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 56);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410829372U), "fileTimeLo",
                  "fileTimeLo", 56);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 56);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 56);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 56);
  sf_mex_assign(&c1_rhs56, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs56, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs56), "rhs", "rhs",
                  56);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs56), "lhs", "lhs",
                  56);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 57);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 57);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 57);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 57);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1393352458U), "fileTimeLo",
                  "fileTimeLo", 57);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 57);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 57);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 57);
  sf_mex_assign(&c1_rhs57, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs57, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs57), "rhs", "rhs",
                  57);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs57), "lhs", "lhs",
                  57);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 58);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 58);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 58);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 58);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410829370U), "fileTimeLo",
                  "fileTimeLo", 58);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 58);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 58);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 58);
  sf_mex_assign(&c1_rhs58, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs58, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs58), "rhs", "rhs",
                  58);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs58), "lhs", "lhs",
                  58);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 59);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.refblas.xgemm"),
                  "name", "name", 59);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 59);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "resolved", "resolved", 59);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410829372U), "fileTimeLo",
                  "fileTimeLo", 59);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 59);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 59);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 59);
  sf_mex_assign(&c1_rhs59, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs59, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs59), "rhs", "rhs",
                  59);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs59), "lhs", "lhs",
                  59);
  sf_mex_destroy(&c1_rhs0);
  sf_mex_destroy(&c1_lhs0);
  sf_mex_destroy(&c1_rhs1);
  sf_mex_destroy(&c1_lhs1);
  sf_mex_destroy(&c1_rhs2);
  sf_mex_destroy(&c1_lhs2);
  sf_mex_destroy(&c1_rhs3);
  sf_mex_destroy(&c1_lhs3);
  sf_mex_destroy(&c1_rhs4);
  sf_mex_destroy(&c1_lhs4);
  sf_mex_destroy(&c1_rhs5);
  sf_mex_destroy(&c1_lhs5);
  sf_mex_destroy(&c1_rhs6);
  sf_mex_destroy(&c1_lhs6);
  sf_mex_destroy(&c1_rhs7);
  sf_mex_destroy(&c1_lhs7);
  sf_mex_destroy(&c1_rhs8);
  sf_mex_destroy(&c1_lhs8);
  sf_mex_destroy(&c1_rhs9);
  sf_mex_destroy(&c1_lhs9);
  sf_mex_destroy(&c1_rhs10);
  sf_mex_destroy(&c1_lhs10);
  sf_mex_destroy(&c1_rhs11);
  sf_mex_destroy(&c1_lhs11);
  sf_mex_destroy(&c1_rhs12);
  sf_mex_destroy(&c1_lhs12);
  sf_mex_destroy(&c1_rhs13);
  sf_mex_destroy(&c1_lhs13);
  sf_mex_destroy(&c1_rhs14);
  sf_mex_destroy(&c1_lhs14);
  sf_mex_destroy(&c1_rhs15);
  sf_mex_destroy(&c1_lhs15);
  sf_mex_destroy(&c1_rhs16);
  sf_mex_destroy(&c1_lhs16);
  sf_mex_destroy(&c1_rhs17);
  sf_mex_destroy(&c1_lhs17);
  sf_mex_destroy(&c1_rhs18);
  sf_mex_destroy(&c1_lhs18);
  sf_mex_destroy(&c1_rhs19);
  sf_mex_destroy(&c1_lhs19);
  sf_mex_destroy(&c1_rhs20);
  sf_mex_destroy(&c1_lhs20);
  sf_mex_destroy(&c1_rhs21);
  sf_mex_destroy(&c1_lhs21);
  sf_mex_destroy(&c1_rhs22);
  sf_mex_destroy(&c1_lhs22);
  sf_mex_destroy(&c1_rhs23);
  sf_mex_destroy(&c1_lhs23);
  sf_mex_destroy(&c1_rhs24);
  sf_mex_destroy(&c1_lhs24);
  sf_mex_destroy(&c1_rhs25);
  sf_mex_destroy(&c1_lhs25);
  sf_mex_destroy(&c1_rhs26);
  sf_mex_destroy(&c1_lhs26);
  sf_mex_destroy(&c1_rhs27);
  sf_mex_destroy(&c1_lhs27);
  sf_mex_destroy(&c1_rhs28);
  sf_mex_destroy(&c1_lhs28);
  sf_mex_destroy(&c1_rhs29);
  sf_mex_destroy(&c1_lhs29);
  sf_mex_destroy(&c1_rhs30);
  sf_mex_destroy(&c1_lhs30);
  sf_mex_destroy(&c1_rhs31);
  sf_mex_destroy(&c1_lhs31);
  sf_mex_destroy(&c1_rhs32);
  sf_mex_destroy(&c1_lhs32);
  sf_mex_destroy(&c1_rhs33);
  sf_mex_destroy(&c1_lhs33);
  sf_mex_destroy(&c1_rhs34);
  sf_mex_destroy(&c1_lhs34);
  sf_mex_destroy(&c1_rhs35);
  sf_mex_destroy(&c1_lhs35);
  sf_mex_destroy(&c1_rhs36);
  sf_mex_destroy(&c1_lhs36);
  sf_mex_destroy(&c1_rhs37);
  sf_mex_destroy(&c1_lhs37);
  sf_mex_destroy(&c1_rhs38);
  sf_mex_destroy(&c1_lhs38);
  sf_mex_destroy(&c1_rhs39);
  sf_mex_destroy(&c1_lhs39);
  sf_mex_destroy(&c1_rhs40);
  sf_mex_destroy(&c1_lhs40);
  sf_mex_destroy(&c1_rhs41);
  sf_mex_destroy(&c1_lhs41);
  sf_mex_destroy(&c1_rhs42);
  sf_mex_destroy(&c1_lhs42);
  sf_mex_destroy(&c1_rhs43);
  sf_mex_destroy(&c1_lhs43);
  sf_mex_destroy(&c1_rhs44);
  sf_mex_destroy(&c1_lhs44);
  sf_mex_destroy(&c1_rhs45);
  sf_mex_destroy(&c1_lhs45);
  sf_mex_destroy(&c1_rhs46);
  sf_mex_destroy(&c1_lhs46);
  sf_mex_destroy(&c1_rhs47);
  sf_mex_destroy(&c1_lhs47);
  sf_mex_destroy(&c1_rhs48);
  sf_mex_destroy(&c1_lhs48);
  sf_mex_destroy(&c1_rhs49);
  sf_mex_destroy(&c1_lhs49);
  sf_mex_destroy(&c1_rhs50);
  sf_mex_destroy(&c1_lhs50);
  sf_mex_destroy(&c1_rhs51);
  sf_mex_destroy(&c1_lhs51);
  sf_mex_destroy(&c1_rhs52);
  sf_mex_destroy(&c1_lhs52);
  sf_mex_destroy(&c1_rhs53);
  sf_mex_destroy(&c1_lhs53);
  sf_mex_destroy(&c1_rhs54);
  sf_mex_destroy(&c1_lhs54);
  sf_mex_destroy(&c1_rhs55);
  sf_mex_destroy(&c1_lhs55);
  sf_mex_destroy(&c1_rhs56);
  sf_mex_destroy(&c1_lhs56);
  sf_mex_destroy(&c1_rhs57);
  sf_mex_destroy(&c1_lhs57);
  sf_mex_destroy(&c1_rhs58);
  sf_mex_destroy(&c1_lhs58);
  sf_mex_destroy(&c1_rhs59);
  sf_mex_destroy(&c1_lhs59);
}

static const mxArray *c1_emlrt_marshallOut(const char * c1_u)
{
  const mxArray *c1_y = NULL;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c1_u)), false);
  return c1_y;
}

static const mxArray *c1_b_emlrt_marshallOut(const uint32_T c1_u)
{
  const mxArray *c1_y = NULL;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", &c1_u, 7, 0U, 0U, 0U, 0), false);
  return c1_y;
}

static real_T c1_power(SFc1_Ali_modelInstanceStruct *chartInstance, real_T c1_a)
{
  real_T c1_b_a;
  real_T c1_ak;
  real_T c1_c_a;
  c1_b_a = c1_a;
  c1_eml_scalar_eg(chartInstance);
  c1_ak = c1_b_a;
  c1_c_a = c1_ak;
  c1_eml_scalar_eg(chartInstance);
  return c1_c_a * c1_c_a;
}

static void c1_eml_scalar_eg(SFc1_Ali_modelInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static real_T c1_rdivide(SFc1_Ali_modelInstanceStruct *chartInstance, real_T
  c1_x, real_T c1_y)
{
  real_T c1_b_x;
  real_T c1_b_y;
  real_T c1_c_x;
  real_T c1_c_y;
  (void)chartInstance;
  c1_b_x = c1_x;
  c1_b_y = c1_y;
  c1_c_x = c1_b_x;
  c1_c_y = c1_b_y;
  return c1_c_x / c1_c_y;
}

static void c1_b_eml_scalar_eg(SFc1_Ali_modelInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c1_eml_xgemm(SFc1_Ali_modelInstanceStruct *chartInstance, real_T
  c1_A[16], real_T c1_B[4], real_T c1_C[4], real_T c1_b_C[4])
{
  int32_T c1_i47;
  int32_T c1_i48;
  real_T c1_b_A[16];
  int32_T c1_i49;
  real_T c1_b_B[4];
  for (c1_i47 = 0; c1_i47 < 4; c1_i47++) {
    c1_b_C[c1_i47] = c1_C[c1_i47];
  }

  for (c1_i48 = 0; c1_i48 < 16; c1_i48++) {
    c1_b_A[c1_i48] = c1_A[c1_i48];
  }

  for (c1_i49 = 0; c1_i49 < 4; c1_i49++) {
    c1_b_B[c1_i49] = c1_B[c1_i49];
  }

  c1_b_eml_xgemm(chartInstance, c1_b_A, c1_b_B, c1_b_C);
}

static const mxArray *c1_e_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_u;
  const mxArray *c1_y = NULL;
  SFc1_Ali_modelInstanceStruct *chartInstance;
  chartInstance = (SFc1_Ali_modelInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_u = *(int32_T *)c1_inData;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", &c1_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static int32_T c1_e_emlrt_marshallIn(SFc1_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  int32_T c1_y;
  int32_T c1_i50;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_i50, 1, 6, 0U, 0, 0U, 0);
  c1_y = c1_i50;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_b_sfEvent;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  int32_T c1_y;
  SFc1_Ali_modelInstanceStruct *chartInstance;
  chartInstance = (SFc1_Ali_modelInstanceStruct *)chartInstanceVoid;
  c1_b_sfEvent = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_sfEvent),
    &c1_thisId);
  sf_mex_destroy(&c1_b_sfEvent);
  *(int32_T *)c1_outData = c1_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

static uint8_T c1_f_emlrt_marshallIn(SFc1_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c1_b_is_active_c1_Ali_model, const char_T *c1_identifier)
{
  uint8_T c1_y;
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_g_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c1_b_is_active_c1_Ali_model), &c1_thisId);
  sf_mex_destroy(&c1_b_is_active_c1_Ali_model);
  return c1_y;
}

static uint8_T c1_g_emlrt_marshallIn(SFc1_Ali_modelInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  uint8_T c1_y;
  uint8_T c1_u0;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_u0, 1, 3, 0U, 0, 0U, 0);
  c1_y = c1_u0;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_b_eml_xgemm(SFc1_Ali_modelInstanceStruct *chartInstance, real_T
  c1_A[16], real_T c1_B[4], real_T c1_C[4])
{
  int32_T c1_i51;
  int32_T c1_i52;
  int32_T c1_i53;
  (void)chartInstance;
  for (c1_i51 = 0; c1_i51 < 4; c1_i51++) {
    c1_C[c1_i51] = 0.0;
    c1_i52 = 0;
    for (c1_i53 = 0; c1_i53 < 4; c1_i53++) {
      c1_C[c1_i51] += c1_A[c1_i52 + c1_i51] * c1_B[c1_i53];
      c1_i52 += 4;
    }
  }
}

static void init_dsm_address_info(SFc1_Ali_modelInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void init_simulink_io_address(SFc1_Ali_modelInstanceStruct *chartInstance)
{
  chartInstance->c1_kp = (real_T (*)[16])ssGetInputPortSignal_wrapper
    (chartInstance->S, 0);
  chartInstance->c1_kd = (real_T (*)[16])ssGetInputPortSignal_wrapper
    (chartInstance->S, 1);
  chartInstance->c1_e = (real_T (*)[4])ssGetInputPortSignal_wrapper
    (chartInstance->S, 2);
  chartInstance->c1_e_dot = (real_T (*)[4])ssGetInputPortSignal_wrapper
    (chartInstance->S, 3);
  chartInstance->c1_q = (real_T (*)[8])ssGetInputPortSignal_wrapper
    (chartInstance->S, 4);
  chartInstance->c1_Tau = (real_T (*)[4])ssGetOutputPortSignal_wrapper
    (chartInstance->S, 1);
  chartInstance->c1_q_dot = (real_T (*)[8])ssGetInputPortSignal_wrapper
    (chartInstance->S, 5);
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

void sf_c1_Ali_model_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2984294100U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2687222538U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(2091818659U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(208109663U);
}

mxArray* sf_c1_Ali_model_get_post_codegen_info(void);
mxArray *sf_c1_Ali_model_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals", "postCodegenInfo" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1, 1, sizeof
    (autoinheritanceFields)/sizeof(autoinheritanceFields[0]),
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("m3aY8v09NIwyuzIIOgo35E");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,6,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(4);
      pr[1] = (double)(4);
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
      pr[1] = (double)(4);
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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(4);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(4);
      pr[1] = (double)(1);
      mxSetField(mxData,3,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,3,"type",mxType);
    }

    mxSetField(mxData,3,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(8);
      pr[1] = (double)(1);
      mxSetField(mxData,4,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,4,"type",mxType);
    }

    mxSetField(mxData,4,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(8);
      pr[1] = (double)(1);
      mxSetField(mxData,5,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,5,"type",mxType);
    }

    mxSetField(mxData,5,"complexity",mxCreateDoubleScalar(0));
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
      pr[0] = (double)(4);
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
    mxArray* mxPostCodegenInfo = sf_c1_Ali_model_get_post_codegen_info();
    mxSetField(mxAutoinheritanceInfo,0,"postCodegenInfo",mxPostCodegenInfo);
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c1_Ali_model_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c1_Ali_model_jit_fallback_info(void)
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

mxArray *sf_c1_Ali_model_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

mxArray* sf_c1_Ali_model_get_post_codegen_info(void)
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

static const mxArray *sf_get_sim_state_info_c1_Ali_model(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"Tau\",},{M[8],M[0],T\"is_active_c1_Ali_model\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c1_Ali_model_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc1_Ali_modelInstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc1_Ali_modelInstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _Ali_modelMachineNumber_,
           1,
           1,
           1,
           0,
           7,
           0,
           0,
           0,
           0,
           2,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"kp");
          _SFD_SET_DATA_PROPS(1,1,1,0,"kd");
          _SFD_SET_DATA_PROPS(2,1,1,0,"e");
          _SFD_SET_DATA_PROPS(3,1,1,0,"e_dot");
          _SFD_SET_DATA_PROPS(4,1,1,0,"q");
          _SFD_SET_DATA_PROPS(5,2,0,1,"Tau");
          _SFD_SET_DATA_PROPS(6,1,1,0,"q_dot");
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
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,590);
        _SFD_CV_INIT_SCRIPT(0,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_SCRIPT_FCN(0,0,"Tau_computed",0,-1,1859);
        _SFD_CV_INIT_SCRIPT(1,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_SCRIPT_FCN(1,0,"B_inverse",0,-1,744);

        {
          unsigned int dimVector[2];
          dimVector[0]= 4;
          dimVector[1]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 4;
          dimVector[1]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 8;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)
            c1_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 8;
          _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_VALUE_PTR(0U, *chartInstance->c1_kp);
        _SFD_SET_DATA_VALUE_PTR(1U, *chartInstance->c1_kd);
        _SFD_SET_DATA_VALUE_PTR(2U, *chartInstance->c1_e);
        _SFD_SET_DATA_VALUE_PTR(3U, *chartInstance->c1_e_dot);
        _SFD_SET_DATA_VALUE_PTR(4U, *chartInstance->c1_q);
        _SFD_SET_DATA_VALUE_PTR(5U, *chartInstance->c1_Tau);
        _SFD_SET_DATA_VALUE_PTR(6U, *chartInstance->c1_q_dot);
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
  return "ktfoFxOo8LIrZmK0msZsIB";
}

static void sf_opaque_initialize_c1_Ali_model(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc1_Ali_modelInstanceStruct*) chartInstanceVar
    )->S,0);
  initialize_params_c1_Ali_model((SFc1_Ali_modelInstanceStruct*)
    chartInstanceVar);
  initialize_c1_Ali_model((SFc1_Ali_modelInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c1_Ali_model(void *chartInstanceVar)
{
  enable_c1_Ali_model((SFc1_Ali_modelInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c1_Ali_model(void *chartInstanceVar)
{
  disable_c1_Ali_model((SFc1_Ali_modelInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c1_Ali_model(void *chartInstanceVar)
{
  sf_gateway_c1_Ali_model((SFc1_Ali_modelInstanceStruct*) chartInstanceVar);
}

static const mxArray* sf_opaque_get_sim_state_c1_Ali_model(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  return get_sim_state_c1_Ali_model((SFc1_Ali_modelInstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
}

static void sf_opaque_set_sim_state_c1_Ali_model(SimStruct* S, const mxArray *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  set_sim_state_c1_Ali_model((SFc1_Ali_modelInstanceStruct*)
    chartInfo->chartInstance, st);
}

static void sf_opaque_terminate_c1_Ali_model(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc1_Ali_modelInstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_Ali_model_optimization_info();
    }

    finalize_c1_Ali_model((SFc1_Ali_modelInstanceStruct*) chartInstanceVar);
    utFree(chartInstanceVar);
    if (crtInfo != NULL) {
      utFree(crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc1_Ali_model((SFc1_Ali_modelInstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c1_Ali_model(SimStruct *S)
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
    initialize_params_c1_Ali_model((SFc1_Ali_modelInstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c1_Ali_model(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_Ali_model_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,1);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,1,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,1,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,1);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 4, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 5, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,1,6);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,1,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 6; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,1);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(2590820365U));
  ssSetChecksum1(S,(1515932696U));
  ssSetChecksum2(S,(1271852803U));
  ssSetChecksum3(S,(1407760124U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c1_Ali_model(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c1_Ali_model(SimStruct *S)
{
  SFc1_Ali_modelInstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc1_Ali_modelInstanceStruct *)utMalloc(sizeof
    (SFc1_Ali_modelInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc1_Ali_modelInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c1_Ali_model;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c1_Ali_model;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c1_Ali_model;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c1_Ali_model;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c1_Ali_model;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c1_Ali_model;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c1_Ali_model;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c1_Ali_model;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c1_Ali_model;
  chartInstance->chartInfo.mdlStart = mdlStart_c1_Ali_model;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c1_Ali_model;
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

void c1_Ali_model_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c1_Ali_model(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c1_Ali_model(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c1_Ali_model(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c1_Ali_model_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
