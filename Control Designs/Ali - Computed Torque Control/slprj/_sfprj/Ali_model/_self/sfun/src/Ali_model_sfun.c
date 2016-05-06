/* Include files */

#include "Ali_model_sfun.h"
#include "Ali_model_sfun_debug_macros.h"
#include "c1_Ali_model.h"
#include "c2_Ali_model.h"
#include "c3_Ali_model.h"
#include "c4_Ali_model.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */
uint32_T _Ali_modelMachineNumber_;

/* Function Declarations */

/* Function Definitions */
void Ali_model_initializer(void)
{
}

void Ali_model_terminator(void)
{
}

/* SFunction Glue Code */
unsigned int sf_Ali_model_method_dispatcher(SimStruct *simstructPtr, unsigned
  int chartFileNumber, const char* specsCksum, int_T method, void *data)
{
  if (chartFileNumber==1) {
    c1_Ali_model_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==2) {
    c2_Ali_model_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==3) {
    c3_Ali_model_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==4) {
    c4_Ali_model_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  return 0;
}

extern void sf_Ali_model_uses_exported_functions(int nlhs, mxArray * plhs[], int
  nrhs, const mxArray * prhs[])
{
  plhs[0] = mxCreateLogicalScalar(0);
}

unsigned int sf_Ali_model_process_check_sum_call( int nlhs, mxArray * plhs[],
  int nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[20];
  if (nrhs<1 || !mxIsChar(prhs[0]) )
    return 0;

  /* Possible call to get the checksum */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"sf_get_check_sum"))
    return 0;
  plhs[0] = mxCreateDoubleMatrix( 1,4,mxREAL);
  if (nrhs>1 && mxIsChar(prhs[1])) {
    mxGetString(prhs[1], commandName,sizeof(commandName)/sizeof(char));
    commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
    if (!strcmp(commandName,"machine")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2916468723U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1091988707U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(2980474420U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(551690433U);
    } else if (!strcmp(commandName,"exportedFcn")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(0U);
    } else if (!strcmp(commandName,"makefile")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(4011605029U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3625523369U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3797628295U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(2361107479U);
    } else if (nrhs==3 && !strcmp(commandName,"chart")) {
      unsigned int chartFileNumber;
      chartFileNumber = (unsigned int)mxGetScalar(prhs[2]);
      switch (chartFileNumber) {
       case 1:
        {
          extern void sf_c1_Ali_model_get_check_sum(mxArray *plhs[]);
          sf_c1_Ali_model_get_check_sum(plhs);
          break;
        }

       case 2:
        {
          extern void sf_c2_Ali_model_get_check_sum(mxArray *plhs[]);
          sf_c2_Ali_model_get_check_sum(plhs);
          break;
        }

       case 3:
        {
          extern void sf_c3_Ali_model_get_check_sum(mxArray *plhs[]);
          sf_c3_Ali_model_get_check_sum(plhs);
          break;
        }

       case 4:
        {
          extern void sf_c4_Ali_model_get_check_sum(mxArray *plhs[]);
          sf_c4_Ali_model_get_check_sum(plhs);
          break;
        }

       default:
        ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(0.0);
      }
    } else if (!strcmp(commandName,"target")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3061339410U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1991824845U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3599338742U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(2357874978U);
    } else {
      return 0;
    }
  } else {
    ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(869050372U);
    ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2030791745U);
    ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1082868614U);
    ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3783099393U);
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_Ali_model_autoinheritance_info( int nlhs, mxArray * plhs[], int
  nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[32];
  char aiChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]) )
    return 0;

  /* Possible call to get the autoinheritance_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_autoinheritance_info"))
    return 0;
  mxGetString(prhs[2], aiChksum,sizeof(aiChksum)/sizeof(char));
  aiChksum[(sizeof(aiChksum)/sizeof(char)-1)] = '\0';

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        if (strcmp(aiChksum, "m3aY8v09NIwyuzIIOgo35E") == 0) {
          extern mxArray *sf_c1_Ali_model_get_autoinheritance_info(void);
          plhs[0] = sf_c1_Ali_model_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 2:
      {
        if (strcmp(aiChksum, "XEkbVaFWOS4tenBGP7etE") == 0) {
          extern mxArray *sf_c2_Ali_model_get_autoinheritance_info(void);
          plhs[0] = sf_c2_Ali_model_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 3:
      {
        if (strcmp(aiChksum, "rw0gTy6bB2ja2XHY3DORmF") == 0) {
          extern mxArray *sf_c3_Ali_model_get_autoinheritance_info(void);
          plhs[0] = sf_c3_Ali_model_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 4:
      {
        if (strcmp(aiChksum, "TCx0HKIQmdHBqCmrIErImD") == 0) {
          extern mxArray *sf_c4_Ali_model_get_autoinheritance_info(void);
          plhs[0] = sf_c4_Ali_model_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_Ali_model_get_eml_resolved_functions_info( int nlhs, mxArray *
  plhs[], int nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[64];
  if (nrhs<2 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the get_eml_resolved_functions_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_eml_resolved_functions_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        extern const mxArray *sf_c1_Ali_model_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c1_Ali_model_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 2:
      {
        extern const mxArray *sf_c2_Ali_model_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c2_Ali_model_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 3:
      {
        extern const mxArray *sf_c3_Ali_model_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c3_Ali_model_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 4:
      {
        extern const mxArray *sf_c4_Ali_model_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c4_Ali_model_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_Ali_model_third_party_uses_info( int nlhs, mxArray * plhs[], int
  nrhs, const mxArray * prhs[] )
{
  char commandName[64];
  char tpChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the third_party_uses_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  mxGetString(prhs[2], tpChksum,sizeof(tpChksum)/sizeof(char));
  tpChksum[(sizeof(tpChksum)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_third_party_uses_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        if (strcmp(tpChksum, "ktfoFxOo8LIrZmK0msZsIB") == 0) {
          extern mxArray *sf_c1_Ali_model_third_party_uses_info(void);
          plhs[0] = sf_c1_Ali_model_third_party_uses_info();
          break;
        }
      }

     case 2:
      {
        if (strcmp(tpChksum, "IQ19MRfQmmznYuV2EbCurE") == 0) {
          extern mxArray *sf_c2_Ali_model_third_party_uses_info(void);
          plhs[0] = sf_c2_Ali_model_third_party_uses_info();
          break;
        }
      }

     case 3:
      {
        if (strcmp(tpChksum, "L8yfMi2qPpifF1yCkuVosE") == 0) {
          extern mxArray *sf_c3_Ali_model_third_party_uses_info(void);
          plhs[0] = sf_c3_Ali_model_third_party_uses_info();
          break;
        }
      }

     case 4:
      {
        if (strcmp(tpChksum, "xriMN9rWXsigfkGwJfQQTE") == 0) {
          extern mxArray *sf_c4_Ali_model_third_party_uses_info(void);
          plhs[0] = sf_c4_Ali_model_third_party_uses_info();
          break;
        }
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;
}

unsigned int sf_Ali_model_jit_fallback_info( int nlhs, mxArray * plhs[], int
  nrhs, const mxArray * prhs[] )
{
  char commandName[64];
  char tpChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the jit_fallback_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  mxGetString(prhs[2], tpChksum,sizeof(tpChksum)/sizeof(char));
  tpChksum[(sizeof(tpChksum)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_jit_fallback_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        if (strcmp(tpChksum, "ktfoFxOo8LIrZmK0msZsIB") == 0) {
          extern mxArray *sf_c1_Ali_model_jit_fallback_info(void);
          plhs[0] = sf_c1_Ali_model_jit_fallback_info();
          break;
        }
      }

     case 2:
      {
        if (strcmp(tpChksum, "IQ19MRfQmmznYuV2EbCurE") == 0) {
          extern mxArray *sf_c2_Ali_model_jit_fallback_info(void);
          plhs[0] = sf_c2_Ali_model_jit_fallback_info();
          break;
        }
      }

     case 3:
      {
        if (strcmp(tpChksum, "L8yfMi2qPpifF1yCkuVosE") == 0) {
          extern mxArray *sf_c3_Ali_model_jit_fallback_info(void);
          plhs[0] = sf_c3_Ali_model_jit_fallback_info();
          break;
        }
      }

     case 4:
      {
        if (strcmp(tpChksum, "xriMN9rWXsigfkGwJfQQTE") == 0) {
          extern mxArray *sf_c4_Ali_model_jit_fallback_info(void);
          plhs[0] = sf_c4_Ali_model_jit_fallback_info();
          break;
        }
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;
}

unsigned int sf_Ali_model_updateBuildInfo_args_info( int nlhs, mxArray * plhs[],
  int nrhs, const mxArray * prhs[] )
{
  char commandName[64];
  char tpChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the updateBuildInfo_args_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  mxGetString(prhs[2], tpChksum,sizeof(tpChksum)/sizeof(char));
  tpChksum[(sizeof(tpChksum)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_updateBuildInfo_args_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        if (strcmp(tpChksum, "ktfoFxOo8LIrZmK0msZsIB") == 0) {
          extern mxArray *sf_c1_Ali_model_updateBuildInfo_args_info(void);
          plhs[0] = sf_c1_Ali_model_updateBuildInfo_args_info();
          break;
        }
      }

     case 2:
      {
        if (strcmp(tpChksum, "IQ19MRfQmmznYuV2EbCurE") == 0) {
          extern mxArray *sf_c2_Ali_model_updateBuildInfo_args_info(void);
          plhs[0] = sf_c2_Ali_model_updateBuildInfo_args_info();
          break;
        }
      }

     case 3:
      {
        if (strcmp(tpChksum, "L8yfMi2qPpifF1yCkuVosE") == 0) {
          extern mxArray *sf_c3_Ali_model_updateBuildInfo_args_info(void);
          plhs[0] = sf_c3_Ali_model_updateBuildInfo_args_info();
          break;
        }
      }

     case 4:
      {
        if (strcmp(tpChksum, "xriMN9rWXsigfkGwJfQQTE") == 0) {
          extern mxArray *sf_c4_Ali_model_updateBuildInfo_args_info(void);
          plhs[0] = sf_c4_Ali_model_updateBuildInfo_args_info();
          break;
        }
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;
}

void sf_Ali_model_get_post_codegen_info( int nlhs, mxArray * plhs[], int nrhs,
  const mxArray * prhs[] )
{
  unsigned int chartFileNumber = (unsigned int) mxGetScalar(prhs[0]);
  char tpChksum[64];
  mxGetString(prhs[1], tpChksum,sizeof(tpChksum)/sizeof(char));
  tpChksum[(sizeof(tpChksum)/sizeof(char)-1)] = '\0';
  switch (chartFileNumber) {
   case 1:
    {
      if (strcmp(tpChksum, "ktfoFxOo8LIrZmK0msZsIB") == 0) {
        extern mxArray *sf_c1_Ali_model_get_post_codegen_info(void);
        plhs[0] = sf_c1_Ali_model_get_post_codegen_info();
        return;
      }
    }
    break;

   case 2:
    {
      if (strcmp(tpChksum, "IQ19MRfQmmznYuV2EbCurE") == 0) {
        extern mxArray *sf_c2_Ali_model_get_post_codegen_info(void);
        plhs[0] = sf_c2_Ali_model_get_post_codegen_info();
        return;
      }
    }
    break;

   case 3:
    {
      if (strcmp(tpChksum, "L8yfMi2qPpifF1yCkuVosE") == 0) {
        extern mxArray *sf_c3_Ali_model_get_post_codegen_info(void);
        plhs[0] = sf_c3_Ali_model_get_post_codegen_info();
        return;
      }
    }
    break;

   case 4:
    {
      if (strcmp(tpChksum, "xriMN9rWXsigfkGwJfQQTE") == 0) {
        extern mxArray *sf_c4_Ali_model_get_post_codegen_info(void);
        plhs[0] = sf_c4_Ali_model_get_post_codegen_info();
        return;
      }
    }
    break;

   default:
    break;
  }

  plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
}

void Ali_model_debug_initialize(struct SfDebugInstanceStruct* debugInstance)
{
  _Ali_modelMachineNumber_ = sf_debug_initialize_machine(debugInstance,
    "Ali_model","sfun",0,4,0,0,0);
  sf_debug_set_machine_event_thresholds(debugInstance,_Ali_modelMachineNumber_,0,
    0);
  sf_debug_set_machine_data_thresholds(debugInstance,_Ali_modelMachineNumber_,0);
}

void Ali_model_register_exported_symbols(SimStruct* S)
{
}

static mxArray* sRtwOptimizationInfoStruct= NULL;
mxArray* load_Ali_model_optimization_info(void)
{
  if (sRtwOptimizationInfoStruct==NULL) {
    sRtwOptimizationInfoStruct = sf_load_rtw_optimization_info("Ali_model",
      "Ali_model");
    mexMakeArrayPersistent(sRtwOptimizationInfoStruct);
  }

  return(sRtwOptimizationInfoStruct);
}

void unload_Ali_model_optimization_info(void)
{
  if (sRtwOptimizationInfoStruct!=NULL) {
    mxDestroyArray(sRtwOptimizationInfoStruct);
    sRtwOptimizationInfoStruct = NULL;
  }
}
