#include "c_backend_api.h"

TVMValue param0[2];
TVMValue param1[3];
TVMValue param2[2];
TVMArray a0[1];
TVMArray a1[1];
TVMArray b0[1];
TVMArray b1[1];
TVMArray b2[1];
TVMArray c1[1];
TVMArray c2[1];


__attribute__((noinline))
void kernel1(int32_t ax0_ax1_fused_ax2_fused, float* placeholder, float* T_layout_trans)
{
    int32_t ax3;
    for (ax3 = 0; ax3 < 64; ++ax3) 
    {
      T_layout_trans[((ax0_ax1_fused_ax2_fused * 64) + ax3)] = placeholder[((ax0_ax1_fused_ax2_fused * 64) + ax3)];
    }
}

__attribute__((noinline))
void kernel2(int32_t ax0_ax1_fused, float* placeholder, float* T_layout_trans)
{
  int32_t ax2;
  int32_t ax3;
  
    for (ax2 = 0; ax2 < 64; ++ax2) {
      for (ax3 = 0; ax3 < 64; ++ax3) {
        T_layout_trans[(((ax0_ax1_fused * 4096) + (ax2 * 64)) + ax3)] = placeholder[(((ax2 * 128) + (ax3 * 2)) + ax0_ax1_fused)];
      }
    }
}

__attribute__((noinline))
void kernel3(int32_t i1_i2_fused, float *data_pad, float* placeholder)
{
  int32_t i3;
  for (i3 = 0; i3 < 66; ++i3) {
    (( float*)data_pad)[((i1_i2_fused * 66) + i3)] = (((((1 <= i1_i2_fused) && (i1_i2_fused < 65)) && (1 <= i3)) && (i3 < 65)) ? placeholder[(((i1_i2_fused * 64) + i3) - 65)] : 0.000000e+00f);
    }
}

__attribute__((noinline))
void kernel4(int32_t n_oc_chunk_fused_oh_fused, float *data_pad, float* placeholder1, float* conv2d_NCHWc)
{
         float2 conv2d_NCHWc_global[16];
    for (ow_outer = 0; ow_outer < 4; ++ow_outer) {
      conv2d_NCHWc_global[0] = ((float2)(0.000000e+00f, 0.000000e+00f));
      conv2d_NCHWc_global[1] = ((float2)(0.000000e+00f, 0.000000e+00f));
      conv2d_NCHWc_global[2] = ((float2)(0.000000e+00f, 0.000000e+00f));
      conv2d_NCHWc_global[3] = ((float2)(0.000000e+00f, 0.000000e+00f));
      conv2d_NCHWc_global[4] = ((float2)(0.000000e+00f, 0.000000e+00f));
      conv2d_NCHWc_global[5] = ((float2)(0.000000e+00f, 0.000000e+00f));
      conv2d_NCHWc_global[6] = ((float2)(0.000000e+00f, 0.000000e+00f));
      conv2d_NCHWc_global[7] = ((float2)(0.000000e+00f, 0.000000e+00f));
      conv2d_NCHWc_global[8] = ((float2)(0.000000e+00f, 0.000000e+00f));
      conv2d_NCHWc_global[9] = ((float2)(0.000000e+00f, 0.000000e+00f));
      conv2d_NCHWc_global[10] = ((float2)(0.000000e+00f, 0.000000e+00f));
      conv2d_NCHWc_global[11] = ((float2)(0.000000e+00f, 0.000000e+00f));
      conv2d_NCHWc_global[12] = ((float2)(0.000000e+00f, 0.000000e+00f));
      conv2d_NCHWc_global[13] = ((float2)(0.000000e+00f, 0.000000e+00f));
      conv2d_NCHWc_global[14] = ((float2)(0.000000e+00f, 0.000000e+00f));
      conv2d_NCHWc_global[15] = ((float2)(0.000000e+00f, 0.000000e+00f));
      for (kh = 0; kh < 3; ++kh) {
        for (kw = 0; kw < 3; ++kw) {
          conv2d_NCHWc_global[0] = (conv2d_NCHWc_global[0] + (((float2)((( float*)data_pad)[((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw)], (( float*)data_pad)[((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw)])) * (( float2*)(placeholder1 + ((kh * 6) + (kw * 2))))[0]));
          conv2d_NCHWc_global[1] = (conv2d_NCHWc_global[1] + (((float2)((( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 1)], (( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 1)])) * (( float2*)(placeholder1 + ((kh * 6) + (kw * 2))))[0]));
          conv2d_NCHWc_global[2] = (conv2d_NCHWc_global[2] + (((float2)((( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 2)], (( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 2)])) * (( float2*)(placeholder1 + ((kh * 6) + (kw * 2))))[0]));
          conv2d_NCHWc_global[3] = (conv2d_NCHWc_global[3] + (((float2)((( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 3)], (( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 3)])) * (( float2*)(placeholder1 + ((kh * 6) + (kw * 2))))[0]));
          conv2d_NCHWc_global[4] = (conv2d_NCHWc_global[4] + (((float2)((( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 4)], (( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 4)])) * (( float2*)(placeholder1 + ((kh * 6) + (kw * 2))))[0]));
          conv2d_NCHWc_global[5] = (conv2d_NCHWc_global[5] + (((float2)((( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 5)], (( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 5)])) * (( float2*)(placeholder1 + ((kh * 6) + (kw * 2))))[0]));
          conv2d_NCHWc_global[6] = (conv2d_NCHWc_global[6] + (((float2)((( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 6)], (( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 6)])) * (( float2*)(placeholder1 + ((kh * 6) + (kw * 2))))[0]));
          conv2d_NCHWc_global[7] = (conv2d_NCHWc_global[7] + (((float2)((( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 7)], (( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 7)])) * (( float2*)(placeholder1 + ((kh * 6) + (kw * 2))))[0]));
          conv2d_NCHWc_global[8] = (conv2d_NCHWc_global[8] + (((float2)((( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 8)], (( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 8)])) * (( float2*)(placeholder1 + ((kh * 6) + (kw * 2))))[0]));
          conv2d_NCHWc_global[9] = (conv2d_NCHWc_global[9] + (((float2)((( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 9)], (( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 9)])) * (( float2*)(placeholder1 + ((kh * 6) + (kw * 2))))[0]));
          conv2d_NCHWc_global[10] = (conv2d_NCHWc_global[10] + (((float2)((( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 10)], (( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 10)])) * (( float2*)(placeholder1 + ((kh * 6) + (kw * 2))))[0]));
          conv2d_NCHWc_global[11] = (conv2d_NCHWc_global[11] + (((float2)((( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 11)], (( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 11)])) * (( float2*)(placeholder1 + ((kh * 6) + (kw * 2))))[0]));
          conv2d_NCHWc_global[12] = (conv2d_NCHWc_global[12] + (((float2)((( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 12)], (( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 12)])) * (( float2*)(placeholder1 + ((kh * 6) + (kw * 2))))[0]));
          conv2d_NCHWc_global[13] = (conv2d_NCHWc_global[13] + (((float2)((( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 13)], (( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 13)])) * (( float2*)(placeholder1 + ((kh * 6) + (kw * 2))))[0]));
          conv2d_NCHWc_global[14] = (conv2d_NCHWc_global[14] + (((float2)((( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 14)], (( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 14)])) * (( float2*)(placeholder1 + ((kh * 6) + (kw * 2))))[0]));
          conv2d_NCHWc_global[15] = (conv2d_NCHWc_global[15] + (((float2)((( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 15)], (( float*)data_pad)[(((((kh * 66) + (n_oc_chunk_fused_oh_fused * 66)) + (ow_outer * 16)) + kw) + 15)])) * (( float2*)(placeholder1 + ((kh * 6) + (kw * 2))))[0]));
        }
      }
      }
      int32_t ow_inner;
          for (ow_inner = 0; ow_inner < 16; ++ow_inner) {
        (( float2*)(conv2d_NCHWc + (((n_oc_chunk_fused_oh_fused * 128) + (ow_outer * 32)) + (ow_inner * 2))))[0] = (( float2*)(( float*)conv2d_NCHWc_global + (ow_inner * 2)))[0];
      }

}


__attribute__((noinline))
void parallel1(float* placeholder, float* T_layout_trans)
{
  int32_t ax0_ax1_fused_ax2_fused;

  #pragma omp parallel for
  for (ax0_ax1_fused_ax2_fused = 0; ax0_ax1_fused_ax2_fused < 64; ++ax0_ax1_fused_ax2_fused)
  {
    kernel1(ax0_ax1_fused_ax2_fused, placeholder, T_layout_trans);
  }
}

__attribute__((noinline))
void parallel2(float* placeholder, float* T_layout_trans)
{
  int32_t ax0_ax1_fused;
  #pragma omp parallel for
  for (ax0_ax1_fused = 0; ax0_ax1_fused < 2; ++ax0_ax1_fused)
  {
    kernel2(ax0_ax1_fused, placeholder, T_layout_trans);
  }
}

__attribute__((noinline))
void parallel3(float *data_pad, float* placeholder)
{
  int32_t i1_i2_fused;
  #pragma omp parallel for
  for (i1_i2_fused = 0; i1_i2_fused < 66; ++i1_i2_fused)
  {
    kernel3(i1_i2_fused, data_pad, placeholder);
  }
}

__attribute__((noinline))
void parallel4(float *data_pad, float* placeholder1, float* conv2d_NCHWc)
{
  int32_t n_oc_chunk_fused_oh_fused;
  #pragma omp parallel for
  for (n_oc_chunk_fused_oh_fused = 0; n_oc_chunk_fused_oh_fused < 64; ++n_oc_chunk_fused_oh_fused)
  {
    kernel4(n_oc_chunk_fused_oh_fused, data_pad, placeholder1, conv2d_NCHWc);
  }
}

int32_t fused_layout_transform_2( void* args,  void* arg_type_ids, int32_t num_args) {
  void* arg0 = (((TVMValue*)args)[0].v_handle);
  void* arg1 = (((TVMValue*)args)[1].v_handle);
  float* placeholder = (float*)(((TVMArray*)arg0)[0].data);
  float* T_layout_trans = (float*)(((TVMArray*)arg1)[0].data);
  
  parallel1(placeholder, T_layout_trans);

  return 0;
}

int32_t fused_layout_transform_1( void* args,  void* arg_type_ids, int32_t num_args) {
  void* arg0 = (((TVMValue*)args)[0].v_handle);
  void* arg1 = (((TVMValue*)args)[1].v_handle);
  float* placeholder = (float*)(((TVMArray*)arg0)[0].data);
  float* T_layout_trans = (float*)(((TVMArray*)arg1)[0].data);

  parallel2(placeholder, T_layout_trans);
  
  return 0;
}

int32_t fused_nn_contrib_conv2d_NCHWc( void* args,  void* arg_type_ids, int32_t num_args) {
  void* arg0 = (((TVMValue*)args)[0].v_handle);
  void* arg1 = (((TVMValue*)args)[1].v_handle);
  void* arg2 = (((TVMValue*)args)[2].v_handle);
  float* placeholder = (float*)(((TVMArray*)arg0)[0].data);
  float* placeholder1 = (float*)(((TVMArray*)arg1)[0].data);
  float* conv2d_NCHWc = (float*)(((TVMArray*)arg2)[0].data);
  
  void* data_pad = TVMBackendAllocWorkspace(1, dev_id, (uint64_t)17424, 2, 32);
  if (data_pad == NULL) {
    return -1;
  }
  
  parallel3(data_pad, placeholder);
  
  parallel4(data_pad, placeholder1, conv2d_NCHWc);

  if (TVMBackendFreeWorkspace(1, dev_id, data_pad) != 0) {
    return -1;
  }
  return 0;
}

int32_t conv(TVMValue* param0, TVMValue param1, TVMValue param2){
  int res1, res2, res3;
  res1 = fused_layout_transform_2(param0, 0, 0);
  res2 = fused_nn_contrib_conv2d_NCHWc(param1, 0, 0);
  res3 = fused_layout_transform_1(param2, 0, 0);
  return res3;
}

int32_t fused_conv2d_wrapper(float* X, float* p0, float* T_fused_layout_transform_1) {

  float out_fused_layout_transform_2[4096];
  float out_conv[8192];
  int32_t res;
  a0[0].data = X;
  a1[0].data = out_fused_layout_transform_2;
  param0[0].v_handle = a0;
  param0[1].v_handle = a1;
  
  b0[0].data = out_fused_layout_transform_2;
  b1[0].data = p0;
  b2[0].data = out_conv;
  param1[0].v_handle = b0;
  param1[1].v_handle = b1;
  param1[2].v_handle = b2;
  
  c0[0].data = out_conv;
  c1[0].data = T_fused_layout_transform_1;
  param2[0].v_handle = c0;
  param2[1].v_handle = c1;
  
  
#ifdef BAMBU_PROFILING
  __builtin_bambu_time_start();
#endif

  res = conv(param0, param1, param2);
    
#ifdef BAMBU_PROFILING
  __builtin_bambu_time_stop();
#endif
    
    return res;
}
