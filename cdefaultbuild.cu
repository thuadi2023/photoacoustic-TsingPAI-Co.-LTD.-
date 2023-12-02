#include "cuda_runtime.h"  
#include "device_launch_parameters.h"    
#include "stdafx.h"
#include "cdefaultbuild.cuh"


__global__ void recKernel(GPUMemory gpu_memory, RecoParam reco_param)
{

	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = j * 512 + i;

	float dx = 0;
	float dy = 0;
	float rr0 = 0;
	float idx_temp = 0;
	float pa_temp = 0;
	float new_a = 0;

	float Vs = reco_param.voice_speed;
	float Fr = reco_param.sampe_fr;
	int Tp = reco_param.time_point;
	float R = reco_param.radius;


	float pa = 0;
	float p = 0;
	float pa1 = 0;
	float pa2 = 0;
	float pa3 = 0;
	float pa4 = 0;
	
	int int_idx = 0;
	
      float pi = 3.14159265359;
	
	float paDiff;
	float paBP;
	float angleCorrection;
	float dOmega;


	for (int iStep = 0; iStep < reco_param.channel; iStep++)
	{
		dx = (k % 512 + 1 - 512 / 2.0)*reco_param.pixel_size - (iStep - reco_param.channel / 2) * reco_param.pitch;
		dy = ((512 - k / 512) - 512 / 2.0)*reco_param.pixel_size - 512 / 2 * reco_param.pixel_size;
		if (fabs(dx) < fabs(dy))
		{
			rr0 = sqrt(dx * dx + dy * dy);
			idx_temp = rr0 / Vs * Fr+74;
			if (idx_temp > (Tp - 1 - 4))
			{
				idx_temp = Tp - 1 - 4;
			}
			else if (idx_temp < 1)
			{
				idx_temp = 1;
			}

			int_idx = int(idx_temp);
			pa1 = gpu_memory.Rawdata[(int_idx)+iStep * Tp  - 1];
			pa2 = gpu_memory.Rawdata[(int_idx)+iStep * Tp ];
			pa23 = gpu_memory.Rawdata[(int_idx)+iStep * Tp + 1 ];
			paDiff = pa2 - pa1;
			paBP = 2 * pa1;
			rr0 = rr0 + 0.001;
			dOmega = 0.25*1e-3*0.25*1e-3* fabs(dy) / (rr0 * rr0 * rr0);
			angleCorrection = dOmega / (2 * pi);
			pa = pa + pa2;
		}
		else
		{
			pa = pa + 0;
		}
	}



	new_a = (((exp(((reco_param.scanDepth[k / 512 / 96 + 2] - reco_param.scanDepth[k / 512 / 96 + 1])/96.0*(k / 512 % 96)+reco_param.scanDepth[k / 512 / 96 + 1])*0.1) + std::exp(j / 150.0)))* exp((reco_param.amp - 50)*0.2))*0.1;
	pa_temp = pa * new_a / 200 + (reco_param.dynamicRange - 50) * 5;
	if (pa_temp > 255)
		pa_temp = 255;
	else if (pa_temp < 0)
		pa_temp = 0;
	gpu_memory.ImageArray[k] = uchar(pa_temp);
}

