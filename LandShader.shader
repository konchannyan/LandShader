//Copyright(c) <2019> <JackyGun twitter@konchannyan>
//Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files(the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and / or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions :
//
//The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Shader "JackyGun/LeadShader/land"
{
	SubShader
	{
		Tags{ "Queue" = "Geometry+1" "RenderType" = "Opaque" }

		Pass
	{
		CGPROGRAM
#pragma target 5.0
#pragma vertex mainVS
#pragma hull mainHS
#pragma domain mainDS
#pragma geometry mainGS
#pragma fragment mainFS
#pragma multi_compile_instancing

#include "UnityCG.cginc"

#define TESS 20
#define GP_INSTANCE_SIZE 12
#define GP_INSTANCE_NUM 12
#define GP_INSTANCE_START 4
#define GP_INSTANCE_END 12
#define RAIL_INSTANCE_START 6
#define RAIL_INSTANCE_END 12

	// Structure
	struct VS_IN
	{
		float4 pos       : POSITION;
		uint  instanceID : SV_InstanceID;
	};

	struct VS_OUT
	{
		float4 pos       : POSITION;
		uint  instanceID : SV_InstanceID;
		uint4 ids0       : IDS0;
		uint4 ids1       : IDS1;
		uint4 ids2       : IDS2;
	};

	struct CONSTANT_HS_OUT
	{
		float Edges[4] : SV_TessFactor;
		float Inside[2] : SV_InsideTessFactor;
	};

	struct HS_OUT
	{
		uint  instanceID : SV_InstanceID;
		uint4 ids0       : IDS0;
		uint4 ids1       : IDS1;
		uint4 ids2       : IDS2;
	};

	struct DS_OUT
	{
		uint  tid        : TID;
		uint  instanceID : SV_InstanceId;
		uint4 ids0	     : IDS0;
		uint4 ids1	     : IDS1;
		uint4 ids2	     : IDS2;
	};

	struct GS_OUT
	{
		float4 vertex : SV_POSITION;
		float3 color0  : COLOR0;
		float3 color1  : COLOR1;
		float3 normal : NORMAL;
		float3 light   : LIGHT;
	};

	float rand(float2 p, float seed)
	{
		return frac(sin(dot(p, fixed2(12.9898, 78.233)) + seed) * 43758.5453);
	}

	float random(float2 st, float seed)
	{
		float2 p = floor(st);
		float2 f = frac(st);
		float2 u = f*f*(3.0 - 2.0*f);

		float v00 = rand(p + float2(0, 0), seed);
		float v10 = rand(p + float2(1, 0), seed);
		float v01 = rand(p + float2(0, 1), seed);
		float v11 = rand(p + float2(1, 1), seed);

		return lerp(lerp(dot(v00, f - float2(0, 0)), dot(v10, f - float2(1, 0)), u.x),
			lerp(dot(v01, f - float2(0, 1)), dot(v11, f - float2(1, 1)), u.x),
			u.y) + 0.5f;
	}

	//float4 PseudoColor(float a, float min, float max)
	//{
	//	float x = (a - min) / (max - min);
	//	x = saturate(x);

	//	float r = (x<0.50) ? 0 : ((x<0.75) ? (x - 0.5) * 4 : 1);
	//	float g = (x<0.25) ? x * 4 : (x<0.75) ? 1 : (1 - x) * 4;
	//	float b = (x<0.25) ? 1 : (x<0.5) ? (0.5 - x) * 4 : 0;

	//	return float4(r, g, b, 1);
	//}

	float3x3 MatRollX(float a)
	{
		return float3x3( 1, 0, 0, 0, cos(a), -sin(a), 0, sin(a), cos(a) );
	}

	float3x3 MatRollY(float a)
	{
		return float3x3( sin(a), 0, cos(a), 0, 1, 0, cos(a), 0, -sin(a) );
	}

	//float3x3 MatRollZ(float a)
	//{
	//	return float3x3( cos(a), -sin(a), 0, sin(a), cos(a), 0, 0, 0, 1 );
	//}

	// Main
	VS_OUT mainVS(VS_IN In)
	{

		VS_OUT Out;
		Out.pos = In.pos;
		Out.instanceID = In.instanceID;
		Out.ids0 = 0;
		Out.ids1 = 0;
		Out.ids2 = 0;

#if defined(UNITY_INSTANCING_ENABLED)

		VS_IN dm;

		uint i, j;

		uint index[GP_INSTANCE_SIZE];
		float scale[GP_INSTANCE_SIZE];

		for (i = 0; i < GP_INSTANCE_SIZE; i++)
		{
			dm.instanceID = i;
			UNITY_SETUP_INSTANCE_ID(dm);
			scale[i] = length(unity_ObjectToWorld[0].xyz);
		}

		// bad sort	
		for (i = 0; i < GP_INSTANCE_SIZE; i++)
		{
			index[i] = i;
		}
		for (i = 0; i < (GP_INSTANCE_SIZE - 1); i++)
		{
			uint b = i;
			for (j = i + 1; j < GP_INSTANCE_SIZE; j++)
			{
				if (scale[index[b]] > scale[index[j]])
				{
					b = j;
				}
			}
			uint s = index[i];
			index[i] = index[b];
			index[b] = s;
		}

		Out.ids0 = uint4(index[ 0], index[ 1], index[ 2], index[ 3]);
		Out.ids1 = uint4(index[ 4], index[ 5], index[ 6], index[ 7]);
		Out.ids2 = uint4(index[ 8], index[ 9], index[10], index[11]);

#endif
		
		return Out;

	}

	CONSTANT_HS_OUT mainCHS()
	{
		CONSTANT_HS_OUT Out;

		float d = distance(unity_ObjectToWorld._14_24_34 / unity_ObjectToWorld._44, _WorldSpaceCameraPos);

		int t = d > 10 ? 0 : TESS + 1;
		Out.Edges[0] = t;
		Out.Edges[1] = t;
		Out.Edges[2] = t;
		Out.Edges[3] = t;
		Out.Inside[0] = t;
		Out.Inside[1] = t;

		return Out;
	}

	[domain("quad")]
	[partitioning("pow2")]
	[outputtopology("point")]
	[outputcontrolpoints(4)]
	[patchconstantfunc("mainCHS")]
	HS_OUT mainHS(InputPatch<VS_OUT, 1> In)
	{
		HS_OUT Out;
		Out.instanceID = In[0].instanceID;
		Out.ids0 = In[0].ids0;
		Out.ids1 = In[0].ids1;
		Out.ids2 = In[0].ids2;
		return Out;
	}

	[domain("quad")]
	DS_OUT mainDS(CONSTANT_HS_OUT In, const OutputPatch<HS_OUT, 4> patch, float2 uv : SV_DomainLocation)
	{
		DS_OUT Out;
		Out.tid = (uint)(uv.x * TESS) + ((uint)(uv.y * TESS) * TESS);
		Out.instanceID = patch[0].instanceID;
		Out.ids0 = patch[0].ids0;
		Out.ids1 = patch[0].ids1;
		Out.ids2 = patch[0].ids2;
		return Out;
	}

#define GS_INSTANCE_NUM 1
	[maxvertexcount(64)]
	[instance(GS_INSTANCE_NUM)]
	void mainGS(point DS_OUT input[1], inout TriangleStream<GS_OUT> outStream, uint gid : SV_GSInstanceID)
	{

#if defined(UNITY_INSTANCING_ENABLED)

		GS_OUT o = (GS_OUT)0;
		
		int i;

		// ids
		uint iid = input[0].instanceID;
		uint tid = input[0].tid;
		uint pid = tid * GS_INSTANCE_NUM + gid;
		uint  id = pid * GP_INSTANCE_NUM + iid;
		
		uint ids[GP_INSTANCE_SIZE] = {
			input[0].ids0.x, input[0].ids0.y, input[0].ids0.z, input[0].ids0.w,
			input[0].ids1.x, input[0].ids1.y, input[0].ids1.z, input[0].ids1.w,
			input[0].ids2.x, input[0].ids2.y, input[0].ids2.z, input[0].ids2.w,
		};

		// fetch
		float3 po[GP_INSTANCE_SIZE] = { float3(0, 0, 0), float3(0, 0, 0), float3(0, 0, 0), float3(0, 0, 0), float3(0, 0, 0), float3(0, 0, 0), float3(0, 0, 0), float3(0, 0, 0), float3(0, 0, 0), float3(0, 0, 0), float3(0, 0, 0), float3(0, 0, 0)};
		for (i = 0; i < GP_INSTANCE_SIZE; i++)
		{
			DS_OUT dm = (DS_OUT)0;
			dm.instanceID = ids[i % GP_INSTANCE_SIZE];
			UNITY_SETUP_INSTANCE_ID(dm);
			po[i] = unity_ObjectToWorld._14_24_34 / unity_ObjectToWorld._44;
		}
		float3 bp = po[0];
		//float3 cp = 0;
		//for (i = GP_INSTANCE_START; i < GP_INSTANCE_END; i++)
		//{
		//	cp += po[i];
		//}
		//cp /= (GP_INSTANCE_END - GP_INSTANCE_START);

		// lsm
		const float spdSM = 0.01f;
		const float radSM = 0.15f;
		float angS = _Time.y * 2 * 3.14159265f * spdSM + 0.00000000f;
		float angM = _Time.y * 2 * 3.14159265f * spdSM + 3.14159265f;
		float3 posS = bp + normalize(float3(cos(angS), sin(angS), cos(angS))) * 1.8f;
		float3 posM = bp + normalize(float3(cos(angM), sin(angM), cos(angM))) * 1.8f;
		o.light = normalize(posS - bp);
		
		// sm
		const uint divSM0 = 12;
		const uint divSM1 = divSM0 / 2;
		float angSM00 = ((id + 0) % divSM0) / (float)divSM0 * 3.14159265f * 2;
		float angSM01 = ((id + 1) % divSM0) / (float)divSM0 * 3.14159265f * 2;
		float angSM10 = ((id / divSM0 + 0)) / (float)divSM1 * 3.14159265f + 3.14159265f / 2;
		float angSM11 = ((id / divSM0 + 1)) / (float)divSM1 * 3.14159265f + 3.14159265f / 2;
		float3 posSM00 = float3(cos(angSM10) * cos(angSM00), sin(angSM10), cos(angSM10) * sin(angSM00)) * 0.15f;
		float3 posSM01 = float3(cos(angSM10) * cos(angSM01), sin(angSM10), cos(angSM10) * sin(angSM01)) * 0.15f;
		float3 posSM10 = float3(cos(angSM11) * cos(angSM00), sin(angSM11), cos(angSM11) * sin(angSM00)) * 0.15f;
		float3 posSM11 = float3(cos(angSM11) * cos(angSM01), sin(angSM11), cos(angSM11) * sin(angSM01)) * 0.15f;
		if (id < divSM0 * divSM1)
		{
			o.color0 = float3(0.995f, 0.125f, 0.023f);
			o.vertex = mul(UNITY_MATRIX_VP, float4(posS + posSM00, 1));
			outStream.Append(o);
			o.vertex = mul(UNITY_MATRIX_VP, float4(posS + posSM01, 1));
			outStream.Append(o);
			o.vertex = mul(UNITY_MATRIX_VP, float4(posS + posSM10, 1));
			outStream.Append(o);
			o.vertex = mul(UNITY_MATRIX_VP, float4(posS + posSM11, 1));
			outStream.Append(o);
			outStream.RestartStrip();
			o.color0 = float3(0.925f, 0.995f, 0.123f);
			o.color1 = o.color0;
			o.normal = normalize(posSM00);
			o.vertex = mul(UNITY_MATRIX_VP, float4(posM + posSM00, 1));
			outStream.Append(o);
			o.normal = normalize(posSM01);
			o.vertex = mul(UNITY_MATRIX_VP, float4(posM + posSM01, 1));
			outStream.Append(o);
			o.normal = normalize(posSM10);
			o.vertex = mul(UNITY_MATRIX_VP, float4(posM + posSM10, 1));
			outStream.Append(o);
			o.normal = normalize(posSM11);
			o.vertex = mul(UNITY_MATRIX_VP, float4(posM + posSM11, 1));
			outStream.Append(o);
			outStream.RestartStrip();
		}

		const float3 normalC[6] = { float3(0, +1, 0), float3(+1, 0, 0), float3(-1, 0, 0), float3(0, 0, +1), float3(0, 0, -1), float3(0, -1, 0) };
		const float3 idxC[6 * 4] =
		{
			float3(-1, +1, +1), float3(+1, +1, +1), float3(-1, +1, -1), float3(+1, +1, -1),
			float3(+1, +1, -1), float3(+1, +1, +1), float3(+1, -1, -1), float3(+1, -1, +1),
			float3(-1, +1, +1), float3(-1, +1, -1), float3(-1, -1, +1), float3(-1, -1, -1),
			float3(+1, +1, +1), float3(-1, +1, +1), float3(+1, -1, +1), float3(-1, -1, +1),
			float3(-1, +1, -1), float3(+1, +1, -1), float3(-1, -1, -1), float3(+1, -1, -1),
			float3(-1, -1, -1), float3(+1, -1, -1), float3(-1, -1, +1), float3(+1, -1, +1),
		};

		// p
		const uint lenP = GP_INSTANCE_END - GP_INSTANCE_START;
		const float sizeP = 0.005f;
		uint idP0 = id / 6 % lenP;
		uint idP1 = id % 6;
		float3 pP = po[idP0 + GP_INSTANCE_START];
		if (id < lenP * 6)
		{
			o.color0 = float3(0.99f, 0.01f, 0.01f);
			o.color1 = o.color0;
			o.normal = 0;
			o.vertex = mul(UNITY_MATRIX_VP, float4(pP + sizeP * idxC[idP1 * 4 + 0], 1));
			outStream.Append(o);
			o.vertex = mul(UNITY_MATRIX_VP, float4(pP + sizeP * idxC[idP1 * 4 + 1], 1));
			outStream.Append(o);
			o.vertex = mul(UNITY_MATRIX_VP, float4(pP + sizeP * idxC[idP1 * 4 + 2], 1));
			outStream.Append(o);
			o.vertex = mul(UNITY_MATRIX_VP, float4(pP + sizeP * idxC[idP1 * 4 + 3], 1));
			outStream.Append(o);
			outStream.RestartStrip();
		}

		// w
		const uint divW = 64;
		const float scaleW = 1.8f;
		const float sW0 = 0.5f;
		const float sW1 = 0.2f;
		const float yW = -2.0f;
		float xW0 = (id % divW + 0) / (float)(divW - 1) - 0.5f;
		float xW1 = (id % divW + 1) / (float)(divW - 1) - 0.5f;
		float zW0 = (id / divW + 0) / (float)(divW - 1) - 0.5f;
		float zW1 = (id / divW + 1) / (float)(divW - 1) - 0.5f;
		float3 posW00 = bp + float3(0, yW, 0) + float3(xW0, 0, zW0) * 2 * scaleW;
		float3 posW01 = bp + float3(0, yW, 0) + float3(xW1, 0, zW0) * 2 * scaleW;
		float3 posW10 = bp + float3(0, yW, 0) + float3(xW0, 0, zW1) * 2 * scaleW;
		float3 posW11 = bp + float3(0, yW, 0) + float3(xW1, 0, zW1) * 2 * scaleW;
		posW00.y += random(posW00.xz + _Time.y * 0.1f, 0) * sW0 + random(posW00.xz * 5 + _Time.y * 0.4f, 0) * sW1;
		posW01.y += random(posW01.xz + _Time.y * 0.1f, 0) * sW0 + random(posW01.xz * 5 + _Time.y * 0.4f, 0) * sW1;
		posW10.y += random(posW10.xz + _Time.y * 0.1f, 0) * sW0 + random(posW10.xz * 5 + _Time.y * 0.4f, 0) * sW1;
		posW11.y += random(posW11.xz + _Time.y * 0.1f, 0) * sW0 + random(posW11.xz * 5 + _Time.y * 0.4f, 0) * sW1;
		if (id < divW * divW)
		{
			o.color0 = float3(0.115f, 0.275f, 0.853f);
			o.color1 = float3(0.555f, 0.625f, 0.923f);
			o.normal = normalize(cross(posW11 - posW10, posW00 - posW11));
			o.vertex = mul(UNITY_MATRIX_VP, float4(posW10, 1));
			outStream.Append(o);
			o.vertex = mul(UNITY_MATRIX_VP, float4(posW11, 1));
			outStream.Append(o);
			o.vertex = mul(UNITY_MATRIX_VP, float4(posW00, 1));
			outStream.Append(o);
			o.vertex = mul(UNITY_MATRIX_VP, float4(posW01, 1));
			outStream.Append(o);
			outStream.RestartStrip();
		}

		// land
		const uint divL = 64;
		const float scaleL = 1.8f;
		const float yL = -1.25f;
		const float radL = 0.85f;
		const float minL = 0.9f;
		float xL0 = (id % divL + 0) / (float)(divL - 1) - 0.5f;
		float xL1 = (id % divL + 1) / (float)(divL - 1) - 0.5f;
		float zL0 = (id / divL + 0) / (float)(divL - 1) - 0.5f;
		float zL1 = (id / divL + 1) / (float)(divL - 1) - 0.5f;
		float3 posL00 = bp + float3(0, yL, 0) + float3(xL0, 0, zL0) * 2 * scaleL;
		float3 posL01 = bp + float3(0, yL, 0) + float3(xL1, 0, zL0) * 2 * scaleL;
		float3 posL10 = bp + float3(0, yL, 0) + float3(xL0, 0, zL1) * 2 * scaleL;
		float3 posL11 = bp + float3(0, yL, 0) + float3(xL1, 0, zL1) * 2 * scaleL;
		posL00.y += distance(posL00.xz, bp.xz) * -0.3f + random(posL00.xz, 0) * 0.2f + random(posL00.xz * 3, 0) * 0.03f;
		posL01.y += distance(posL01.xz, bp.xz) * -0.3f + random(posL01.xz, 0) * 0.2f + random(posL01.xz * 3, 0) * 0.03f;
		posL10.y += distance(posL10.xz, bp.xz) * -0.3f + random(posL10.xz, 0) * 0.2f + random(posL10.xz * 3, 0) * 0.03f;
		posL11.y += distance(posL11.xz, bp.xz) * -0.3f + random(posL11.xz, 0) * 0.2f + random(posL11.xz * 3, 0) * 0.03f;
		for (i = GP_INSTANCE_START; i < GP_INSTANCE_END; i++)
		{
			float dL;
			dL = sqrt(distance(posL00.y, po[i].y)) * radL;
			posL00.y = lerp(posL00.y, po[i].y, min(minL, sqrt(saturate((dL - distance(po[i].xz, posL00.xz)) / dL))));
			dL = sqrt(distance(posL01.y, po[i].y)) * radL;
			posL01.y = lerp(posL01.y, po[i].y, min(minL, sqrt(saturate((dL - distance(po[i].xz, posL01.xz)) / dL))));
			dL = sqrt(distance(posL10.y, po[i].y)) * radL;
			posL10.y = lerp(posL10.y, po[i].y, min(minL, sqrt(saturate((dL - distance(po[i].xz, posL10.xz)) / dL))));
			dL = sqrt(distance(posL11.y, po[i].y)) * radL;
			posL11.y = lerp(posL11.y, po[i].y, min(minL, sqrt(saturate((dL - distance(po[i].xz, posL11.xz)) / dL))));
		}
		if (id < divL * divL)
		{
			o.color0 = float3(0.315f, 0.875f, 0.253f);
			o.color1 = float3(0.213f, 0.923f, 0.318f);
			o.normal = normalize(cross(posL11 - posL10, posL00 - posL11));
			o.vertex = mul(UNITY_MATRIX_VP, float4(posL10, 1));
			outStream.Append(o);
			o.vertex = mul(UNITY_MATRIX_VP, float4(posL11, 1));
			outStream.Append(o);
			o.vertex = mul(UNITY_MATRIX_VP, float4(posL00, 1));
			outStream.Append(o);
			o.vertex = mul(UNITY_MATRIX_VP, float4(posL01, 1));
			outStream.Append(o);
			outStream.RestartStrip();
		}

		// cloud
		const uint divC = 16;
		const float yC = 1.6f;
		const float scaleC = 1.8f;
		const float3 spdC = float3(0.2f, 0.0f, 0.1f) * 5.0f;
		uint idC0 = id / 6;
		uint idC1 = id % 6;
		float seedC = (float)idC0 / (divC * divC);
		float3 sizeC = (float3(rand(float2(seedC, 0.1f), 0), rand(float2(seedC, 0.2f), 0), rand(float2(seedC, 0.3f), 0)) + 0.3f) * (rand(float2(seedC, 0.4f), 0) + 0.3f) * 0.2f;
		float3 posCb = (float3(random(float2(seedC, seedC) * 300, 0), random(float2(seedC, seedC) * 300, 1), random(float2(seedC, seedC) * 300, 2))) * 2 * float3(4.0f, 0.2f, 4.0f) + float3(-4.0f, 0.0, -4.0f) + float3(4.0f, 0.0f, 4.0f) * 2 * frac(_Time.x);
		posCb = (frac((posCb + 4) * 0.125f) - 0.5f) * 2 * 2 + float3(0, yC, 0);

		if (id < divC * divC * 6)
		{
			o.color0 = float3(0.99f, 0.99f, 0.99f);
			o.color1 = o.color0;
			o.normal = normalize(lerp(float3(0, 1, 0), normalC[idC1], 0.25f));
			float3 pC;
			pC = posCb + sizeC * idxC[idC1 * 4 + 0];
			pC = (saturate((pC + 2) * 0.25f) - 0.5f) * 2 * scaleC;
			o.vertex = mul(UNITY_MATRIX_VP, float4(bp + pC, 1));
			outStream.Append(o);
			pC = posCb + sizeC * idxC[idC1 * 4 + 1];
			pC = (saturate((pC + 2) * 0.25f) - 0.5f) * 2 * scaleC;
			o.vertex = mul(UNITY_MATRIX_VP, float4(bp + pC, 1));
			outStream.Append(o);
			pC = posCb + sizeC * idxC[idC1 * 4 + 2];
			pC = (saturate((pC + 2) * 0.25f) - 0.5f) * 2 * scaleC;
			o.vertex = mul(UNITY_MATRIX_VP, float4(bp + pC, 1));
			outStream.Append(o);
			pC = posCb + sizeC * idxC[idC1 * 4 + 3];
			pC = (saturate((pC + 2) * 0.25f) - 0.5f) * 2 * scaleC;
			o.vertex = mul(UNITY_MATRIX_VP, float4(bp + pC, 1));
			outStream.Append(o);
			outStream.RestartStrip();
		}

		// rail
		const uint lenR = RAIL_INSTANCE_END - RAIL_INSTANCE_START;
		float3 aR[lenR];
		float3 bR[lenR];
		float3 cR[lenR];
		float3 dR[lenR];
		float3 wR[lenR];
		for (i = 0; i < lenR; i++)
		{
			aR[i] = bR[i] = cR[i] = dR[i] = wR[i] = 0;
		}
		for (i = 0; i < lenR; i++)
		{
			aR[i] = po[i % lenR + RAIL_INSTANCE_START];
		}
		for (i = 0; i < lenR; i++) {
			cR[i] = 3.0 * (aR[(i + lenR - 1) % lenR] - 2.0 * aR[i] + aR[(i + 1) % lenR]);
		}
		for (i = 0; i < lenR; i++) {
			float tdR = 1.0 / (4.0 - wR[(i + lenR - 1) % lenR]);
			cR[i] = (cR[i] - cR[(i + lenR - 1) % lenR]) * tdR;
			wR[i] = 1.0 * tdR;
		}
		for (i = lenR - 1; i >= 1; i--) {
			cR[i] = cR[i] - cR[(i + 1) % lenR] * wR[i];
		}
		for (i = 0; i < lenR; i++) {
			dR[i] = (cR[(i + 1) % lenR] - cR[i]) / 3.0;
			bR[i] = aR[(i + 1) % lenR] - aR[i] - cR[i] - dR[i];
		}

		const uint divR = 32;
		uint iR = id % divR;
		uint lR = id / divR % lenR;
		uint ibR = id / divR / lenR % 2;

		float eR0 = (float)(iR + 0) / divR;
		float eR1 = (float)(iR + 1) / divR;

		float3 pR0 = mad(mad(mad(dR[lR], eR0, cR[lR]), eR0, bR[lR]), eR0, aR[lR]);
		float3 pR1 = mad(mad(mad(dR[lR], eR1, cR[lR]), eR1, bR[lR]), eR1, aR[lR]);

		float3 dfR = pR1 - pR0;
		float3 vR = normalize(float3(dfR.z, 0, -dfR.x)) * 0.01f;
		vR *= ibR ? -1 : +1;
		if (id < divR * lenR * 2)
		{
			o.color0 = float3(0.4f, 0.4f, 0.4f);
			o.color1 = float3(0.6f, 0.6f, 0.6f);
			o.normal = float3(0, 1, 0);

			o.vertex = mul(UNITY_MATRIX_VP, float4(pR0 - vR, 1));
			outStream.Append(o);
			o.vertex = mul(UNITY_MATRIX_VP, float4(pR1 - vR, 1));
			outStream.Append(o);
			o.vertex = mul(UNITY_MATRIX_VP, float4(pR0 + vR, 1));
			outStream.Append(o);
			o.vertex = mul(UNITY_MATRIX_VP, float4(pR1 + vR, 1));
			outStream.Append(o);
			outStream.RestartStrip();
		}

		// train
		const float3 sizeT = float3(0.02f, 0.02f, 0.125f);
		uint idT0 = id / 6;
		uint idT1 = id % 6;
		float tT0 = (_Time.y * 0.25f) + 0.00f;
		float tT1 = (_Time.y * 0.25f) + 0.01f;
		float eT0 = frac(tT0);
		float eT1 = frac(tT1);
		uint iT0 = (uint)tT0 % lenR;
		uint iT1 = (uint)tT1 % lenR;
		float3 pT0 = mad(mad(mad(dR[iT0], eT0, cR[iT0]), eT0, bR[iT0]), eT0, aR[iT0] + float3(0, 0.02f, 0)); 
		float3 pT1 = mad(mad(mad(dR[iT1], eT1, cR[iT1]), eT1, bR[iT1]), eT1, aR[iT1] + float3(0, 0.02f, 0)); 
		float3 pT = (pT0 + pT1) * 0.5f;
		float3 dT = pT1 - pT0;
		float aTy = atan2(-dT.z, dT.x);
		float aTx = atan2(-dT.y, length(dT.xz));
		if (id < 6)
		{
			float3x3 matTy = MatRollY(aTy);
			float3x3 matTx = MatRollX(aTx);
			float3x3 matT = mul(matTy, matTx);
			float3 vT0 = mul(matT, sizeT * idxC[idT1 * 4 + 0]);
			float3 vT1 = mul(matT, sizeT * idxC[idT1 * 4 + 1]);
			float3 vT2 = mul(matT, sizeT * idxC[idT1 * 4 + 2]);
			float3 vT3 = mul(matT, sizeT * idxC[idT1 * 4 + 3]);
			
			o.color0 = float3(0.4f, 0.6f, 0.8f);
			o.color1 = o.color0;
			o.normal = normalize(cross(vT0 - vT1, vT2 - vT1));
			float3 pC;
			pC = pT + vT0;
			o.vertex = mul(UNITY_MATRIX_VP, float4(pC, 1));
			outStream.Append(o);
			pC = pT + vT2;
			o.vertex = mul(UNITY_MATRIX_VP, float4(pC, 1));
			outStream.Append(o);
			pC = pT + vT1;
			o.vertex = mul(UNITY_MATRIX_VP, float4(pC, 1));
			outStream.Append(o);
			pC = pT + vT3;
			o.vertex = mul(UNITY_MATRIX_VP, float4(pC, 1));
			outStream.Append(o);
			outStream.RestartStrip();
		}

#endif

	}

	float4 mainFS(GS_OUT i) : SV_Target
	{
		if (!any(i.normal))
		{
			return float4(i.color0, 1);
		}

		return float4(lerp(i.color0, i.color1, smoothstep(0.925f, 0.975f, i.normal.y)) * saturate(saturate(dot(i.normal, i.light) + 0.1f) + 0.1f), 1);
	}
		ENDCG
	}
	}
}
