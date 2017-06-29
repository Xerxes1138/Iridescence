/*Iridescence shader made by @xerxes1138, based on the code from : https://belcour.github.io/blog/research/2017/05/01/brdf-thin-film.html
//
//
//
//"A Practical Extension to Microfacet Theory for the Modeling of Varying Iridescence
//Laurent Belcour, Pascal Barla
//ACM Transactions on Graphics (proc. of SIGGRAPH 2017)
//
//May 2017"
//
*/
Shader "Xerxes1138/Iridescence"
{
	Properties
	{
		_Dinc("Dinc", Range(0.0, 10.0)) = 0.570
		_eta2("eta2", Range(1.0, 5.0)) = 1.8
		_eta3("eta3", Range(1.0, 5.0)) = 1.08
		_kappa3("kappa3", Range(0.0, 5.0)) = 0.51
		_alpha("alpha", Range(0.01, 1.0)) = 0.07
		//_IBLTex ("IBL", Cube) = "black" {} // Used to test IS reference, 64 to 128 samples are enough with filtered importance sampling, cube face size is hardcoded at 128px
	}
	SubShader
	{
		Tags { "RenderType"="Opaque" }
		LOD 100

		Pass // One pass forward shader
		{
			Tags {"LightMode"="ForwardBase"}

			CGPROGRAM
			#pragma target 5.0
			#pragma vertex vert
			#pragma fragment frag

			#pragma multi_compile_fwdbase nolightmap nodirlightmap nodynlightmap novertexlight
			#pragma multi_compile_fog
			
			#include "UnityCG.cginc"
			#include "Lighting.cginc"
			#include "AutoLight.cginc"

			// Shader parameters
			samplerCUBE _IBLTex;

			half4		_IBLTex_HDR;

			half		_Dinc,
						_eta2,
						_eta3,
						_kappa3,
						_alpha;

			// Common constants
			#define PI 3.14159265358979323846
		
			// XYZ to CIE 1931 RGB color space (using neutral E illuminant)
			#define XYZ_TO_RGB float3x3(2.3706743, -0.5138850, 0.0052982,-0.9000405, 1.4253036, -0.0146949, -0.4706338, 0.0885814, 1.0093968)

			// Square functions for cleaner code
			float sqr(float x) {return x*x;}
			float2 sqr(float2 x) {return x*x;}

			// Depolarization functions for natural light
			float depol (float2 polV){ return 0.5 * (polV.x + polV.y); }
			float3 depolColor (float3 colS, float3 colP){ return 0.5 * (colS + colP); }

			// GGX distribution function
			float GGX(float NdotH, float a) 
			{
				float a2 = sqr(a);
				float d = sqr(sqr(NdotH) * (a2 - 1.0) + 1.0);
				return a2 / (PI * d + 1e-7);
			}

			// Smith GGX geometric functions
			float smithG1_GGX(float NdotV, float a) 
			{
				float a2 = sqr(a);
				return 2.0 / (1.0 + sqrt(1.0 + a2 * (1.0-sqr(NdotV)) / sqr(NdotV) ));
			}

			float smithG_GGX(float NdotL, float NdotV, float a) 
			{
				return smithG1_GGX(NdotL, a) * smithG1_GGX(NdotV, a);
			}

			// Bruce Walter, Stephen R. Marschner, Hongsong Li, and Kenneth E. Torrance. Microfacet models forrefraction through rough surfaces. In Proceedings of the 18th Eurographics conference on RenderingTechniques, EGSR'07
			float G_GGX(float Roughness, float NdotL, float NdotV)
			{
				float m = Roughness;
				float m2 = m * m;

				float G_L = 1.0f / (NdotL + sqrt(m2 + (1 - m2) * NdotL * NdotL));
				float G_V = 1.0f / (NdotV + sqrt(m2 + (1 - m2) * NdotV * NdotV));
				float G = G_L * G_V;
	
				return G;
			}

			// Fresnel equations for dielectric/dielectric interfaces.
			void fresnelDielectric(in float ct1, in float n1, in float n2, out float2 R, out float2 phi) 
			{
				float st1  = (1 - ct1 * ct1); // Sinus theta1 'squared'
				float nr  = n1 / n2;

				if(sqr(nr) * st1 > 1) // Total reflection
				{
					float2 R = 1.0;

					float2 var = float2(-sqr(nr) *  sqrt(st1 - 1.0 / sqr(nr)) / ct1, -sqrt(st1 - 1.0 / sqr(nr)) / ct1);

					phi = 2.0 * atan(var);
				} 
				else // Transmission & Reflection
				{
					float ct2 = sqrt(1 - sqr(nr) * st1);

					float2 r = float2
					(
						(n2 * ct1 - n1 * ct2) / (n2 * ct1 + n1 * ct2),
						(n1 * ct1 - n2 * ct2) / (n1 * ct1 + n2 * ct2)
					);

					phi.x = (r.x < 0.0) ? PI : 0.0;
					phi.y = (r.y < 0.0) ? PI : 0.0;

					R = sqr(r);
				}
			}

			// Fresnel equations for dielectric/conductor interfaces.
			void fresnelConductor(in float ct1, in float n1, in float n2, in float k, out float2 R, out float2 phi) 
			{
				if (k == 0)// use dielectric formula to avoid numerical issues
				{ 
					fresnelDielectric(ct1, n1, n2, R, phi);
					return;
				}

				float A = sqr(n2) * (1.0 - sqr(k)) - sqr(n1) * (1.0 - sqr(ct1));
				float B = sqrt(sqr(A) + sqr(2.0 * sqr(n2) * k));
				float U = sqrt((A + B) / 2.0);
				float V = sqrt((B - A) / 2.0);

				R.y = (sqr(n1 * ct1 - U) + sqr(V)) / (sqr(n1 * ct1 + U) + sqr(V));

				float2 var1 = float2(2.0 * n1 * V * ct1, sqr(U) + sqr(V) - sqr(n1* ct1));
				phi.y = atan2(var1.x, var1.y) + PI;

				R.x = ( sqr(sqr(n2)*(1-sqr(k))*ct1 - n1*U) + sqr(2*sqr(n2)*k*ct1 - n1*V) ) / ( sqr(sqr(n2)*(1-sqr(k))*ct1 + n1*U) + sqr(2*sqr(n2)*k*ct1 + n1*V) );

				float2 var2 = float2(2*n1*sqr(n2)*ct1 * (2*k*U - (1-sqr(k))*V), sqr(sqr(n2)*(1+sqr(k))*ct1) - sqr(n1)*(sqr(U)+sqr(V))) ;
				phi.x = atan2(var2.x, var2.y);
			}

			// Evaluation XYZ sensitivity curves in Fourier space
			float3 evalSensitivity(float opd, float shift)
			{
				// Use Gaussian fits, given by 3 parameters: val, pos and var
				float phase = 2.0 * PI * opd * 1e-6;
				float3 val = float3(5.4856e-13, 4.4201e-13, 5.2481e-13);
				float3 pos = float3(1.6810e+06, 1.7953e+06, 2.2084e+06);
				float3 var = float3(4.3278e+09, 9.3046e+09, 6.6121e+09);
				float3 xyz = val * sqrt(2.0 * PI * var) * cos(pos * phase + shift) * exp(-var * phase * phase);
				xyz.x += 9.7470e-14 * sqrt(2.0 * PI * 4.5282e+09) * cos(2.2399e+06 * phase + shift) * exp(-4.5282e+09 * phase * phase);
				return xyz / 1.0685e-7;
			}

			// Main function expected by BRDF Explorer
			float3 BRDF(float3 L, float3 V, float3 H, float3 N) 
			{
				float Dinc = _Dinc;
				float eta2 = max(_eta2, 1.000277);
				float eta3 = max(_eta3, 1.000277);
				float kappa3 = max(_kappa3, 1e-5);
				float alpha = max(_alpha, 0.05);

				// Force eta_2 -> 1.0 when Dinc -> 0.0
				float eta_2 = lerp(1.0, eta2, smoothstep(0.0, 0.03, Dinc));

				// Compute dot products
				float NdotL = /*saturate*/(dot(N,L));
				float NdotV = /*saturate*/(dot(N,V));
				if (NdotL < 0 || NdotV < 0) return 0.0;
				//H = normalize(L + V);
				float NdotH = /*saturate*/(dot(N,H));
				float VdotH = /*saturate*/(dot(V,H));
				float cosTheta1 = /*saturate*/(dot(H,L));
				float cosTheta2 = sqrt(1.0 - sqr(1.0 / eta_2) * (1.0 - sqr(cosTheta1)));

				// First interface
				float2 R12, phi12;
				fresnelDielectric(cosTheta1, 1.0, eta_2, R12, phi12);
				float2 R21 = R12;
				float2 T121 = 1.0 - R12;
				float2 phi21 = PI - phi12;

				// Second interface
				float2 R23, phi23;
				fresnelConductor(cosTheta2, eta_2, eta3, kappa3, R23, phi23);

				// Phase shift
				float OPD = Dinc * cosTheta2;
				float2 phi2 = phi21 + phi23;

				// Compound terms
				float3 I = 0.0;
				float2 R123 = R12 * R23;
				float2 r123 = sqrt(R123);
				float2 Rs   = sqr(T121)*R23 / (1-R123);

				// Reflectance term for m=0 (DC term amplitude)
				float2 C0 = R12 + Rs;
				float3 S0 = evalSensitivity(0.0, 0.0);
				I += depol(C0) * S0;

				// Reflectance term for m>0 (pairs of diracs)
				float2 Cm = Rs - T121;
				for (int m=1; m<=3; ++m)
				{
					Cm *= r123;
					float3 SmS = 2.0 * evalSensitivity(m*OPD, m*phi2.x);
					float3 SmP = 2.0 * evalSensitivity(m*OPD, m*phi2.y);
					I += depolColor(Cm.x*SmS, Cm.y*SmP);
				}

				// Convert back to RGB reflectance
				I = saturate(mul(I, XYZ_TO_RGB)); 

				// Microfacet BRDF formula
				float D = GGX(NdotH, alpha);
				float G = smithG_GGX(NdotL, NdotV, alpha);
				return (D*G*I) / (4.0 * NdotL * NdotV);
			}

			float4 TangentToWorld(float3 N, float4 H)
			{
				float3 UpVector = abs(N.z) < 0.999 ? float3(0,0,1) : float3(1,0,0);
				float3 T = normalize( cross( UpVector, N ) );
				float3 B = cross( N, T );
			 
				return float4((T * H.x) + (B * H.y) + (N * H.z), H.w);
			}

			float calcLOD(int cubeSize, float pdf, int NumSamples)
			{
				float lod = (0.5 * log2( (cubeSize*cubeSize)/float(NumSamples) ) + 2.0) - 0.5*log2(pdf); 
				return lod;
			}

			// Brian Karis, Epic Games "Real Shading in Unreal Engine 4"
			float4 ImportanceSampleGGX(float2 Xi, float Roughness)
			{
				float m = Roughness;
				float m2 = m * m;
		
				float Phi = 2 * PI * Xi.x;
				 
				float CosTheta = sqrt((1.0 - Xi.y) / (1.0 + (m2 - 1.0) * Xi.y));
				float SinTheta = sqrt(max(1e-5, 1.0 - CosTheta * CosTheta));
				 
				float3 H;
				H.x = SinTheta * cos(Phi);
				H.y = SinTheta * sin(Phi);
				H.z = CosTheta;
		
				float d = (CosTheta * m2 - CosTheta) * CosTheta + 1;
				float D = m2 / (PI * d * d);
				float pdf = D * CosTheta;

				return float4(H, pdf); 
			}

			uint ReverseBits32(uint bits)
			{
				#if 0 // Shader model 5
					return reversebits(bits);
				#else
					bits = ( bits << 16) | ( bits >> 16);
					bits = ((bits & 0x00ff00ff) << 8) | ((bits & 0xff00ff00) >> 8);
					bits = ((bits & 0x0f0f0f0f) << 4) | ((bits & 0xf0f0f0f0) >> 4);
					bits = ((bits & 0x33333333) << 2) | ((bits & 0xcccccccc) >> 2);
					bits = ((bits & 0x55555555) << 1) | ((bits & 0xaaaaaaaa) >> 1);
					return bits;
				#endif
			}

			float RadicalInverse_VdC(uint bits) 
			{
				return float(ReverseBits32(bits)) * 2.3283064365386963e-10; // 0x100000000
			}

			float2 Hammersley2d(uint i, uint maxSampleCount)
			{
				return float2(float(i) / float(maxSampleCount), RadicalInverse_VdC(i));
			}

			float3 Integrate_GGXIridescence(float Roughness, float3 N, float3 V, float3 R, uint NumSamples, int cubeSize )
			{
				float3 SpecularLighting = 0.0;
				#if 1 // No IS
					Unity_GlossyEnvironmentData IBLData;
					IBLData.roughness = Roughness;
					IBLData.reflUVW = R;

					float3 SampleColor = Unity_GlossyEnvironment (UNITY_PASS_TEXCUBE(unity_SpecCube0), unity_SpecCube0_HDR, IBLData);

					float Dinc = _Dinc;
					float eta2 = max(_eta2, 1.000277);
					float eta3 = max(_eta3, 1.000277);
					float kappa3 = max(_kappa3, 1e-5);
					float alpha = max(_alpha, 0.05);

					// Force eta_2 -> 1.0 when Dinc -> 0.0
					float eta_2 = lerp(1.0, eta2, smoothstep(0.0, 0.03, Dinc));

					float cosTheta1 = dot(N, V);

					float cosTheta2 = sqrt(1.0 - sqr(1.0 / eta_2) * (1.0 - sqr(cosTheta1)));

					// First interface
					float2 R12, phi12;
					fresnelDielectric(cosTheta1, 1.0, eta_2, R12, phi12);
					float2 R21 = R12;
					float2 T121 = 1.0 - R12;
					float2 phi21 = PI - phi12;

					// Second interface
					float2 R23, phi23;
					fresnelConductor(cosTheta2, eta_2, eta3, kappa3, R23, phi23);

					// Phase shift
					float OPD = Dinc * cosTheta2;
					float2 phi2 = phi21 + phi23;

					// Compound terms
					float3 I = 0.0;
					float2 R123 = R12 * R23;
					float2 r123 = sqrt(R123);
					float2 Rs   = sqr(T121)*R23 / (1-R123);

					// Reflectance term for m=0 (DC term amplitude)
					float2 C0 = R12 + Rs;
					float3 S0 = evalSensitivity(0.0, 0.0);
					I += depol(C0) * S0;

					// Reflectance term for m>0 (pairs of diracs)
					float2 Cm = Rs - T121;
					for (int m = 1; m <= 3; ++m)
					{
						Cm *= r123;
						float3 SmS = 2.0 * evalSensitivity(m * OPD, m * phi2.x);
						float3 SmP = 2.0 * evalSensitivity(m * OPD, m * phi2.y);
						I += depolColor(Cm.x * SmS, Cm.y * SmP);
					}

					// Convert back to RGB reflectance
					I = saturate(mul(I, XYZ_TO_RGB)); 

					SpecularLighting = SampleColor * I;
				#else // IS
					for( uint i = 0; i < NumSamples; i++ )
					{
						float2 Xi = Hammersley2d( i, NumSamples );
			
						float4 H = TangentToWorld(N, ImportanceSampleGGX(Xi, Roughness));

						float3 L = 2.0 * dot( V, H ) * H - V;

						float NoV = saturate( dot( N, V ) );      
	 					float NoL = saturate( dot( N, L ) );
						float NoH = saturate( dot( N, H ) );
						float VoH = saturate( dot( V, H ) );

						if( NoL > 0 )
						{            
							float3 SampleColor = DecodeHDR(texCUBElod(_IBLTex, float4(L, calcLOD(cubeSize, H.w, NumSamples))), _IBLTex_HDR);				
	
							float Dinc = _Dinc;
							float eta2 = max(_eta2, 1.000277);
							float eta3 = max(_eta3, 1.000277);
							float kappa3 = max(_kappa3, 1e-5);
							float alpha = max(_alpha, 0.05);

							// Force eta_2 -> 1.0 when Dinc -> 0.0
							float eta_2 = lerp(1.0, eta2, smoothstep(0.0, 0.03, Dinc));

							// Compute dot products
							float NdotL = dot(N, L);
							float NdotV = dot(N, V);

							//if (NdotL < 0 || NdotV < 0) return 0.0;

							float NdotH = dot(N, H);
							float VdotH = dot(V, H);
	
							float cosTheta1 = dot(H, L);		

							float cosTheta2 = sqrt(1.0 - sqr(1.0 / eta_2) * (1.0 - sqr(cosTheta1)));

							// First interface
							float2 R12, phi12;
							fresnelDielectric(cosTheta1, 1.0, eta_2, R12, phi12);
							float2 R21 = R12;
							float2 T121 = 1.0 - R12;
							float2 phi21 = PI - phi12;

							// Second interface
							float2 R23, phi23;
							fresnelConductor(cosTheta2, eta_2, eta3, kappa3, R23, phi23);

							// Phase shift
							float OPD = Dinc * cosTheta2;
							float2 phi2 = phi21 + phi23;

							// Compound terms
							float3 I = 0.0;
							float2 R123 = R12 * R23;
							float2 r123 = sqrt(R123);
							float2 Rs   = sqr(T121)*R23 / (1-R123);

							// Reflectance term for m=0 (DC term amplitude)
							float2 C0 = R12 + Rs;
							float3 S0 = evalSensitivity(0.0, 0.0);
							I += depol(C0) * S0;

							// Reflectance term for m>0 (pairs of diracs)
							float2 Cm = Rs - T121;
							for (int m = 1; m <= 3; ++m)
							{
								Cm *= r123;
								float3 SmS = 2.0 * evalSensitivity(m * OPD, m * phi2.x);
								float3 SmP = 2.0 * evalSensitivity(m * OPD, m * phi2.y);
								I += depolColor(Cm.x * SmS, Cm.y * SmP);
							}

							// Convert back to RGB reflectance
							I = clamp(mul(I, XYZ_TO_RGB), 0.0, 1.0); 

							float G = G_GGX(alpha, NdotL, NdotV);

							SpecularLighting += SampleColor * I * 4.0 * G * NdotL * VdotH / NdotH;	
						}  
					}
					SpecularLighting = SpecularLighting / NumSamples;
				#endif

				return SpecularLighting;
			}

			struct VertexInput
			{
				float4 vertex : POSITION;
				float3 normal : NORMAL;
			};

			struct VertexOutput
			{
				float4 vertex : SV_POSITION;
				float3 worldPos : TEXCOORD0;
				float3 normal : TEXCOORD1;
			};

			VertexOutput vert (VertexInput v)
			{
				VertexOutput o;
				o.vertex = UnityObjectToClipPos(v.vertex);
				o.worldPos = mul(unity_ObjectToWorld, v.vertex);
				o.normal = UnityObjectToWorldNormal(v.normal.xyz);
				UNITY_TRANSFER_FOG(o,o.vertex);
				return o;
			}

			float4 frag (VertexOutput i) : SV_Target
			{
				float3 L = _WorldSpaceLightPos0.xyz;
				float3 V = normalize(i.worldPos - _WorldSpaceCameraPos.xyz);
				float3 H = normalize(L + (-V));
				float3 N = normalize(i.normal);
				float3 R = reflect(V, N);

                float3 IBL = Integrate_GGXIridescence(_alpha, N, -V, R, 64 /*NumSamples*/, 128 /*cubeface size*/);

				float3 brdf = IBL + BRDF(L, -V, H, N) * saturate(dot(N, L)) * _LightColor0.rgb;

				return float4(brdf, 1.0);
			}
			ENDCG
		}
	}
}
