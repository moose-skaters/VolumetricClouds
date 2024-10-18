//
// Function to remap a value from one range to another. It is slightly cheaper than SetRange
//
#define ValueRemapFuncionDef(DATA_TYPE) \
	DATA_TYPE ValueRemap(DATA_TYPE inValue, DATA_TYPE inOldMin, DATA_TYPE inOldMax, DATA_TYPE inMin, DATA_TYPE inMax) \
	{ \
		DATA_TYPE old_min_max_range = (inOldMax - inOldMin); \
		DATA_TYPE clamped_normalized = saturate((inValue - inOldMin) / old_min_max_range); \
		return inMin + (clamped_normalized*(inMax - inMin)); \
	}

ValueRemapFuncionDef(float)
ValueRemapFuncionDef(float2)
ValueRemapFuncionDef(float3)
ValueRemapFuncionDef(float4)

//
// Function to erode a value given an erosion amount. A simplified version of SetRange.
//
float ValueErosion(float inValue, float inOldMin)
{
	// derrived from Set-Range, this function uses the oldMin to erode or inflate the input value. - inValues inflate while + inValues erode
	float old_min_max_range = (1.0 - inOldMin);
	float clamped_normalized = saturate((inValue - inOldMin) / old_min_max_range);
	return (clamped_normalized);
}

//
// Funcion to get a fractional height value within the cloud rendering domain. 
//
float GetFractionFromValue(float inValue, float inMin, float inMax)
{
	return saturate((inValue - inMin) / (inMax - inMin));
}

//
// Function to get mipmap level for voxel clouds
// 
float GetVoxelCloudMipLevel(CloudRenderingRaymarchInfo inRaymarchInfo, float inBaseMipLevel)
{
	// Apply Distance based Mip Offset
	float mipmap_level = cUseVoxelFineDetailMipMaps ? log2(1.0 + abs(inRaymarchInfo.mDistance * cVoxelFineDetailMipMapDistanceScale)) + inBaseMipLevel : inBaseMipLevel;

	// return result
	return mipmap_level;
}

//
// Func to upres low res voxel data
//
float GetUprezzedVoxelCloudDensity(CloudRenderingRaymarchInfo inRaymarchInfo, float3 inSamplePosition, float inDimensionalProfile, float inType, float inDensityScale, float inMipLevel, bool inHFDetails)
{
	// Apply wind offset 
	inSamplePosition -= float3(cCloudWindOffset.x, cCloudWindOffset.y, 0.0) * low_cloud_animation_speed;
	
	// Sample noise
	float mipmap_level = GetVoxelCloudMipLevel(inRaymarchInfo, inMipLevel);
	float4 noise = Cloud3DNoiseTextureC.SampleLOD(Cloud3DNoiseSamplerC, inSamplePosition * 0.01 , mipmap_level);

	// Define wispy noise
	float wispy_noise = lerp(noise.r, noise.g, inDimensionalProfile);

	// Define billowy noise 
	float billowy_type_gradient = pow(inDimensionalProfile, 0.25);
	float billowy_noise = lerp(noise.b * 0.3, noise.a * 0.3, billowy_type_gradient);

	// Define Noise composite - blend to wispy as the density scale decreases.
	float noise_composite = lerp(wispy_noise, billowy_noise, inType);

	// Get the hf noise which is to be applied nearby - First, get the distance from the sample to camera and only do the work within a distance of 150 meters. 
	if (inHFDetails)
	{
		// Get the hf noise by folding the highest frequency billowy noise. 
		float hhf_noise = saturate(lerp(1.0 - pow(abs(abs(noise.g * 2.0 - 1.0) * 2.0 - 1.0), 4.0), pow(abs(abs(noise.a * 2.0 - 1.0) * 2.0 - 1.0), 2.0), inType));
	
		// Apply the HF nosie near camera.
		float hhf_noise_distance_range_blender = ValueRemap(inRaymarchInfo.mDistance, 50.0, 150.0, 0.9, 1.0);
		noise_composite = lerp(hhf_noise, noise_composite, hhf_noise_distance_range_blender);
	}

	// Composote Noises and use as a Value Erosion
	float uprezzed_density = ValueErosion(inDimensionalProfile, noise_composite);

	// Modify User density scale
	float powered_density_scale = pow(saturate(inDensityScale), 4.0);

	// Apply User Density Scale Data to Result
	uprezzed_density *= powered_density_scale; 
		
	// Sharpen result
	uprezzed_density = pow(uprezzed_density, lerp(0.3, 0.6, max(EPSILON, powered_density_scale)));
	if (inHFDetails)
	{
		float hhf_noise_distance_range_blender = GetFractionFromValue(inRaymarchInfo.mDistance, 50.0, 150.0);
		uprezzed_density = pow(uprezzed_density, lerp(0.5, 1.0, hhf_noise_distance_range_blender)) * lerp(0.666, 1.0, hhf_noise_distance_range_blender);
	} //can try some

	// Return result with softened edges
	return uprezzed_density;
}

//
// Get voxel cloud density
//
VoxelCloudDensitySamples GetVoxelCloudDensitySamples(CloudRenderingRaymarchInfo inRaymarchInfo, float3 inSamplePosition, float inMipLevel, bool inHFDetails)
{
	// Init Output
	VoxelCloudDensitySamples density_samples;

	// Get Sample Coordinates
	float3 sample_coord = (inSamplePosition - cVoxelCloudBoundsMin) / (cVoxelCloudBoundsMax - cVoxelCloudBoundsMin);
	float3 modeling_data = VoxelCloudModelingDataTexture.SampleLOD(VoxelDataSampler, sample_coord, inMipLevel).rgb;

	// Assign data
	float dimensional_profile = modeling_data[0];
	float type = modeling_data[1];
	float density_scale = modeling_data[2];
	
	// Upres Density Data only if dimensional profile is nonzero. When rendering the reflection dome we don't up-res as fine details will hardly 
	// be visible in the reflections anyways.
	if (dimensional_profile > 0.0)
	{
		density_samples.mProfile = dimensional_profile * density_scale;

		// Scale down density close to camera so that aloy is visible inside of clouds and we mitigate the worst animated noise artifacts.
		if(!cRenderDome)
			density_samples.mFull = GetUprezzedVoxelCloudDensity(inRaymarchInfo, inSamplePosition, dimensional_profile, type, density_scale, inMipLevel, inHFDetails) * ValueRemap(inRaymarchInfo.mDistance, 10.0, 120.0, 0.25, 1.0);
		else
			density_samples.mFull = density_samples.mProfile;
	}
	else 		
	{
		density_samples.mFull = density_samples.mProfile = 0.0;
	}

	// Return result
	return density_samples;
}
