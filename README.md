# HVVRText3D

HVVRText3D integrates fast high quality analytic text rendering into HVVR,
leveraging capabilities in Thorax True-Type.

It is a purely analytic anti-aliasing technique, which computes the 
coverage of vector primitives, found in True-Type fonts (pre-processed 
by Thorax), over pixel footprints directly.

The techniques we employ in the math make the computation extremely rapid.
The crux of the novel implementation is in:
libraries/hvvr/raycaster/grid_sampler.h

The individual projects that are integrated, HVVR and Thorax, have their 
respective readme.md text inlined below.

The HVVR project was forked by myself in July 2018, and the fork was 
converted into this upload. Do not expect it to adhere to any expectations
of, or match, HVVR proper. Of course licensing has been preserved. See 
NOTICE.md and LICENSE in the root directory of this project.

Thorax was downloaded from github by myself in July 2018. Of course licensing has 
been preserved. See NOTICE.md and LICENSE2.

My own code leveraging HVVR and Thorax e.g. thorax_bridge,
grid, grid_types, grid_builders, and grid_sampler, you will find have been
licensed under BSD Three Clause, the same license as HVVR, which you can
again find in LICENSE.

Most of the usage instructions for HVVR apply to HVVRText3D with the exception,
you should set textviewer as the startup project instead of modelviewer.

# Thorax Truetype

Thorax Truetype is a library for loading and rendering TrueType fonts.  This
project grew out of a pedagogical exercise to have an easy to understand
TrueType loader.  My desire to do so came from stb_truetype being difficult
to understand due to a lack of C++ features or type punning (the address
computations required to index into the truetype tables are pretty
incomprehensible as inline C code) and the size of the opentype loader.

The loader consists of two files:
  thorax_truetype.h
  thorax_truetype.cpp

In addition to the font loader, this project contains a font rasterizer in the
renderer folder.  The rasterizer uses a novel coverage algorithm and is capable
of generated nicely anti-aliased renderings.

The files in the renderer:
  build.h         : interface for managing text scenes
  bvh.cpp         : implementation of bounding volume hierarchy builders
  drawcontext.h   : context interface for drawing text to a bitmap
  drawcontext.cpp : context implementation for drawing text to a bitmap
  types.h         : a variety of helper datatypes

A more detailed description of the rasterizer will come later.


# HVVR
HVVR (Hierarchical Visibility for Virtual Reality) is an optimized software raycaster. It implements a hybrid CPU/GPU raycaster, suited for real-time rendering of effects such as lens distortion.

## Examples
modelviewer.vcxproj (vs2015/hvvr.sln) implements a simple model viewer, with optional #defines to enable depth of field and foveated rendering.

## Requirements
HVVR requires:
* Windows
* CUDA-compatible GPU and CUDA SDK (tested on GeForce GTX 1080 and Quadro P6000)
* CPU with AVX2 support
* (optional) FBX SDK for importing content

## Building HVVR
0. Install CUDA SDK (and optionally FBX SDK)
1. Open vs2015/hvvr.sln
2. Select the "Develop/x64" configuration
3. Set "modelviewer" as the startup project
4. Build solution

By default, raycaster.vcxproj is configured to build all .cu files to the compute_61,sm_61 target. If your GPU needs a different code gen target, you can change this in the raycaster project properties.

If enabling the FBX import path, you'll need to install the FBX SDK, and you'll need to add it to your Microsoft.Cpp.x64.user "VC++ Directories" paths. For example, append "C:\Program Files\Autodesk\FBX\FBX SDK\2018.1.1\include;" to your "Include Directories" and "C:\Program Files\Autodesk\FBX\FBX SDK\2018.1.1\lib\vs2015\x64\;" to your "Library Directories".

## Running HVVR
1. Download pre-built scene data files from: https://github.com/facebookresearch/HVVR/releases/latest
2. Copy pre-built scene data files (for example, sponza.bin) to HVVR/libraries/hvvr/samples_shared/data/scenes
3. Set gSceneSelect in modelviewer.cpp to match the scene data file you want to load

## Importing content
To convert and load your own scenes:
1. Enable FBX import in HVVR/libraries/hvvr/samples_shared/model_import.h by setting MODEL_IMPORT_ENABLE_FBX to 1
2. Change the value of scenePath in gOnInit() in HVVR/projects/hvvr_samples/modelviewer/modelviewer.cpp to point at your FBX file
3. (Optionally, you can also convert the FBX file to .bin for faster loading in the future, see the modelconvert project)

## Full documentation
Some installations of the CUDA SDK and Visual Studio integration appear to have broken build dependency tracking. If you're seeing behavior after code modification that looks like stale code from header files referenced in .cu files, try cleaning the solution and building from scratch (or install a newer CUDA Toolkit).

There are several #define macros to be aware of, which control the behavior of the raycaster:
* DISABLE_MOVEMENT - modelviewer.cpp, disables camera movement and rotation
* ENABLE_VSYNC - modelviewer.cpp, throttle to monitor refresh rate
* ENABLE_DEPTH_OF_FIELD - modelviewer.cpp
* ENABLE_FOVEATED - modelviewer.cpp
* GAZE_CURSOR_MODE - modelviewer.cpp, controls how the foveated gaze point is set
* RT_WIDTH and RT_HEIGHT - modelviewer.cpp, sets the render target size in pixels and window client rect size
* MODEL_IMPORT_ENABLE_FBX - model_import.h, enables FBX model loading (requires FBX SDK)
* MAX_TRI_INDICES_TO_INTERSECT - magic_constants.h, this controls the size of the intermediate buffer which contains the list of triangles to intersect after BVH traversal. You may need to increase this to avoid a crash/assert in more complex scenes.
* COLOR_MODE_MSAA_RATE - raycaster_common.h, controls the number of subsamples per ray, valid values are 1 to 32
* SUPERSHADING_MODE - kernel_constants.h, allows switching between MSAA and SSAA shading
* JITTER_SAMPLES - kernel_constants.h, enables temporal sample position jitter when combined with foveated rendering
* FOVEATED_TRIANGLE_FRUSTA_TEST_DISABLE - kernel_constants.h, bug workaround to prevent incorrect triangle culling when foveated rendering is enabled
* ENABLE_HACKY_WIDE_FOV, HACKY_WIDE_FOV_W, and HACKY_WIDE_FOV_H - kernel_constants.h, forces rendering to a specified wide FoV

Enabling foveated rendering:
1. Set ENABLE_FOVEATED to 1
2. (Optional) Set GAZE_CURSOR_MODE to GAZE_CURSOR_MODE_MOUSE
3. Set JITTER_SAMPLES to 1 to enable temporal sample jitter
4. Set FOVEATED_TRIANGLE_FRUSTA_TEST_DISABLE to 1

Enabling depth of field:
1. Set ENABLE_DEPTH_OF_FIELD to 1
2. (Optional) Adjust "lensRadius" and "focalDistance" variables in modelviewer.cpp

Enabling wide field of view: (not tested in combination with other config options, yet)
1. Set ENABLE_HACKY_WIDE_FOV to 1
2. Set HACKY_WIDE_FOV_W and HACKY_WIDE_FOV_H to desired horizontal and vertical fields of view (in degrees)

Benchmarking vs other raytracers:
1. Set DISABLE_MOVEMENT to 1 to ensure consistent output
2. Enable export of the scene and rays in libraries/hvvr/raycaster/render.cpp
3. Import the scene and rays to https://github.com/anankervis/RaytraceBenchmarking

## Join the HVVR community
See the CONTRIBUTING file for how to help out.

## License
HVVR is BSD-licensed.