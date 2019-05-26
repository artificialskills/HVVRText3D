/**
* Copyright (c) 2018-present, Apollo Ellis
* All rights reserved.
*
* This source code is licensed under the BSD-style license found in the
* LICENSE file in the root directory of this source tree. An additional grant
* of patent rights can be found in the PATENTS file in the same directory.
*/


#pragma once
#include <string>
#include <float.h>

#include "vector_math.h"
#include "gpu_buffer.h"
#include "thorax/grid_types.h"

namespace hvvr{
struct TextGridGPU {
	int hres, vres;
	float lowest_x;
	float lowest_y;
	float highest_x;
	float highest_y;
	Grid::TextGridCell *cells;
	Grid::GlyphPtr *glyph_ptrs;
	Grid::GlyphGrid *glyph_grids;
	Grid::GlyphGridCell *glyph_grid_cells;
	Grid::Shape *shapes;
	Grid::shape_ptr *shape_ptrs;
};

void UploadGrid(Grid::TextGrid *tgrid);
}