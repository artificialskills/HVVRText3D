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

#include "grids.h"
#include "renderer/drawcontext.h"

Grid::TextGrid* GetTextGridCPU(thorax_ai::Array<thorax_ai::Shape> &shapes, thorax_ai::Array<Grid::GlyphGrid> &glyphGrids, thorax_ai::Array<Grid::GlyphGridCell> &gridCells,
	thorax_ai::Array<Grid::shape_ptr> &shapePtrs);

Grid::TextGrid* InitializeText(std::string path, std::string filename);