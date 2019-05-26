/**
* Copyright (c) 2018-present, Apollo Ellis
* All rights reserved.
*
* This source code is licensed under the BSD-style license found in the
* LICENSE file in the root directory of this source tree. An additional grant
* of patent rights can be found in the PATENTS file in the same directory.
*/

#include "grid_buffers.h"

namespace hvvr {

	TextGridGPU gDeviceTextGrid;
	GPUBuffer<Grid::TextGridCell> cells;
	GPUBuffer<Grid::GlyphPtr> glyph_ptrs;
	GPUBuffer<Grid::GlyphGrid> glyph_grids;
	GPUBuffer<Grid::GlyphGridCell> glyph_grid_cells;
	GPUBuffer<Grid::Shape> shapes;
	GPUBuffer<Grid::shape_ptr> shape_ptrs;
	
	void UploadGrid(Grid::TextGrid *tgrid)
	{
		gDeviceTextGrid.hres = tgrid->hres;
		gDeviceTextGrid.vres = tgrid->vres;
		gDeviceTextGrid.lowest_x = tgrid->lowest_x;
		gDeviceTextGrid.lowest_y = tgrid->lowest_y;
		gDeviceTextGrid.highest_x = tgrid->highest_x;
		gDeviceTextGrid.highest_y = tgrid->highest_y;

		cells.resizeDestructive(tgrid->cell_count);
		glyph_ptrs.resizeDestructive(tgrid->glyph_ptr_count);
		glyph_grids.resizeDestructive(tgrid->glyph_grid_count);
		glyph_grid_cells.resizeDestructive(tgrid->glyph_grid_cell_count);
		shapes.resizeDestructive(tgrid->shape_count);
		shape_ptrs.resizeDestructive(tgrid->shape_ptr_count);

		cells.upload(tgrid->cells);
		glyph_ptrs.upload(tgrid->glyph_ptrs);
		glyph_grids.upload(tgrid->glyph_grids);
		glyph_grid_cells.upload(tgrid->glyph_grid_cells);
		shapes.upload(tgrid->shapes);
		shape_ptrs.upload(tgrid->shape_ptrs);
		
		gDeviceTextGrid.cells = cells.data();
		gDeviceTextGrid.glyph_ptrs = glyph_ptrs.data();
		gDeviceTextGrid.glyph_grids = glyph_grids.data();
		gDeviceTextGrid.glyph_grid_cells = glyph_grid_cells.data();
		gDeviceTextGrid.shapes = shapes.data();
		gDeviceTextGrid.shape_ptrs = shape_ptrs.data();
	}
}