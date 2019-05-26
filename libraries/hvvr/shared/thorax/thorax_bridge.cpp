/**
 * Copyright (c) 2018-present, Apollo Ellis
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree. An additional grant
 * of patent rights can be found in the PATENTS file in the same directory.
 */

#include <iostream>
#include "thorax_bridge.h"
#include "renderer/build.h"

Grid::TextGrid* GetTextGridCPU(thorax_ai::Array<thorax_ai::Shape> &shapes, thorax_ai::Array<Grid::GlyphGrid> &glyphGrids, thorax_ai::Array<Grid::GlyphGridCell> &glyphGridCells,
	thorax_ai::Array<Grid::shape_ptr> &shapePtrs){
	
	Grid::Shape *cpuShapes = new Grid::Shape[shapes.size];
	for (int i = 0; i < shapes.size; i++) {
		cpuShapes[i].FromShape(shapes[i].p0v1p2v3.data.m256_f32);
	}
	Grid::TextGrid *textGrid = new Grid::TextGrid;
	textGrid->shapes = cpuShapes;
	textGrid->glyph_grids = glyphGrids.data;
	textGrid->glyph_grid_cells = glyphGridCells.data;
	textGrid->shape_ptrs = shapePtrs.data;
	textGrid->glyph_grid_count = (int)glyphGrids.size;
	textGrid->glyph_grid_cell_count = (int)glyphGridCells.size;
	textGrid->shape_count = (int)shapes.size;
	textGrid->shape_ptr_count = (int)shapePtrs.size;
	return textGrid;
}

Grid::TextGrid* InitializeText(std::string path, std::string filename)
{
	std::string textfilename = "Times\ New\ Roman.ttf";
	const Font *font;
	{
		FILE *f = new FILE();
		fopen_s(&f, (path + textfilename).c_str(), "rb");

		fseek(f, 0, SEEK_END);
		int size = ftell(f);
		rewind(f);
		const char* data = new char[size];
		fread_s((void*)data, size, sizeof(char), size, f);
		font = Font::AsFont(data);
		fclose(f);
	}
	thorax_ai::FontRenderInfo *ft_render_info = new thorax_ai::FontRenderInfo();
	ft_render_info->Initialize(font);

	std::string words = "To be, or not to be, that is the question :";

	unsigned points[4096];
	int hres = 0;
	int vresmax = 0;
	int hresmax = 0;
	int counts = 0;
	for (int i = 0; i < words.length(); i++)
	{
		if (words[i] == '\n') {
			vresmax++;
			hres = 0;
		}
		else {
			hres++;
			hresmax = hresmax < hres ? hres : hresmax;
		}
		points[i] = words[i];
		counts++;
	}

	int refs = (int)ft_render_info->LayoutGlyphs((Grid::GridRef*)NULL, 0, (const unsigned *)points, counts);
	Grid::GridRef *glyphRefs = new Grid::GridRef[refs];
	int numRefs = ft_render_info->LayoutGlyphs(glyphRefs, 0, (const unsigned *)points, counts);

	Grid::TextGrid *tgrid = GetTextGridCPU(ft_render_info->shapes, ft_render_info->glyph_grids,
		ft_render_info->glyph_grid_cells, ft_render_info->shape_ptrs);
	GridBuild(tgrid, glyphRefs, numRefs, hresmax, vresmax);
	return tgrid;
}