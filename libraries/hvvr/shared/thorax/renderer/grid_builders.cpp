/**
* Copyright (c) 2018-present, Apollo Ellis
* All rights reserved.
*
* This source code is licensed under the BSD-style license found in the
* LICENSE file in the root directory of this source tree. An additional grant
* of patent rights can be found in the PATENTS file in the same directory.
*/

#include "../grids.h"
#include "drawcontext.h"
#include <float.h>
#include <limits.h>
#include <map>
#include <vector>
#include <algorithm>

using namespace Grid;

thorax_ai::GlyphGridBuilder::GlyphGridBuilder(const thorax_ai::Shape* shapes, size_t numShapes, size_t firstShape)
	: lowest_x(FLT_MAX),lowest_y(FLT_MAX),highest_x(-FLT_MAX),highest_y(-FLT_MAX)
{
	std::vector<ShapeEvent> horizontalEvents;
	std::vector<ShapeEvent> verticalEvents;
	for (int i = 0; i < numShapes; i++) {
		Box box = shapes[i].Bound();
		horizontalEvents.push_back(ShapeEvent(i, -box.data[0], box.data[2]));
		verticalEvents.push_back(ShapeEvent(i, -box.data[1], box.data[3]));
		lowest_x = -box.data[0] < lowest_x ? -box.data[0] : lowest_x;
		lowest_y = -box.data[1] < lowest_y ? -box.data[1] : lowest_y;
		highest_x = box.data[2] > highest_x ? box.data[2] : highest_x;
		highest_y = box.data[3] > highest_y ? box.data[3] : highest_y;
	}

	std::sort(horizontalEvents.begin(), horizontalEvents.end(), EventCompare);
	std::sort(verticalEvents.begin(), verticalEvents.end(), EventCompare);

	size_t averageOverlap = 0;
	for (unsigned int i = 0; i < horizontalEvents.size(); i++) {
		float endValue = horizontalEvents[i].end;
		size_t overLap = 0;
		for (unsigned int j = i; j < horizontalEvents.size(); j++) {
			if (i == j)
				continue;
			if (endValue <= horizontalEvents[j].start)
				break;
			overLap++;
		}
		averageOverlap += overLap;
	}
	averageOverlap = 1 + averageOverlap / horizontalEvents.size();
	hres = 1 + horizontalEvents.size() / averageOverlap;
	
	averageOverlap = 0;
	for (unsigned int i = 0; i < verticalEvents.size(); i++) {
		float endValue = verticalEvents[i].end;
		size_t overLap = 0;
		for (unsigned int j = i; j < verticalEvents.size(); j++) {
			if (i == j)
				continue;
			if (endValue <= verticalEvents[j].start)
				break;
			overLap++;
		}
		averageOverlap += overLap;
	}
	averageOverlap = 1 + averageOverlap / verticalEvents.size();
	vres = 1 + verticalEvents.size() / averageOverlap;

	//create grid
	int height =  (int)vres;
	int width = (int)hres;
	cellCount = height*width;
	shapeOffsets = new GlyphGridCell[cellCount];
	totalRefs = 0;

	//count the number of shapes in each cell
	memset(shapeOffsets, 0, sizeof(GlyphGridCell)*cellCount);
	for (int n = 0;n < numShapes;n++)
	{
		Shape shape = shapes[n];

		Box box = shape.Bound();
		float x_min = -box.data[0];
		float y_min = -box.data[1];
		float x_max = box.data[2];
		float y_max = box.data[3];
		float start_x = (float)width*((x_min - lowest_x) / (highest_x - lowest_x));
		float start_y = (float)height*((y_min - lowest_y) / (highest_y - lowest_y));
		float end_x = (float)width*((x_max - lowest_x) / (highest_x - lowest_x));
		float end_y = (float)height*((y_max - lowest_y) / (highest_y - lowest_y));
		for (int i = int(start_y); i <= int(end_y) && i < height; i++)
		{
			for (int j = int(start_x); j <= int(end_x) && j < width; j++)
			{
				shapeOffsets[i*width + j].count++;
				totalRefs++;
			}
		}
	}
	//prefix sum
	int last_val = 0;
	for (int i = 0; i < cellCount; i++)
	{
		int temp = last_val;
		last_val = shapeOffsets[i].count + last_val;
		shapeOffsets[i].offset = temp;
	}
	//cells inc indexes in the packed array and increment that index at each step
	//it starts off identical to offset
	int *cells_inc = new int[cellCount];
	for (int i = 0; i < cellCount; i++)
	{
		cells_inc[i] = shapeOffsets[i].offset;
	}
	//these are the data of the grid they point to shapes in the globa shape array
	shapePtrs = new shape_ptr[totalRefs]; 
	for (int n = 0; n < numShapes; n++)
	{
		Shape shape = shapes[n];
		 
		Box box = shape.Bound();
		float x_min = -box.data[0];
		float y_min = -box.data[1];
		float x_max = box.data[2];
		float y_max = box.data[3];
		float start_x = (float)width*((x_min - lowest_x) / (highest_x - lowest_x));
		float start_y = (float)height*((y_min - lowest_y) / (highest_y - lowest_y));
		float end_x = (float)width*((x_max - lowest_x) / (highest_x - lowest_x));
		float end_y = (float)height*((y_max - lowest_y) / (highest_y - lowest_y));
		for (int i = int(start_y); i <= int(end_y) && i < height; i++)
		{
			for (int j = int(start_x); j <= int(end_x) && j < width; j++)
			{
				int cell = i*width + j;
				//put in the shape ptr an increment the placement offset for the next shape ptr for that cell
				shapePtrs[cells_inc[cell]].ptr = n + (int)firstShape;
				shapePtrs[cells_inc[cell]++].id = n;
			}
		}
	}
}


void GridBuild(TextGrid *grid, const GridRef *input, size_t num_refs, int horizontal_count, int vertical_count)
{
	grid->lowest_x = FLT_MAX;
	grid->lowest_y = FLT_MAX;
	grid->highest_x = -FLT_MAX;
	grid->highest_y = -FLT_MAX;
	float &lowest_x = grid->lowest_x;
	float &lowest_y = grid->lowest_y;
	float &highest_x = grid->highest_x;
	float &highest_y = grid->highest_y;
	for (int i = 0; i < num_refs; i++)
	{
		GridRef gridref = input[i];
		GlyphGrid glyph_grid = grid->glyph_grids[gridref.nodeIndex];

		Point lpos = invert(gridref.objectFromParent) * Point(glyph_grid.lowest_x, glyph_grid.lowest_y);
		Point hpos = invert(gridref.objectFromParent) * Point(glyph_grid.highest_x, glyph_grid.highest_y);
		lowest_x = lpos.x < lowest_x ? lpos.x : lowest_x;
		lowest_y = lpos.y < lowest_y ? lpos.y : lowest_y;
		highest_x = lpos.x > highest_x ? lpos.x : highest_x;
		highest_y = lpos.y > highest_y ? lpos.y : highest_y;
		lowest_x = hpos.x < lowest_x ? hpos.x : lowest_x;
		lowest_y = hpos.y < lowest_y ? hpos.y : lowest_y;
		highest_x = hpos.x > highest_x ? hpos.x : highest_x;
		highest_y = hpos.y > highest_y ? hpos.y : highest_y;
	}

	int width = (int)((highest_x - lowest_x) / horizontal_count);
	int height = (int)((highest_y - lowest_y) / vertical_count);

	int grid_index_count = 0;
	grid->cells = new TextGridCell[horizontal_count*vertical_count];
	memset(grid->cells, 0, sizeof(TextGridCell)*horizontal_count*vertical_count);
	for (int n = 0; n < num_refs; n++)
	{
		GridRef gridref = input[n];
		GlyphGrid glyph_grid = grid->glyph_grids[gridref.nodeIndex];

		Point lpos = invert(gridref.objectFromParent) * Point(glyph_grid.lowest_x, glyph_grid.lowest_y);
		Point hpos = invert(gridref.objectFromParent) * Point(glyph_grid.highest_x, glyph_grid.highest_y);
		lpos.x /= width;
		lpos.y /= height;
		hpos.x /= width;
		hpos.y /= height;
		for (int i = int(lpos.y); i <= int(hpos.y) && i < vertical_count; i++)
		{
			for (int j = int(lpos.x); j <= int(hpos.x) && j < horizontal_count; j++)
			{
				grid_index_count++;
				grid->cells[i*horizontal_count + j].count++;
			}
		}
	}

	//genius level implementation of prefix sum
	int last_val = 0;
	for (int i = 0; i < horizontal_count*vertical_count; i++)
	{
		int temp = last_val;
		last_val = grid->cells[i].count + last_val;
		grid->cells[i].offset = temp;
	}

	int *cells_inc = new int[horizontal_count*vertical_count];
	for (int i = 0; i < horizontal_count*vertical_count; i++)
	{
		cells_inc[i] = grid->cells[i].offset;
	}
	grid->glyph_ptrs = new GlyphPtr[grid_index_count];
	for (int n = 0; n < num_refs; n++)
	{
		GridRef gridref = input[n];
		GlyphGrid glyph_grid = grid->glyph_grids[gridref.nodeIndex];

		Point lpos = invert(gridref.objectFromParent) * Point(glyph_grid.lowest_x, glyph_grid.lowest_y);
		Point hpos = invert(gridref.objectFromParent) * Point(glyph_grid.highest_x, glyph_grid.highest_y); 
		lpos.x /= width;
		lpos.y /= height;
		hpos.x /= width;
		hpos.y /= height;

		Matrix2x3 objFrPa_i = invert(gridref.objectFromParent);

		for (int i = int(lpos.y); i <= int(hpos.y) && i < vertical_count; i++)
		{
			for (int j = int(lpos.x); j <= int(hpos.x) && j < horizontal_count; j++)
			{
				int cell = i*horizontal_count + j;
				grid->glyph_ptrs[cells_inc[cell]].ptr = gridref.nodeIndex;
				grid->glyph_ptrs[cells_inc[cell]].transform_to_local = gridref.objectFromParent;
				grid->glyph_ptrs[cells_inc[cell]++].instance = n;
			}
		}
	}
	grid->hres = horizontal_count;
	grid->vres = vertical_count;
	grid->cell_count = horizontal_count * vertical_count;
	grid->glyph_ptr_count = grid_index_count;
}