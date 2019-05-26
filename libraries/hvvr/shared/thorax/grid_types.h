/**
 * Copyright (c) 2018-present, Apollo Ellis
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree. An additional grant
 * of patent rights can be found in the PATENTS file in the same directory.
 */

#pragma once

#include "cuda_decl.h"

namespace Grid
{
	struct Box {
		float data[4];
	};

	struct Shape {
		CUDA_HOST_DEVICE
		void FromShape(float p0V1P2V3[8]){
			for(int i = 0 ; i < 8; i++)
				p0v1p2v3[i] = p0V1P2V3[i];
		}
		CUDA_HOST_DEVICE 
		bool IsTrapazoid() const { return p0v1p2v3[1] != p0v1p2v3[5]; }
		CUDA_HOST_DEVICE
		bool OutBounds() const {
			float x1 = p0v1p2v3[0];
			float y1 = p0v1p2v3[1];
			float v0X = p0v1p2v3[2];
			float v0Y = p0v1p2v3[3];
			float x2 = p0v1p2v3[4];
			float y2 = p0v1p2v3[5];
			float v2X = p0v1p2v3[6];
			float v2Y = p0v1p2v3[7];
			bool out = x1 > 1 && (x1 + v0X) > 1 && (x2 + v2X) > 1 && x2 > 1;
			out |= x1 < 0 && (x1 + v0X) < 0 && (x2 + v2X) < 0 && x2 < 0;
			out |= y1 > 1 && (y1 + v0Y) > 1 && (y2 + v2Y) > 1 && y2 > 1;
			out |= y1 < 0 && (y1 + v0Y) < 0 && (y2 + v2Y) < 0 && y2 < 0;
			return out;
		}
		float p0v1p2v3[8];
	};

	struct Vector {
		CUDA_HOST_DEVICE Vector() = default;
		CUDA_HOST_DEVICE Vector(float x, float y) : x(x), y(y) {}
		CUDA_HOST_DEVICE Vector operator+(const Vector& a) const { return Vector(x + a.x, y + a.y); }
		CUDA_HOST_DEVICE Vector operator-(const Vector& a) const { return Vector(x - a.x, y - a.y); }
		CUDA_HOST_DEVICE Vector operator*(float a) const { return Vector(x * a, y * a); }
		float x, y;
	};

	struct Point {
		CUDA_HOST_DEVICE Point() = default;
		CUDA_HOST_DEVICE Point(float x, float y) : x(x), y(y) {}
		CUDA_HOST_DEVICE explicit Point(const Vector& a) : x(a.x), y(a.y) {}
		CUDA_HOST_DEVICE Vector operator-(const Point& a) const { return Vector(x - a.x, y - a.y); }
		float x, y;
	};

	struct Matrix2x3 {
		CUDA_HOST_DEVICE Matrix2x3() = default;
		CUDA_HOST_DEVICE Matrix2x3(float x_x, float y_x, float w_x, float x_y, float y_y, float w_y)
			: x(x_x, x_y), y(y_x, y_y), w(w_x, w_y) {}
		CUDA_HOST_DEVICE Matrix2x3(const Vector& x, const Vector& y, const Point& w) : x(x), y(y), w(w) {}
		CUDA_HOST_DEVICE Vector operator*(const Vector& a) const { return Vector(x.x * a.x + y.x * a.y, x.y * a.x + y.y * a.y); }
		CUDA_HOST_DEVICE Point operator*(const Point& a) const { return Point(x.x * a.x + y.x * a.y + w.x, x.y * a.x + y.y * a.y + w.y); }
		CUDA_HOST_DEVICE Matrix2x3 operator*(const Matrix2x3& a) const { return Matrix2x3((*this) * a.x, (*this) * a.y, (*this) * a.w); }
		CUDA_HOST_DEVICE Shape operator*(const Shape& a) const {
			float result[8];
			result[0] = x.x*a.p0v1p2v3[0] + y.x*a.p0v1p2v3[1] + w.x;
			result[1] = x.y*a.p0v1p2v3[0] + y.y*a.p0v1p2v3[1] + w.y;
			result[2] = x.x*a.p0v1p2v3[2] + y.x*a.p0v1p2v3[3];
			result[3] = x.y*a.p0v1p2v3[2] + y.y*a.p0v1p2v3[3];
			result[4] = x.x*a.p0v1p2v3[4] + y.x*a.p0v1p2v3[5] + w.x;
			result[5] = x.y*a.p0v1p2v3[4] + y.y*a.p0v1p2v3[5] + w.y;
			result[6] = x.x*a.p0v1p2v3[6] + y.x*a.p0v1p2v3[7];
			result[7] = x.y*a.p0v1p2v3[6] + y.y*a.p0v1p2v3[7];
			Shape shape;
			shape.FromShape(result);
			return shape;
		}
		Vector x, y;
		Point w;
	};

	CUDA_HOST_DEVICE static inline Matrix2x3 invert(const Matrix2x3& a) {
		auto scale = 1.0f / (a.x.x * a.y.y - a.x.y * a.y.x);
		auto x = Vector(a.y.y, -a.x.y) * scale;
		auto y = Vector(-a.y.x, a.x.x) * scale;
		auto w = Point(x * -a.w.x - y * a.w.y);
		return Matrix2x3(x, y, w);
	}

	struct GridRef {
		Matrix2x3 objectFromParent;
		int codeindex;
		unsigned nodeIndex;
	};

	struct GlyphGridCell {
		int offset;
		int count;
	};

	struct TextGridCell : public GlyphGridCell {
		int offset;
		int count;
	};

	struct shape_ptr {
		int ptr;
		int id;
	};

	struct GlyphGrid {
		size_t hres, vres;
		float lowest_x;
		float lowest_y;
		float highest_x;
		float highest_y;
		int ptr_fixup;
		int first_cell;
		bool null;
	};

	struct GlyphPtr {
		int instance;
		int ptr;
		int offset;
		int count;
		Matrix2x3 transform_to_local;
	};

	struct TextGrid {
		int hres, vres;
		float lowest_x;
		float lowest_y;
		float highest_x;
		float highest_y;
		TextGridCell *cells;
		GlyphPtr *glyph_ptrs;
		GlyphGrid *glyph_grids;
		GlyphGridCell *glyph_grid_cells;
		Shape *shapes;
		shape_ptr *shape_ptrs;
		int cell_count;
		int glyph_ptr_count;
		int glyph_grid_count;
		int glyph_grid_cell_count;
		int shape_count;
		int shape_ptr_count;
	};
}