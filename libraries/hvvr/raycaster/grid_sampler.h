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
#include "vector_math.h"
#include "grid_buffers.h"

namespace hvvr {

#define SAMPLES_U 4
#define SAMPLES_V 4
	struct Range {
		CUDA_DEVICE Range() = default;
		CUDA_DEVICE Range(float lower, float upper) : lower(lower), upper(upper) {}
		CUDA_DEVICE static Range Make(float a, float b) { return Range(fminf(a, b), fmaxf(a, b)); }
		float lower, upper;
	};

	CUDA_DEVICE static Range Cubic(const Range& a) {
		auto cubic = [](float x) {
			return x;
		};
		return Range(cubic(a.lower), cubic(a.upper));
	}

	CUDA_DEVICE static Range Clamp(const Range& x, const Range& bounds = Range(0.0, 1.0)) {
		return Range(fminf(fmaxf(x.lower, bounds.lower), bounds.upper), fminf(fmaxf(x.upper, bounds.lower), bounds.upper));
	}

	struct LinearEqn {
		CUDA_DEVICE LinearEqn(float a, float b, float i_a, const Range& solve) : a(a), b(b), i_a(i_a), solve(solve) {}
		CUDA_DEVICE LinearEqn(float A, float B) {
			a = A;
			b = B;
			// TODO: fix the case where A == 0
			i_a = -1.0f / A;
			auto i_b = i_a * b;
			solve = Range::Make(i_b, i_b - i_a);
		}
		CUDA_DEVICE Range operator()(const Range& t) const { return Range(a * t.lower + b, a * t.upper + b); }
		CUDA_DEVICE LinearEqn Apply(float x) const {
			//effectively a shift of the t values for offsets not at 0 1, not used in shapesOnAPlane
			return LinearEqn(a, b + x, i_a, Range(i_a * x + solve.lower, i_a * x + solve.upper));
		}
		CUDA_DEVICE Range Solve() const { return solve; }
		float a, b;
		float i_a;
		Range solve;
	};

	struct TrapezoidEqns {
		CUDA_DEVICE TrapezoidEqns(const Grid::Shape& a)
			//y1 - y0, y0
			: y0(a.p0v1p2v3[3], a.p0v1p2v3[1]),
			//x1 - x0, x0
			x0(a.p0v1p2v3[2], a.p0v1p2v3[0]),
			//x2 - x3, x3
			x1(-a.p0v1p2v3[6], a.p0v1p2v3[4] + a.p0v1p2v3[6]) {}
		LinearEqn y0, x0, x1;
	};

	struct QuadraticEqn {
		CUDA_DEVICE QuadraticEqn(float a, float b, float c, float i_a, float i_b, float apex_sign, const Range& solve)
			: a(a), b(b), c(c), i_a(i_a), i_b(i_b), apex_sign(apex_sign), solve(solve) {}
		CUDA_DEVICE QuadraticEqn(float A, float B, float C) {
			//(y2 - y0) - 2 * (y1 - y0)
			a = float(A);
			//2 * (y1 - y0)
			b = float(B);
			//y0
			c = float(C);
			// TODO: fix the case where A == 0
			i_a = float(-1.0f / A);
			// TODO: write comment about the B == 0 branch.
			i_b = 0.5f * B * i_a;
			apex_sign = B == 0 || i_b < 0.0 ? -1.0f : 1.0f;
			auto i_c = i_b * i_b + i_a * C;
			solve = Range::Make(i_c, i_c - i_a);
		}
		CUDA_DEVICE float operator()(float t) const { return (a * t + b) * t + c; }
		CUDA_DEVICE Range operator()(const Range& t) const { return Range(operator()(t.lower), operator()(t.upper)); }
		CUDA_DEVICE float SolveApex() const { return i_b; }
		CUDA_DEVICE float GetApexSign() const { return apex_sign; }
		CUDA_DEVICE Range SolveDelta() const {
			return Range(sqrtf(fmaxf(0.0f, solve.lower)), sqrtf(fmaxf(0.0f, solve.upper)));
		}
		CUDA_DEVICE QuadraticEqn Apply(float x) const {
			return QuadraticEqn(a, b, c + x, i_a, i_b, apex_sign, Range(i_a * x + solve.lower, i_a * x + solve.upper));
		}
		float a, b, c;
		float i_a, i_b;
		float apex_sign;
		Range solve;
	};

	struct CurveEqns {
		CUDA_DEVICE CurveEqns(const Grid::Shape& a)
			//x0 - x2, x2
			: x0(-a.p0v1p2v3[6], a.p0v1p2v3[0] + a.p0v1p2v3[6]),
			//y0 - y2, y2
			y0(-a.p0v1p2v3[7], a.p0v1p2v3[1] + a.p0v1p2v3[7]),
			//(x2 - x0) -  2 * (x1 - x0), 2 * (x1 - x0),  x0
			x1(a.p0v1p2v3[6] - 2 * a.p0v1p2v3[2], 2 * a.p0v1p2v3[2], a.p0v1p2v3[0]),
			//(y2 - y0) - 2 * (y1 - y0), 2 * (y1 - y0), y0
			y1(a.p0v1p2v3[7] - 2 * a.p0v1p2v3[3], 2 * a.p0v1p2v3[3], a.p0v1p2v3[1]) {}
		LinearEqn x0, y0;
		QuadraticEqn x1, y1;
	};

	
	CUDA_DEVICE float blend(float a, float b, int mask)
	{
		return mask > 0 ? a : b;
	}

	CUDA_DEVICE float IntegrateTrap(const TrapezoidEqns& eqns) {
		//construct t bottom and t top
		auto y0 = eqns.y0;
		//clamp t bottom and t top to 0 1
		auto ty0 = Clamp(y0.Solve());
		//construct y of t bottom and y of t top and round with cubic
		auto yy0 = Cubic(y0(ty0));

		float area;
		{
			//construct t left and t right for the left side of the trap
			auto x0 = eqns.x0;
			//clamp t left and t right to t bottom and t top
			//effectively clamping t left to the square or to the point x0, and symteric ops for t right
			//tx0 lower will be the left most intersection of the left side of the trap
			//upper the right most intersection of the left side of the trap
			auto tx0 = Clamp(x0.Solve(), ty0);
			//y of t left and y of t right
			auto yx0 = Cubic(y0(tx0));
			//prepare to find midpoint of y lower y upper by adding together will multiple by 0.5 later
			auto yx0_mid = yx0.lower + yx0.upper;
			//x of t lower x of t upper
			auto xx0 = Cubic(Clamp(x0(tx0)));
			//the multiply subtract is the half the height of the left trap in y
			//so this is the rectangle of that distance times x of t lower
			//it's a positive area
			area = xx0.lower * (0.5f * yx0_mid - yy0.lower);
			//the negative multiple add is also half the height of the left trap in y
			//so this is the rectangle of that distance times x of t upper - the above area
			//flipping the top triangle of the area over it the full area of the triangle left of the right most intersection
			//it's a positive area
			area = xx0.upper * (-0.5f* yx0_mid + yy0.upper) + area;
		}
		{
			//construct t right and t left for the right side of the trap
			//the equation points left not right this time so lower it the right most point
			auto x1 = eqns.x1.Apply(0);
			//clamp t right and t left to t bottom and t top
			auto tx1 = Clamp(x1.Solve(), ty0);
			//y of t right and y of t left
			auto yx1 = Cubic(y0(tx1));
			//prepare to find midpoint of y lower y upper
			auto yx1_mid = yx1.lower + yx1.upper;
			//x of t lower x of t upper
			auto xx1 = Cubic(Clamp(x1(tx1)));
			//area of half height of trap mult right most x intersection or x3
			//it's a negative area
			area = -xx1.lower*(0.5f * yx1_mid - yy0.lower) + area;
			//area of half height of trap mult left most x intersection or x2
			//it's a negative area
			//adding the parts together yields the a negative area of the entire trapezoid
			area = -xx1.upper*(-0.5f* yx1_mid + yy0.upper) + area;
		}
		return area;
	}

	CUDA_DEVICE float IntegrateCurve(const CurveEqns& eqns) {
		//shifted t top t bottom for the bottom line
		auto y0 = eqns.y0;
		//clamp t bottom and t top to 0 1
		auto ty0 = Clamp(y0.Solve());
		//y of t bottom and y of t top
		auto yy0 = Cubic(y0(ty0));

		float area;
		{
			//construct shifted t left and t right for bottom segment
			auto x0 = eqns.x0;
			//clamp t left and t right to t bottom and t top
			//effectively clamping t left to the square or to the point x0, and symteric ops for t right
			auto tx0 = Clamp(x0.Solve(), ty0);
			//y of t left and y of t right
			auto yx0 = Cubic(y0(tx0));
			//y lower t + y of upper t
			auto yx0_mid = yx0.lower + yx0.upper;
			//x of t min x of t max
			auto xx0 = Cubic(Clamp(x0(tx0)));
			//the msub is the distance to the middle of the clipped trap in y
			//so this is the rectangle of that distance times x tmin
			area = xx0.lower * (0.5f * yx0_mid - yy0.lower);
			//the nmadd is also the distance to the middle of the clipped trap in y negated
			//so this is the rectangle of that distance times x tmax - the above area, flipping the triangle top over makes is correct
			//it's a positive area if the curve is specified properly if it's inverted then this is a negative area
			area = xx0.upper * (-0.5f * yx0_mid + yy0.upper) + area;
		}
		{
			//shift the curve
			auto y1 = eqns.y1;
			//solve for y component of apex
			auto ty1_apex = y1.SolveApex();
			//solve for the delta from apex to the two y hit points on top and bottom of box
			auto ty1_delta = y1.SolveDelta();
			//get the sign of the apex
			auto ty1_apex_sign = y1.GetApexSign();
			//if the apex is negative get tapex - tdelta.lower negated otherwise get tapex - tdelta.upper
			//also if the apex is negative get tdelta.upper negated otherwise get tapex - tdelta.lower
			//the gives us t of the hit points on the top and bottom of box
			//this pushed the hit points forward for apex behind the x0 point and backward for point in front
			auto ty1 = Clamp(Range(ty1_apex - (ty1_apex_sign * blend(ty1_delta.upper, ty1_delta.lower, ty1_apex_sign)),
				ty1_apex - (ty1_apex_sign * blend(ty1_delta.lower, ty1_delta.upper, ty1_apex_sign))));


			auto x1 = eqns.x1;
			auto tx1_apex = x1.SolveApex();
			auto tx1_delta = x1.SolveDelta();
			//because we don't have local y apex but do have local x apex we need all the points here
			// -b/2a + blah - 1 , -b/2a - blah - 1 ...  "" +  blah - 0 "" -  blah - 0
			auto tx1_lower = Clamp(Range(tx1_apex - tx1_delta.upper, tx1_apex - tx1_delta.lower), ty1);
			auto tx1_upper = Clamp(Range(tx1_apex + tx1_delta.lower, tx1_apex + tx1_delta.upper), ty1);

			//solve for areas these are only positive if the curve is specified inverted
			auto yx1_lower = Cubic(y1(tx1_lower));
			auto yx1_lower_mid = yx1_lower.lower + yx1_lower.upper;
			auto xx1_lower = Cubic(Clamp(x1(tx1_lower)));
			area = xx1_lower.lower * (0.5f * yx1_lower_mid - yy0.upper) + area;
			area = xx1_lower.upper * (-0.5f * yx1_lower_mid + yx1_lower.upper) + area;

			auto yx1_upper = Cubic(y1(tx1_upper));
			auto yx1_upper_mid = yx1_upper.lower + yx1_upper.upper;
			auto xx1_upper = Cubic(Clamp(x1(tx1_upper)));
			area = xx1_upper.lower * (0.5f * yx1_upper_mid - yx1_lower.upper) + area;
			area = xx1_upper.upper * (-0.5f * yx1_upper_mid + yy0.lower) + area;
		}
		return area;
	}

	CUDA_DEVICE float IntegrateShape(Grid::Shape shape, bool thisIsATrapezoid)
	{
		if (thisIsATrapezoid)
		{
			return IntegrateTrap(TrapezoidEqns(shape));
		}
		else
		{
			return IntegrateCurve(CurveEqns(shape));
		}
	}

	CUDA_DEVICE matrix3x3 rqd(matrix3x3 A) {
		//c*a21 - y*a22 == 0
		//c^2*a21^2 = y^2*a22^2
		//c^2 + y^2 = 1
		//c^2*a21^2 = (1 - c^2)*a22^2
		//c^2*a21^2 = a22^2 - a22^2*c^2
		//c^2*a21^2 + a22^2*c^2 = a22^2
		//c^2*(a21^2 + a22^2) = (a22^2)
		//c^2 = (a21^2 + a22^2)/(a21^2 + a22^2)
		//c = +/- a22/sqrt((a21^2 + a22^2)
		//y = +/- sqrt(1 - c^2)
		float a11 = A.m0[0];
		float a12 = A.m0[1];
		float a21 = A.m1[0];
		float a22 = A.m1[1];
		float c = a22 / sqrtf(a21*a21 + a22*a22);
		float y = sqrtf(1 - c*c);
		matrix3x3 Q = matrix3x3::identity();
		Q.m0[0] = c;
		Q.m0[1] = -y;
		Q.m1[0] = y;
		Q.m1[1] = c;
		//transpose is inverse of the rotation
		//QR = A
		//R = Q^t A
		return invert(transpose(Q) * A);
	}

	CUDA_DEVICE matrix3x3 ComputeViewFromSample(vector2 x, vector2 y) {
		if (fabs(y.x) > fabs(x.x))
		{
			vector2 temp = x;
			x = y;
			y = temp;
		}
		float xy = x.y;
		float yy = y.y;
		// TODO: still buggy for smaller distortions
		if (fabs(xy) * 30 < fabs(x.x)) xy = 0; // Add dead-zone for highly anisotripic patterns.
													   // TODO: guarantee that out.y.y > 0
		auto scale = 1 / sqrtf(xy * xy + yy * yy);
		xy *= scale;
		yy *= scale;
		matrix3x3 rMat = matrix3x3::identity();
		rMat.m0[0] = x.x * yy - y.x *x.y;
		vector2 rVec = x * xy + y * yy;
		rMat.m1[0] = rVec.x;
		rMat.m1[1] = rVec.y;
		return invert(rMat);
	}

	CUDA_DEVICE matrix3x3 RQDecomp(matrix3x3 &A) {
		//-y*a11 + c*a21 == 0
		//c^2*a21^2 = y^2*a11^2
		//c^2 + y^2 = 1
		//c^2*a21^2 = (1 - c^2)*a11^2
		//c^2*a21^2 = a11^2 - a11^2*c^2
		//c^2*a21^2 + a11^2*c^2 = a11^2
		//c^2*(a21^2 + a11^2) = (a11^2)
		//c^2 = (a21^2 + a11^2)/(a21^2 + a11^2)
		//c = a11/sqrt(a21^2 + a11^2)
		//y = 1 - a11/sqrt(a21^2 + a11^2)
		//y = ((a21^2 + a11^2) - a11^2)/sqrt(a21^2 + a11^2)
		//y = a21^2/sqrt(a21^2 + a11^2)
		float a11 = A.m0[0];
		float a12 = A.m0[1];
		float a21 = A.m1[0];
		float a22 = A.m1[1];
		float c = a11 / sqrtf(a21*a21 + a11*a11);
		float y = a21 / sqrtf(a21*a21 + a11*a11);
		matrix3x3 Q = matrix3x3::identity();
		Q.m0[0] = c;
		Q.m0[1] = y;
		Q.m1[0] = -y;
		Q.m1[1] = c;
		//transpose is inverse of the rotation
		//QR = A
		//R = Q^t A
		return invert(A*Q);
	}

	CUDA_DEVICE float sample_grid_sop(TextGridGPU textGrid, vector2 uv, vector2 dUVdX, vector2 dUVdY) {
		float coverage = 0.0f;
		//uv.y = 1 - uv.y;
		dUVdX.y *= 1.2;
		dUVdY.y *= 1.2;
		matrix3x3 R = RQDecomp(matrix3x3(vector3(dUVdX.x, dUVdX.y, 0), vector3(dUVdY.x, dUVdY.y, 0),vector3(0, 0, 1)));
		float toS = 1.0f / (textGrid.highest_x - textGrid.lowest_x);
		float toT = 1.0f / (textGrid.highest_y - textGrid.lowest_y);
		Grid::Matrix2x3 toUV(Grid::Vector(toS, 0), Grid::Vector(0, toT), Grid::Point(0, 0));
		Grid::Matrix2x3 toOrigin(Grid::Vector(1, 0), Grid::Vector(0, 1), Grid::Point(-uv.x, -uv.y));
		Grid::Matrix2x3 toIntegrate(Grid::Vector(R.m0.x, 0), Grid::Vector(R.m0.y, R.m1.y), Grid::Point(0.5, 0.5));
		//max delta in x and y
		float max_dx = fmaxf(fabsf(dUVdX.x), fabsf(dUVdY.x));
		float max_dy = fmaxf(fabsf(dUVdX.y), fabsf(dUVdY.y));
		//starting and ending barys
		float bary_min_x = fminf(1.0f, fmaxf(0.0f, uv.x - max_dx));
		float bary_max_x = fminf(1.0f, fmaxf(0.0f, uv.x + max_dx));
		float bary_min_y = fminf(1.0f, fmaxf(0.0f, uv.y - max_dy));
		float bary_max_y = fminf(1.0f, fmaxf(0.0f, uv.y + max_dy));
		//starting and ending layout space coords
		float layout_min_x = bary_min_x * (textGrid.highest_x - textGrid.lowest_x) + textGrid.lowest_x;
		float layout_max_x = bary_max_x * (textGrid.highest_x - textGrid.lowest_x) + textGrid.lowest_x;
		float layout_min_y = bary_min_y * (textGrid.highest_y - textGrid.lowest_y) + textGrid.lowest_y;
		float layout_max_y = bary_max_y * (textGrid.highest_y - textGrid.lowest_y) + textGrid.lowest_y;
		//starting ending grid idxs
		int grid_min_x = bary_min_x * textGrid.hres;
		int grid_max_x = fminf(bary_max_x * textGrid.hres + 1, textGrid.hres);
		int grid_min_y = bary_min_y * textGrid.vres;
		int grid_max_y = fminf(bary_max_y * textGrid.vres + 1, textGrid.vres);

		int glyphsSeen[8] = { -1,-1,-1,-1, -1,-1,-1,-1 };
		int glyphsCounted = 0;
		for(int x = grid_min_x; x < grid_max_x; x++)
		{
			for(int y = grid_min_y; y < grid_max_y; y++)
			{

				int gg_index = y * textGrid.hres + x;
				int offset_to_glyph_grids = textGrid.cells[gg_index].offset;
				int glyph_grid_count = textGrid.cells[gg_index].count;

				for (int g = 0; g < glyph_grid_count; g++)
				{
					Grid::GlyphPtr glyph_ptr = textGrid.glyph_ptrs[offset_to_glyph_grids + g];
					int c;
					for (c = 0; c < glyphsCounted; c++) {
						if (glyphsSeen[c] == glyph_ptr.instance)
							break;
					}
					if (glyphsSeen[c] == glyph_ptr.instance)
						continue;
					glyphsSeen[glyphsCounted++] = glyph_ptr.instance;
					Grid::GlyphGrid glyph_grid = textGrid.glyph_grids[glyph_ptr.ptr];
					Grid::Matrix2x3 transform_to_local = glyph_ptr.transform_to_local;
					//get the minimum point in the local space
					Grid::Point local_space_pnt = transform_to_local*Grid::Point(layout_min_x, layout_min_y);
					float localx = (local_space_pnt.x - glyph_grid.lowest_x) / (glyph_grid.highest_x - glyph_grid.lowest_x);
					float localy = (local_space_pnt.y - glyph_grid.lowest_y) / (glyph_grid.highest_y - glyph_grid.lowest_y);

					Grid::Point local_space_pnt_max = transform_to_local*Grid::Point(layout_max_x, layout_max_y);
					float localmaxx = (local_space_pnt_max.x - glyph_grid.lowest_x) / (glyph_grid.highest_x - glyph_grid.lowest_x);
					float localmaxy = (local_space_pnt_max.y - glyph_grid.lowest_y) / (glyph_grid.highest_y - glyph_grid.lowest_y);

					int loc_gridx = fminf((float)glyph_grid.hres - 1, fmaxf(0, (int)(localx*(float)glyph_grid.hres)));
					int loc_gridy = fminf((float)glyph_grid.vres - 1, fmaxf(0, (int)(localy*(float)glyph_grid.vres)));

					int loc_gridmaxx = fminf((float)glyph_grid.hres, fmaxf(0, (int)(localmaxx*(float)glyph_grid.hres + 1)));
					int loc_gridmaxy = fminf((float)glyph_grid.vres, fmaxf(0, (int)(localmaxy*(float)glyph_grid.vres + 1)));
					Grid::Matrix2x3 toLayout(Grid::Vector(1, 0), Grid::Vector(0, 1), Grid::Point(-transform_to_local.w.x - textGrid.lowest_x, -transform_to_local.w.y - textGrid.lowest_y));
					long long hash[8] = { 0,0,0,0,0,0,0,0 };

					for (int xx = loc_gridx; xx < loc_gridmaxx; xx++)
					{
						for (int yy = loc_gridy; yy < loc_gridmaxy; yy++)
						{
							Grid::GlyphGridCell gg_cell = textGrid.glyph_grid_cells[glyph_grid.first_cell + yy * glyph_grid.hres + xx];
							for (int c = 0; c < gg_cell.count; c++)
							{
								int shape_ptr = textGrid.shape_ptrs[gg_cell.offset + c + glyph_grid.ptr_fixup].ptr;
								long long shape_id = textGrid.shape_ptrs[gg_cell.offset + c + glyph_grid.ptr_fixup].id;
								long long idmod = shape_id % 64ll;
								int field = shape_id / 64ll;
								if (hash[field] & (1ll << idmod))
									continue;
								hash[field] |= (1ll << idmod);
								Grid::Shape shape = textGrid.shapes[shape_ptr];
								shape = toLayout*shape;
								shape = toUV*shape;
								shape = toOrigin*shape;
								shape = toIntegrate*shape;
								if(!shape.OutBounds())
									coverage += IntegrateShape(shape, shape.IsTrapazoid());
							}
						}
					}
				}
			}
		}
		return 1.0f - fminf(1,fabs(coverage));
	}
}