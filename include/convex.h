//	-------------------------------------------------------------------------
//	Smallest circumscribing k-gon. https://github.com/atyuwen/smallest-k-gon.
//	Contact : atyuwen@gmail.com
//	Author : Yuwen Wu (https://atyuwen.github.io/)
//	License : Distributed under the MIT License.
//	-------------------------------------------------------------------------

#pragma once

#include <cmath>
#include <vector>
#include <algorithm>

namespace ayw
{
	constexpr float BIG_FLOAT = 1e12f;
	constexpr float EPS_FLOAT = 1e-6f;

	struct float2
	{
		float x, y;
		float2(float _x = 0, float _y = 0) : x(_x), y(_y){}
		float length_sqr() const {return x * x + y * y;}
		float length() const {return std::sqrt(length_sqr());}
		float2& operator += (const float2 &rhs) {x += rhs.x; y += rhs.y; return *this;}
		float2& operator -= (const float2 &rhs) {x -= rhs.x; y -= rhs.y; return *this;}
		float2& operator *= (float rhs) {x *= rhs, y *= rhs; return *this;}
		friend float2 operator - (const float2 &op) {return float2(-op.x, -op.y);}
		friend float2 operator + (const float2 &lhs, const float2 &rhs) {float2 ret = lhs; return ret += rhs;}
		friend float2 operator - (const float2 &lhs, const float2 &rhs) {float2 ret = lhs; return ret -= rhs;}
		friend float2 operator * (const float2 &lhs, float rhs) {float2 ret = lhs; return ret *= rhs;}
		friend float2 operator * (float lhs, const float2 &rhs) {return rhs * lhs;}
		friend float dot(const float2 &lhs, const float2 &rhs){return lhs.x * rhs.x + lhs.y * rhs.y;}
		friend float cross(const float2 &lhs, const float2 &rhs){return lhs.x * rhs.y - lhs.y * rhs.x;}
		friend bool operator < (const float2 &lhs, const float2 &rhs){return lhs.y < rhs.y || (lhs.y == rhs.y && lhs.x < rhs.x);}
	};

	// calculate the intersection of two non-collinear segments [a, b] and [c, d]
	// return parameters(s, t) so that a + ab * s = c + cd * t;
	static float2 intersection(const float2 &a, const float2 &b, const float2 &c, const float2 &d)
	{
		float2 ab = b - a, cd = d - c, ac = c - a;
		float num = cross(ab, ac);
		float denom = cross(cd, ab);
		if (denom == 0) return float2(-BIG_FLOAT, -BIG_FLOAT);
		float t = num / denom;
		float s = 0;
		if (std::abs(ab.x) > std::abs(ab.y)) s = (ac.x + cd.x * t) / ab.x;
		else s = (ac.y + cd.y * t) / ab.y;
		return float2(s, t);
	}

	struct convex
	{
		std::vector<float2> vertices;

		// build the convex hull via Graham scan algorithm
		template <typename float2_iterator>
		bool build(float2_iterator begin, float2_iterator end)
		{
			struct polar_less
			{
				float2 ref;
				polar_less(const float2 &_ref) : ref(_ref){}

				bool operator() (const float2 &lhs, const float2 &rhs) const
				{
					float2 ld = lhs - ref, rd = rhs - ref;
					float cross_prod = cross(ld, rd);
					return cross_prod > 0 || (cross_prod == 0 && ld.length_sqr() < rd.length_sqr());
				}
			};

			if (end - begin < 3) return false;
			std::vector<float2> points(begin, end);
			std::iter_swap(points.begin(), std::min_element(points.begin(), points.end()));
			polar_less less_op(points.front());
			std::sort(points.begin() + 1, points.end(), less_op);
			vertices.clear();
			vertices.insert(vertices.end(), points.begin(), points.begin() + 3);
			for (int i = 3; i != points.size(); ++i)
			{
				float2 &cur = points[i];
				while (true)
				{
					if (vertices.size() < 2) {vertices.emplace_back(cur); break;}
					float2 prev_top = vertices.back() - *(vertices.end() - 2);
					float2 top_cur = cur - vertices.back();
					if (cross(prev_top, top_cur) > EPS_FLOAT) {vertices.emplace_back(cur); break;}
					vertices.pop_back();
				}
			}
			if (vertices.size() < 3) {vertices.clear(); return false;}
			return true;
		}

		bool contains(const float2& point) const
		{
			for (int i = 0; i != vertices.size(); ++i)
			{
				float2 v1 = vertices[i] - point;
				float2 v2 = vertices[i == vertices.size() - 1 ? 0 : i + 1] - point;
				if (cross(v1, v2) < -EPS_FLOAT) return false;
			}
			return true;
		}

		float area() const
		{
			float ret = 0;
			for (int i = 0; i != vertices.size(); ++i)
			{
				int j = (i + 1) % vertices.size();
				ret += cross(vertices[i], vertices[j]);
			}
			return ret / 2;
		}

		// clip the convex by rect [(0,0),(1,1)] via Sutherland–Hodgman algorithm
		void clip()
		{
			struct clipper
			{
				enum edge{LEFT, BOTTOM, RIGHT, TOP};

				std::vector<float2> &vertices;
				float left, bottom, right, top;

				clipper(std::vector<float2> &vertices, float left, float bottom, float right, float top)
					: vertices(vertices), left(left), bottom(bottom), right(right), top(top) {}

				bool is_inside(float2 p, edge e)
				{
					if (e == LEFT) return p.x >= left;
					if (e == BOTTOM) return p.y >= bottom;
					if (e == RIGHT) return p.x <= right;
					/* if (e == TOP) */ return p.y <= top;
				}

				float2 intersection(float2 p, float2 q, edge e)
				{
					if (e == LEFT) return float2(left, p.y + (q.y - p.y) * (p.x - left) / (p.x - q.x));
					if (e == BOTTOM) return float2(p.x + (q.x - p.x) * (p.y - bottom) / (p.y - q.y), bottom);
					if (e == RIGHT) return float2(right, p.y + (q.y - p.y) * (p.x - right) / (p.x - q.x));
					/* if (e == TOP) */ return float2(p.x + (q.x - p.x) * (p.y - top) / (p.y - q.y), top);
				}

				void run()
				{
					std::vector<float2> cliped_vertices;
					for (int i = 0; i != 4; ++i)
					{
						edge e = static_cast<edge>(i);
						cliped_vertices.clear();
						for (int j = 0; j != vertices.size(); ++j)
						{
							float2 &p = vertices[j == 0 ? vertices.size() - 1 : j - 1];
							float2 &q = vertices[j];
							bool p_inside = is_inside(p, e);
							bool q_inside = is_inside(q, e);
							if (p_inside)
							{
								if (q_inside) cliped_vertices.emplace_back(q);
								else cliped_vertices.emplace_back(intersection(p, q, e));
							}
							else if (q_inside)
							{
								cliped_vertices.emplace_back(intersection(p, q, e));
								cliped_vertices.emplace_back(q);
							}
						}
						using std::swap;
						swap(vertices, cliped_vertices);
					}
				}
			};

			if (vertices.size() < 3) return;
			clipper clipper_01(vertices, 0, 0, 1, 1);
			clipper_01.run();
			build(vertices.begin(), vertices.end());
		}

		// reduce the number of vertices with minimal area addition
		// [Alok Aggarwal et al.] "Minimum area circumscribing Polygons"
		void simplify(unsigned int num_vertices)
		{
			struct simplifier
			{
				std::vector<float2> &vertices;
				std::vector<std::vector<float> > a_cut;
				std::vector<std::vector<float> > a_complement;

				std::vector<std::vector<int> > m_balanced;
				std::vector<std::vector<float> > a_balanced;
				std::vector<std::vector<std::vector<int> > > m_flushed;
				std::vector<std::vector<std::vector<float> > > a_flushed;

				simplifier(std::vector<float2> &v, int k)
					: vertices(v)
					, a_cut(v.size(), std::vector<float>(v.size(), 0.0f))
					, a_complement(v.size(), std::vector<float>(v.size(), 0.0f))
					, m_balanced(v.size(), std::vector<int>(v.size(), -1))
					, a_balanced(v.size(), std::vector<float>(v.size(), BIG_FLOAT))
					, m_flushed(k - 3, std::vector<std::vector<int> >(v.size(), std::vector<int>(v.size(), -1)))
					, a_flushed(k - 3, std::vector<std::vector<float> >(v.size(), std::vector<float>(v.size(), BIG_FLOAT)))
				{
					initialize();
				}

				int wrap(int i)
				{
					int n = static_cast<int>(vertices.size());
					int j = i % n;
					return j < 0 ? j + n : j;
				}

				void calc_balanced(int i, int j, int l, int r)
				{
					float min_extra_area = BIG_FLOAT;
					int choose = -1;
					for (int k = l; k != wrap(r + 1) && k != j; k = wrap(k + 1))
					{
						float2 dummy_1, dummy_2;
						float extra_area = extra_area_balanced(i, j, k, dummy_1, dummy_2);
						if (extra_area < min_extra_area)
						{
							min_extra_area = extra_area;
							choose = k;
						}
					}
					m_balanced[i][j] = choose;
					a_balanced[i][j] = min_extra_area;
				}

				void calc_flushed(int h, int i, int j, int l, int r)
				{
					float min_extra_area = BIG_FLOAT;
					int choose = -1;
					for (int k = l; k != wrap(r + 1) && k != j; k = wrap(k + 1))
					{
						float extra_area = BIG_FLOAT;
						int h1 = (h + 1) / 2 - 1;
						int h2 = h - (h + 1) / 2 - 1;
						if (wrap(k - i) > h1 + 1 && wrap(j - k) > h2 + 1)
						{
							float area1 = h1 < 0 ? a_complement[i][k] : a_flushed[h1][i][k];
							float area2 = h2 < 0 ? a_complement[k][j] : a_flushed[h2][k][j];
							extra_area = area1 + area2;
						}
						if (extra_area < min_extra_area)
						{
							min_extra_area = extra_area;
							choose = k;
						}
					}
					m_flushed[h][i][j] = choose;
					a_flushed[h][i][j] = min_extra_area;
				}

				float extra_area_balanced(int i, int j, int k, /*out*/ float2 &p, /*out*/ float2 &q)
				{
					float2 a = vertices[i], b = vertices[wrap(i + 1)];
					float2 c = vertices[wrap(j + 1)], d = vertices[wrap(j)];
					float2 e = vertices[k];
					float2 ab = b - a, cd = d - c;
					float2 next = vertices[wrap(k + 1)] - e;
					float2 prev = e - vertices[wrap(k - 1)];
					float2 opt = next;
					if (k != wrap(i + 1) && std::abs(cross(ab, cd)) > EPS_FLOAT)
					{
						float t = cross(ab, 2 * e - a - c) / cross(ab, cd);
						opt = c + cd * t - e;
						if (cross(prev, opt) < 0) opt = prev;
						if (cross(next, opt) > 0) opt = next;
					}

					float2 s1 = intersection(a, b, e, e - opt);
					float2 s2 = intersection(c, d, e, e + opt);
					if (s1.x < 1 - EPS_FLOAT || s2.x < 1 - EPS_FLOAT) return BIG_FLOAT;

					p = a + ab * s1.x;
					q = c + cd * s2.x;
					float quad_area = (cross(p - b, q - b) + cross(q - b, d - b)) / 2;
					return quad_area - a_cut[i][j];
				}

				void recurse_balanced(int i, int j1, int j2)
				{
					if (wrap(j2 - j1) < 2) return;
					int m = wrap(j1 + wrap(j2 - j1) / 2);
					int l = m_balanced[i][j1] < 0 ? wrap(i + 1) : m_balanced[i][j1]; 
					int r = m_balanced[i][j2] < 0 ? wrap(m - 1) : m_balanced[i][j2];
					calc_balanced(i, m, l, r);
					recurse_balanced(i, j1, m);
					recurse_balanced(i, m, j2);
				}

				void recurse_flushed(int h, int i, int j1, int j2)
				{
					if (wrap(j2 - j1) < 2) return;
					int m = wrap(j1 + wrap(j2 - j1) / 2);
					int l = m_flushed[h][i][j1] < 0 ? wrap(i + 1) : m_flushed[h][i][j1]; 
					int r = m_flushed[h][i][j2] < 0 ? wrap(m - 1) : m_flushed[h][i][j2];
					calc_flushed(h, i, m, l, r);
					recurse_flushed(h, i, j1, m);
					recurse_flushed(h, i, m, j2);
				}

				void reconstruct_flushed(int h, int i, int j, std::vector<int> &out_edges)
				{
					if (h < 0 || wrap(j - i) <= h + 1) return;
					int h1 = (h + 1) / 2 - 1;
					int h2 = h - (h + 1) / 2 - 1;
					int k = m_flushed[h][i][j];
					reconstruct_flushed(h1, i, k, out_edges);
					out_edges.emplace_back(k);
					reconstruct_flushed(h2, k, j, out_edges);
				}

				void initialize()
				{
					// init cut area
					for (int i = 0; i != vertices.size(); ++i)
					{
						for (int j = wrap(i + 3); j != i; j = wrap(j + 1))
						{
							float2 v1 = vertices[wrap(j - 1)] - vertices[wrap(i + 1)];
							float2 v2 = vertices[j] - vertices[wrap(j - 1)];
							float area_augment = cross(v1, v2) / 2;
							a_cut[i][j] = a_cut[i][wrap(j - 1)] + area_augment;
						}
					}

					// init complement area
					for (int i = 0; i != vertices.size(); ++i)
					{
						for (int j = wrap(i + 2); j != i; j = wrap(j + 1))
						{
							float2 a = vertices[i], b = vertices[wrap(i + 1)];
							float2 c = vertices[wrap(j + 1)], d = vertices[j];
							float2 s = intersection(a, b, c, d);
							if (s.x < 1 - EPS_FLOAT || s.y < 1 - EPS_FLOAT) a_complement[i][j] = BIG_FLOAT;
							else
							{
								float2 e = a + (b - a) * s.x;
								float2 v1 = d - e, v2 = b - e;
								a_complement[i][j] = cross(v1, v2) / 2 - a_cut[i][j];
							}
						}
					}

					// compute optimal one-sided chains
					for (int i = 0; i != vertices.size(); ++i)
					{
						int j1 = wrap(i + 2), j2 = wrap(i - 2);
						calc_balanced(i, j1, wrap(i + 1), wrap(j1 - 1));
						calc_balanced(i, j2, wrap(i + 1), wrap(j2 - 1));
						recurse_balanced(i, j1, j2);
					}

					// compute optimal flush chains
					for (int h = 0; h != m_flushed.size(); ++h)
					{
						for (int i = 0; i != vertices.size(); ++i)
						{
							int j1 = wrap(i + 2), j2 = wrap(i - 2);
							calc_flushed(h, i, j1, wrap(i + 1), wrap(j1 - 1));
							calc_flushed(h, i, j2, wrap(i + 1), wrap(j2 - 1));
							recurse_flushed(h, i, j1, j2);
						}
					}
				}

				void run()
				{
					// find the optimal cut pair (s, t)
					float min_extra_area = BIG_FLOAT;
					int s = -1, t = -1;
					int h = (int)a_flushed.size() - 1;
					for (int i = 0; i != vertices.size(); ++i)
					{
						for (int j = 0; j != vertices.size(); ++j)
						{
							float extra_area = a_balanced[i][j] + a_flushed[h][j][i];
							if (extra_area < min_extra_area)
							{
								min_extra_area = extra_area;
								s = i; t = j;
							}
						}
					}

					// reconstruct the optimal k-gon
					std::vector<float2> new_vertices(2);
					extra_area_balanced(s, t, m_balanced[s][t], new_vertices[0], new_vertices[1]);

					std::vector<int> flushed_edges;
					flushed_edges.emplace_back(t);
					reconstruct_flushed(h, t, s, flushed_edges);
					flushed_edges.emplace_back(s);

					for (int i = 0; i < static_cast<int>(flushed_edges.size()) - 1; ++i)
					{
						int e1 = flushed_edges[i];
						int e2 = flushed_edges[i + 1];
						float2 a = vertices[e1], b = vertices[wrap(e1 + 1)];
						float2 c = vertices[e2], d = vertices[wrap(e2 + 1)];
						float2 p = intersection(a, b, c, d);
						new_vertices.emplace_back(a + (b - a) * p.x);
					}

					using std::swap;
					swap(vertices, new_vertices);
				}
			};

			if (num_vertices < 4 || vertices.size() < num_vertices) return;
			simplifier simpifier_inst(vertices, num_vertices);
			simpifier_inst.run();
			build(vertices.begin(), vertices.end());
		}

	};

}
