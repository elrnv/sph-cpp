// Copyright (c) 2011, Regents of the University of Utah
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the <organization> nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "marchingTet.H"

double MarchingTet::catmullrom(const SlArray3D<double> &phi, SlVector3 &x){
	double value = 0;
	int i, j, k, l, m, n, p, q, r;
	SlVector3 z = x - lc_;
	i = (int)floor(z[0]/h_);
	j = (int)floor(z[1]/h_);
	k = (int)floor(z[2]/h_);
	SlVector3 y(i * h_, j * h_, k * h_);
	SlVector3 d = (z - y) / h_;
	double w[3][4];
	double u;
	for (l = 0; l < 3; l++) {
		u = d[l];
		w[l][0] = 0.5*(-u*u*u+2*u*u-u);
		w[l][1] = 0.5*(3*u*u*u-5*u*u+2);
		w[l][2] = 0.5*(-3*u*u*u+4*u*u+u);
		w[l][3] = 0.5*(u*u*u-u*u);
	}
	for (n = -1; n < 3; n++) {
		for (m = -1; m < 3; m++){
			for (l = -1; l < 3; l++) {
				p = i + l;
				q = j + m;
				r = k + n;
				if(p < 0) p = 0;
				if(q < 0) q = 0;
				if(r < 0) r = 0;
				if(p >= nx_) p = nx_ - 1;
				if(q >= ny_) q = ny_ - 1;
				if(r >= nz_) r = nz_ - 1;
				value += w[0][l+1] * w[1][m+1] * w[2][n+1] * phi(p, q, r);
			}
		}
	}
	return value;
}

SlVector3 MarchingTet::createVertexCR(const SlArray3D<double> &phi, SlVector3 &a, SlVector3 &b, double aPhi, double bPhi) {
	if (fabs(bPhi - aPhi) < 1e-8) return ((a + b) * 0.5);
	else {
		double mu = -aPhi / (bPhi - aPhi);
		if(mu < 1E-08) return (a + 1E-08);
		else if((1 - mu) < 1E-08) return (b + 1E-08);
		else {
			SlVector3 &left = a, &right = b, midPoint = (a + b) * 0.5;
			while(sqrMag(right - left) > 1E-10) {
				midPoint = (a + b) * 0.5;
				if((catmullrom(phi, left) * catmullrom(phi, midPoint)) < 0.0) right = midPoint;
				else if((catmullrom(phi, right) * catmullrom(phi, midPoint)) < 0.0) left = midPoint;
				else break;
			}
			return midPoint;
		}
	}
}

SlVector3 MarchingTet::createVertex(SlVector3 &a, SlVector3 &b, double aPhi, double bPhi){
	if (fabs(bPhi - aPhi) < 1e-8) return ((a + b) * 0.5);
	else {
		double mu = -aPhi / (bPhi - aPhi);
		if(mu < 1E-08) return (a + 1E-08);
		else if((1 - mu) < 1E-08) return (b + 1E-08);
		else return (a + mu * (b - a));
	}
}

MarchingTet::MarchingTet(int nx, int ny, int nz, double h, const SlVector3 &lc) {
	nx_ = nx;
	ny_ = ny;
	nz_ = nz;
	h_ = h;
	lc_ = lc;
	xFace.allocate(nx, ny - 1, nz - 1);
	yFace.allocate(nx - 1, ny, nz - 1);
	zFace.allocate(nx - 1, ny - 1, nz);
	xEdge.allocate(nx, ny - 1, nz - 1);
	yEdge.allocate(nx - 1, ny, nz - 1);
	zEdge.allocate(nx - 1, ny - 1, nz);
	for(int i = 0; i < nx; i++)
		for(int j = 0; j < ny - 1; j++)
			for(int k = 0; k < nz - 1; k++) {
        xFace(i, j, k) = xEdge(i, j, k) = -1;
			}
	for(int i = 0; i < nx - 1; i++)
		for(int j = 0; j < ny; j++)
			for(int k = 0; k < nz - 1; k++) {
        yFace(i, j, k) = yEdge(i, j, k) = -1;
			}
	for(int i = 0; i < nx - 1; i++)
		for(int j = 0; j < ny - 1; j++)
			for(int k = 0; k < nz; k++) {
        zFace(i, j, k) = zEdge(i, j, k) = -1;
			}
}

bool MarchingTet::doTet(int e1, int e2, int e3, int e4, int e5, int e6, double val0, double val1,
                        double val2, double val3, std::vector<SlTri> &triangles, std::vector<SlVector3> &meshPts) {
	int index(0);
	if (val0 > 0.0) index += 8;
	if (val1 > 0.0) index += 4;
	if (val2 > 0.0) index += 2;
	if (val3 > 0.0) index += 1;
	switch (index) {
	case 1:
		triangles.push_back(SlTri(e5,e3,e6));
		break;
	case 2:
		triangles.push_back(SlTri(e2,e4,e6));
		break;
	case 3:
		triangles.push_back(SlTri(e3,e4,e5));
		triangles.push_back(SlTri(e3,e2,e4));
		break;
	case 4:
		triangles.push_back(SlTri(e1,e5,e4));
		break;
	case 5:
		triangles.push_back(SlTri(e3,e4,e1));
		triangles.push_back(SlTri(e3,e6,e4));
		break;
	case 6:
		triangles.push_back(SlTri(e1,e6,e2));
		triangles.push_back(SlTri(e1,e5,e6));
		break;
	case 7:
		triangles.push_back(SlTri(e1,e3,e2));
		break;
	case 8:
		triangles.push_back(SlTri(e1,e2,e3));
		break;
	case 9:
		triangles.push_back(SlTri(e1,e6,e5));
		triangles.push_back(SlTri(e1,e2,e6));
		break;
	case 10:
		triangles.push_back(SlTri(e1,e6,e3));
		triangles.push_back(SlTri(e1,e4,e6));
		break;
	case 11:
		triangles.push_back(SlTri(e1,e4,e5));
		break;
	case 12:
		triangles.push_back(SlTri(e3,e4,e2));
		triangles.push_back(SlTri(e3,e5,e4));
		break;
	case 13:
		triangles.push_back(SlTri(e6,e4,e2));
		break;
	case 14:
		triangles.push_back(SlTri(e5,e6,e3));
		break;
	}
	return true;
}

void MarchingTet::buildTriangleMesh(const SlArray3D<double> &phi, std::vector<SlTri> &triangles, std::vector<SlVector3> &meshPts){
	triangles.clear();
	meshPts.clear();
	//Traverse x faces to look for intersections.
	int i, j, k;
	SlVector3 a;
	for(i = 0, a[0] = lc_[0]; i < nx_; i++, a[0] = a[0] + h_) {
		for(j = 0, a[1] = lc_[1]; j < ny_ - 1; j++, a[1] = a[1] + h_) {
			for(k = 0, a[2] = lc_[2]; k < nz_ - 1; k++, a[2] = a[2] + h_) {
				double bVal = phi(i, j, k + 1), cVal = phi(i, j + 1, k);
				if(bVal * cVal < 0.0) {
					SlVector3 b = a, c = a;
					c[1] += h_;
					b[2] += h_;
					meshPts.push_back(createVertexCR(phi, b, c, bVal, cVal));
					//meshPts.push_back(createVertex(b, c, bVal, cVal));
					xFace(i, j, k) = (int)meshPts.size() - 1;
				}
			}
		}
	}
	
	//Traverse y faces to look for intersections.
	for(i = 0, a[0] = lc_[0]; i < nx_ - 1; i++, a[0] = a[0] + h_) {
		for(j = 0, a[1] = lc_[1]; j < ny_; j++, a[1] = a[1] + h_) {
			for(k = 0, a[2] = lc_[2] ; k < nz_ - 1; k++, a[2] = a[2] + h_) {
				double bVal = phi(i, j, k + 1), cVal = phi(i + 1, j, k);
				if (bVal * cVal < 0.0) {
					SlVector3 b = a, c = a;
					c[0] += h_;
					b[2] += h_;
					meshPts.push_back(createVertexCR(phi, b, c, bVal, cVal));
					//meshPts.push_back(createVertex(b, c, bVal, cVal));
					yFace(i, j, k) = (int)meshPts.size() - 1;
				}
			}
		}
	}
	
	//Traverse z faces to look for intersections.
	for(i = 0, a[0] = lc_[0]; i < nx_ - 1; i++, a[0] = a[0] + h_) {
		for(j = 0, a[1] = lc_[1]; j < ny_ - 1; j++, a[1] = a[1] + h_) {
			for(k = 0, a[2] = lc_[2]; k < nz_; k++, a[2] = a[2] + h_) {
				double bVal = phi(i + 1, j, k), cVal = phi(i, j + 1, k);
				if (bVal * cVal < 0.0) {
					SlVector3 b = a, c = a;
					c[1] += h_;
					b[0] += h_;
					meshPts.push_back(createVertexCR(phi, b, c, bVal, cVal));
					//meshPts.push_back(createVertex(b, c, bVal, cVal));
					zFace(i, j, k) = (int)meshPts.size() - 1;
				}
			}
		}
	}
	
	//Traverse x Edges to look for intersections
	for(i = 0, a[0] = lc_[0]; i < nx_ - 1; i++, a[0] = a[0] + h_) {
		for(j = 0, a[1] = lc_[1]; j < ny_; j++, a[1] = a[1] + h_) {
			for(k = 0, a[2] = lc_[2]; k < nz_; k++, a[2] = a[2] + h_) {
				double bVal = phi(i, j, k), cVal = phi(i + 1, j, k);
				if (bVal * cVal < 0.0) {
					SlVector3 b = a, c = a;
					c[0] += h_;
					meshPts.push_back(createVertexCR(phi, b, c, bVal, cVal));
					//meshPts.push_back(createVertex(b, c, bVal, cVal));
					xEdge(i, j, k) = (int)meshPts.size() - 1;
				}
			}
		}
	}
	
	//Traverse y Edges to look for intersections
	for(i = 0, a[0] = lc_[0]; i < nx_; i++, a[0] = a[0] + h_) {
		for(j = 0, a[1] = lc_[1]; j < ny_ - 1; j++, a[1] = a[1] + h_) {
			for(k = 0, a[2] = lc_[2]; k < nz_; k++, a[2] = a[2] + h_) {
				double bVal = phi(i, j, k), cVal = phi(i, j + 1, k);
				if (bVal * cVal < 0.0) {
					SlVector3 b = a, c = a;
					c[1] += h_;
					meshPts.push_back(createVertexCR(phi, b, c, bVal, cVal));
					//meshPts.push_back(createVertex(b, c, bVal, cVal));
					yEdge(i, j, k) = (int)meshPts.size() - 1;
				}
			}
		}
	}
	
	//Traverse z Edges to look for intersections
	for(i = 0, a[0] = lc_[0]; i < nx_; i++, a[0] = a[0] + h_) {
		for(j = 0, a[1] = lc_[1]; j < ny_; j++, a[1] = a[1] + h_) {
			for(k = 0, a[2] = lc_[2]; k < nz_ - 1; k++, a[2] = a[2] + h_){
				double bVal = phi(i, j, k), cVal = phi(i, j, k + 1);
				if (bVal * cVal < 0.0) {
					SlVector3 b = a, c = a;
					c[2] += h_;
					meshPts.push_back(createVertexCR(phi, b, c, bVal, cVal));
					//meshPts.push_back(createVertex(b, c, bVal, cVal));
					zEdge(i, j, k) = (int)meshPts.size() - 1;
				}
			}
		}
	}
	
	for(i = 0, a[0] = lc_[0]; i < nx_ - 2; i++, a[0] = a[0] + h_) {
		for(j = 0, a[1] = lc_[1]; j < ny_ - 2; j++, a[1] = a[1] + h_) {
			for(k = 0, a[2] = lc_[2]; k < nz_ - 2; k++, a[2] = a[2] + h_){
				int ip1 = i + 1, jp1 = j + 1, kp1 = k + 1, e04, e26, e15, e37, e01, e23, e45, e67, e02,
					e46, e13, e57, e12, e56, e24, e35, e14, e36, e16(-1);
				double  val0 = phi(i, j, k), val1 = phi(i, j, kp1), val2 = phi(i, jp1, k), val3 = phi(i, jp1, kp1),
					val4 = phi(ip1, j, k), val5 = phi(ip1, j, kp1), val6 = phi(ip1, jp1, k), val7 = phi(ip1, jp1, kp1);
				if (val1 * val6 < 0.0) {
					SlVector3 b = a, c = a;
					b[2] += h_;
					c[0] += h_;
					c[1] += h_;
					meshPts.push_back(createVertexCR(phi, b, c, val1, val6));
					//meshPts.push_back(createVertex(b, c, val1, val6));
					e16 = (int)meshPts.size() - 1;
				}
				e12 = xFace(i, j, k); e56 = xFace(ip1, j, k);
				e14 = yFace(i, j, k); e36 = yFace(i, jp1, k);
				e24 = zFace(i, j, k); e35 = zFace(i, j, kp1);
				e04 = xEdge(i, j, k); e26 = xEdge(i, jp1, k); e15 = xEdge(i, j, kp1); e37 = xEdge(i, jp1, kp1);
				e02 = yEdge(i, j, k);
				e46 = yEdge(ip1, j, k);
				e13 = yEdge(i, j, kp1);
				e57 = yEdge(ip1, j, kp1);
				e01 = zEdge(i, j, k); e23 = zEdge(i, jp1, k); e45 = zEdge(ip1, j, k); e67 = zEdge(ip1, jp1, k);
				doTet(e02, e04, e01, e24, e12, e14, val0, val2, val4, val1, triangles, meshPts);
				doTet(e26, e16, e46, e12, e24, e14, val6, val2, val1, val4, triangles, meshPts);
				doTet(e26, e36, e16, e23, e12, e13, val6, val2, val3, val1, triangles, meshPts);
				doTet(e46, e16, e56, e14, e45, e15, val6, val4, val1, val5, triangles, meshPts);
				doTet(e16, e36, e56, e13, e15, e35, val6, val1, val3, val5, triangles, meshPts);
				doTet(e36, e67, e56, e37, e35, e57, val6, val3, val7, val5, triangles, meshPts);
			}
		}
	}
	//std::cout<<"vertices: "<<meshPts.size()<<" triangles: "<<triangles.size()<<std::endl;
}


