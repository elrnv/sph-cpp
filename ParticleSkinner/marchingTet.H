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

#ifndef marchingTet_H
#define marchingTet_H

#include "slUtil.H"
#include "slArray.H"


class MarchingTet {
public:
	MarchingTet(int nx, int ny, int nz, double h, const SlVector3 &lc);
	void buildTriangleMesh(const SlArray3D<double> &phi, std::vector<SlTri> &triangles, std::vector<SlVector3> &meshPts);
	
protected:
	int nx_, ny_, nz_;
	double h_;
	SlVector3 lc_, uc_;
	SlArray3D<int> xFace, yFace, zFace, xEdge, yEdge, zEdge;

	SlVector3 createVertex(SlVector3 &a, SlVector3 &b, double aPhi, double bPhi);
	SlVector3 createVertexCR(const SlArray3D<double> &phi, SlVector3 &a, SlVector3 &b, double aPhi, double bPhi);
	double catmullrom(const SlArray3D<double> &phi, SlVector3 &x);
	bool doTet(int e1, int e2, int e3, int e4, int e5, int e6, double val0, double val1, double val2, double val3,
						 std::vector<SlTri> &triangles, std::vector<SlVector3> &meshPts) ;
};



#endif
