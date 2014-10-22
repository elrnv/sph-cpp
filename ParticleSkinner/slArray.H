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

#ifndef SLARRAY_H
#define SLARRAY_H
#include <iostream>
#include "slIO.H"

template <typename T>
class SlArray3D {
	T *data;
	unsigned int _nx, _ny, _nz, nynz;
public:
	inline T &operator()(const unsigned int &i, const unsigned int &j, const unsigned int &k);
	inline T &operator()(const unsigned int &i, const unsigned int &j, const unsigned int &k) const;
	/*inline*/ const T &operator=(const T &d);
	inline const SlArray3D<T> &operator=(const SlArray3D<T> &d);
	inline SlArray3D();
	inline SlArray3D(const unsigned int &nx, const unsigned int &ny, const unsigned int &nz);
	inline ~SlArray3D();
	inline bool allocate(const unsigned int &nx, const unsigned int &ny, const unsigned int &nz);
	unsigned int nx() const {return _nx;};
	unsigned int ny() const {return _ny;};
	unsigned int nz() const {return _nz;};
};

template <typename T> std::istream &operator>>(std::istream &strm, SlArray3D<T> &v);
template <typename T> std::ostream &operator<<(std::ostream &strm, const SlArray3D<T> &v);

template <typename T>
inline bool SlArray3D<T>::allocate(const unsigned int &x, const unsigned int &y, const unsigned int &z) {
	if (data!=NULL && _nx==x && _ny==y && _nz==z) return true;
	_nx = x; 
	_ny = y;
	_nz = z;
	nynz = _ny*_nz;
	if (data) delete []data;
	data = new T[_nx*nynz];
	return true;
}

template<class T>
inline SlArray3D<T> operator*(const SlArray3D<T> &d, T x) {
	SlArray3D<T> y(d.nx(), d.ny(), d.nz());
	for (unsigned int i=0; i<d.nx(); i++)
		for (unsigned int j=0; j<d.ny(); j++)
			for (unsigned int k=0; k<d.nz(); k++)
				y(i,j,k) = d(i,j,k)*x;
	return y;
}

template <class T>
/*inline*/ const T &SlArray3D<T>::operator=(const T &x) {
	T *d = data;
	for (unsigned int i=0; i<_nx; i++) {
		for (unsigned int j=0; j<_ny; j++) {
			for (unsigned int k=0; k<_nz; k++, d++) {
				(*d) = x;
			}
		}
	}
	return x;
}


template <class T>
inline const SlArray3D<T> &SlArray3D<T>::operator=(const SlArray3D<T> &x) {
	T *d = data;
	T *y = x.data;
	for (unsigned int i=0; i<_nx; i++) {
		for (unsigned int j=0; j<_ny; j++) {
			for (unsigned int k=0; k<_nz; k++, d++, y++) {
				(*d) = (*y);
			}
		}
	}
	return x;
}

template <class T>
inline T &SlArray3D<T>::operator()(const unsigned int &i, const unsigned int &j, const unsigned int &k) {
	return data[i*nynz+j*_nz+k];
}

template <class T>
inline T &SlArray3D<T>::operator()(const unsigned int &i, const unsigned int &j, const unsigned int &k) const {
	return data[i*nynz+j*_nz+k];
}

template <class T>
inline SlArray3D<T>::SlArray3D() {
	data = NULL;
}

template <class T>
inline SlArray3D<T>::SlArray3D(const unsigned int &x, const unsigned int &y, const unsigned int &z) {
	data = NULL;
	allocate(x,y,z);
}

template <class T>
inline SlArray3D<T>::~SlArray3D() {
	if (data) delete [] data;
}

template <class T>
std::istream &operator>>(std::istream &strm, SlArray3D<T> &v) {
	unsigned int nx,ny,nz,i,j,k;
	std::ios::fmtflags orgFlags = strm.setf(std::ios::skipws);

	eatStr("[", strm);
	strm >> nx;
	eatStr(",", strm);
	strm >> ny;
	eatStr(",", strm);
	strm >> nz;
	eatStr("]", strm);
	
	eatStr("[", strm);
	v.allocate(nx, ny, nz);
	for (i=0; i<nx; i++) {
		for (j=0; j<ny; j++) {
			for (k=0; k<nz; k++) {
				strm >> v(i,j,k);
			}
		}
	}
	eatStr("]", strm);
	strm.flags(orgFlags);
	return strm;
}

template <typename T>
std::ostream &operator<<(std::ostream &strm, const SlArray3D<T> &v) {
	strm << "[";
	strm << v.nx() << "," << v.ny() << "," <<v.nz() << "]";
	unsigned int i, j, k;
	strm<<"["<<std::endl;
	for (i=0; i<v.nx(); i++) {
		for (j=0; j<v.ny(); j++) {
			for (k=0; k<v.nz(); k++) {
				strm << v(i,j,k) << " ";
			}
		}
	}
	strm<<"]\n";
	return strm;
}
#endif

