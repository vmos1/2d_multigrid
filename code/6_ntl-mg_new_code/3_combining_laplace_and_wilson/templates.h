#pragma once

struct site{ 
    int x,y;
    };

/*  ### Define 3D,2D and 1D array templates ## */ 
template<typename T>
struct Array3D {
    std::shared_ptr<T> ptr;
    int N0, N1, N2, size;

    Array3D() = default;
    Array3D(int N0_, int N1_, int N2_)
     : N0(N0_) ,N1(N1_), N2(N2_), size(N0_*N1_*N2_) {
        ptr = std::shared_ptr<T>(new T[size], [](T const * p) { delete[] p; });
        }
    Array3D(const Array3D&) = default;
    Array3D(Array3D&&) = default;
    Array3D& operator=(const Array3D&) = default;
    Array3D& operator=(Array3D&&) = default;

    T& operator()(const int& n0, const int& n1, const int& n2) 
    {
        return ptr.get()[n0 + N0 * n1 + N0 * N1 * n2];
    }
};

template<typename T>
struct Array2D {
    std::shared_ptr<T> ptr;
    int N0, N1, size;

    Array2D() = default;
    Array2D(int N0_, int N1_)
     : N0(N0_) ,N1(N1_), size(N0_*N1_) {
        ptr = std::shared_ptr<T>(new T[size], [](T const * p) { delete[] p; });
        }
    Array2D(const Array2D&) = default;
    Array2D(Array2D&&) = default;
    Array2D& operator=(const Array2D&) = default;
    Array2D& operator=(Array2D&&) = default;

    T& operator()(const int& n0, const int& n1) 
    {
        return ptr.get()[n0 + N0 * n1];
    }
};

template<typename T>
struct Array1D {
    std::shared_ptr<T> ptr;
    int N0, size;

    Array1D() = default;
    Array1D(int N0_)
     : N0(N0_) , size(N0_) {
        ptr = std::shared_ptr<T>(new T[size], [](T const * p) { delete[] p; });
        }
    Array1D(const Array1D&) = default;
    Array1D(Array1D&&) = default;
    Array1D& operator=(const Array1D&) = default;
    Array1D& operator=(Array1D&&) = default;

    T& operator()(const int& n0) 
    {
        return ptr.get()[n0];
    }
};

// Matrix of arbitrary size:
typedef Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> ColorMatrix;
// Vector of arbitrary size:
typedef Eigen::Matrix<complex<double>, Eigen::Dynamic, 1> ColorVector;

// typedef Array3D<Complex> CArr3D;
// typedef Array2D<Complex> CArr2D;
// typedef Array1D<Complex> CArr1D;

typedef Array1D<ColorVector> VArr1D;
typedef Array2D<ColorVector> VArr2D;
typedef Array2D<ColorMatrix> MArr2D;
typedef Array1D<ColorMatrix> MArr1D;
