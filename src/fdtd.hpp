#ifndef FDTD_HPP
#define FDTD_HPP

#include <cstdint>
#include <arrayfire.h>
#include <vector>

#include "Shape2D.hpp"
#include "material.hpp"

// ~Speed of light
#define SOL 300000000.0

class FDTD
{

public:
    FDTD();

    FDTD( uint32_t dim_x, uint32_t dim_y, uint32_t dim_z, float minWavelength );

    ~FDTD();

    void timeStep( void );

    void setupSourceNode( uint32_t x, uint32_t y, uint32_t z );
    void setupSourceNode( void );

    uint32_t getXDim( void );
    uint32_t getYDim( void );
    uint32_t getZDim( void );

    af::array getYZPlane( uint32_t x );
    af::array getYZPlane( void );
    af::array getAvgPlane( void );

    double getCurrentTime( void );

    void add2dAntenna( Shape2D outline, float thickness, Vec3 centerPoint, Material mat );
    void setMinimumWavelength( float minWavelength );
    void setMaxFreq( float maxFreq );
    void setTimeAvg( void );

private:

    void updateEField( void );
    void updateHField( void );
    void enforceABC( void );
    void enforcePEC( void );
    void applySourceNode( void );
    void applyArrayNodes( void );

private:

    uint64_t mCurTimeStep;

    uint32_t mDepthX;
    uint32_t mDepthY;
    uint32_t mDepthZ; 

    uint32_t mSourceLocX;
    uint32_t mSourceLocY;
    uint32_t mSourceLocZ;
    float mSourceAmplitude;

    float mABCcoeff;
    float mFSI;
    float mCourantLimit;
    // mDeltaS = delta_x,y,z as calculated by courant condition taken from minimum wavelength (min lambda / 20)
    // delta_t = delta_x / (c*sqrt(3))
    double mDeltaS;
    double mDeltaT;
    // min wavelength
    float mMinWavelength;

    bool mTimeAvg;
    af::array mAvgPlane;

    // E field definitions
    af::array e_x;
    af::array e_y;
    af::array e_z;
    // H field definitions
    af::array h_x;
    af::array h_y;
    af::array h_z;
    // E field coefficients
    af::array coeff_e_xe;
    af::array coeff_e_xh;
    af::array coeff_e_ye;
    af::array coeff_e_yh;
    af::array coeff_e_ze;
    af::array coeff_e_zh;
    // H field coefficients
    af::array coeff_h_xe;
    af::array coeff_h_xh;
    af::array coeff_h_ye;
    af::array coeff_h_yh;
    af::array coeff_h_ze;
    af::array coeff_h_zh;
    // ABC Planes
    af::array abc_xp0_ty;
    af::array abc_xp0_tz;
    af::array abc_xp1_ty;
    af::array abc_xp1_tz;

    af::array abc_yp0_tx;
    af::array abc_yp0_tz;
    af::array abc_yp1_tx;
    af::array abc_yp1_tz;

    af::array abc_zp0_tx;
    af::array abc_zp0_ty;
    af::array abc_zp1_tx;
    af::array abc_zp1_ty;

};


#endif
