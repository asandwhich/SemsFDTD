#include "fdtd.hpp"


// courant number = S_c = c*dt/dx
// time step / spacial discretization
// c*dt <= dx (propagation between steps)

// in three dimensions courant limit is 1 / sqrt( 3 )

// when number is anything besides unity, the
// fields become dispersive.
// Field arrangement:
//      ez------hx
//     /|      / |
//    hy+-----+  |
//    | |_____|__ey
//    |/      | /
//    ex______hz


// stability condition:
// delta_t = n / ( 2 * c * sqrt( 1/delta_x^2 + 1/delta_y^2 + 1/delta_z^2 ) ) where n = dim_num
// cell size delta_x = lambda_min / 20
// f = 10 GHz, lambda = 3 cm, delta_x = 1.5 mm
// 50 cm / delta_x = 333.33
//
// C_Hxh = ( 1 - sigma_n*delta_t/ 2*mu ) / ( 1 + sigma_m*delta_t/ 2*mu ) = C_Hyh = C_Hzh
// C_Hxe = 1/(1 + sigma_m*delta_t / 2*mu ) * delta_t/(mu*delta) = C_Hye = C_Hze
// Same as above, but swap E for H, epsilon for mu, and sigma for sigma_m
//
// 4 GB limit around 300 x 300 x300
// want 6-8 GHz, lambda = 3.75 cm
// delta_x = 1.8 mm
// 1.5 m / delta_x = 800
// 0.3 m / delta_x = 60
// 0.86 m / delta_x = ~475
//
// delta_x = delta_y = delta_z
// then delta_t = delta_x / (c*sqrt(3)) = 0.0018 m / sqrt(3) * 3e8 m/s = 3.4641016*10^-12 = ~3.45 ps
// 1 / delta_t = ~288.7 GHz

FDTD::FDTD()
    : FDTD( 100, 100, 100, 0.0375 )
{
}


FDTD::FDTD( uint32_t dim_x, uint32_t dim_y, uint32_t dim_z, float minWavelength )
    : mDepthX( dim_x )
    , mDepthY( dim_y )
    , mDepthZ( dim_z )
    , mCurTimeStep( 0 )
    , mMinWavelength( minWavelength )
    , mTimeAvg( false )
    , e_x( dim_x - 1, dim_y, dim_z )
    , e_y( dim_x, dim_y - 1, dim_z )
    , e_z( dim_x, dim_y, dim_z - 1 )
    , h_x( dim_x, dim_y - 1, dim_z - 1 )
    , h_y( dim_x - 1, dim_y, dim_z - 1 )
    , h_z( dim_x - 1, dim_y - 1, dim_z )
    , coeff_e_xe( dim_x - 1, dim_y, dim_z )
    , coeff_e_xh( dim_x - 1, dim_y, dim_z )
    , coeff_e_ye( dim_x, dim_y - 1, dim_z )
    , coeff_e_yh( dim_x, dim_y - 1, dim_z )
    , coeff_e_ze( dim_x, dim_y, dim_z - 1 )
    , coeff_e_zh( dim_x, dim_y, dim_z - 1 )
    , coeff_h_xe( dim_x, dim_y - 1, dim_z - 1 )
    , coeff_h_xh( dim_x, dim_y - 1, dim_z - 1 )
    , coeff_h_ye( dim_x - 1, dim_y, dim_z - 1 )
    , coeff_h_yh( dim_x - 1, dim_y, dim_z - 1 )
    , coeff_h_ze( dim_x - 1, dim_y - 1, dim_z )
    , coeff_h_zh( dim_x - 1, dim_y - 1, dim_z )
    , abc_xp0_ty( 1, dim_y - 1, dim_z )
    , abc_xp0_tz( 1, dim_y, dim_z - 1 )
    , abc_xp1_ty( 1, dim_y - 1, dim_z )
    , abc_xp1_tz( 1, dim_y, dim_z - 1 )
    , abc_yp0_tx( dim_x - 1, 1, dim_z )
    , abc_yp0_tz( dim_x, 1, dim_z - 1 )
    , abc_yp1_tx( dim_x - 1, 1, dim_z )
    , abc_yp1_tz( dim_x, 1, dim_z - 1 )
    , abc_zp0_tx( dim_x - 1, dim_y, 1 )
    , abc_zp0_ty( dim_x, dim_y - 1, 1 )
    , abc_zp1_tx( dim_x - 1, dim_y, 1 )
    , abc_zp1_ty( dim_x, dim_y - 1, 1 )
    , mAvgPlane( 1, dim_y, dim_z ) // TODO: make flexible, current assumes e_x
{
    // Courant limit in 3d space is 1/sqrt(3)
    mCourantLimit = 1.0f / sqrt( 3.0f );

    // Initialize ABC params
    mABCcoeff = ( mCourantLimit - 1.0f ) / ( mCourantLimit + 1.0f );
    mFSI = 377.0f;

    mDeltaS = (double) mMinWavelength / 20;
    mDeltaT = ( mDeltaS * (double)mCourantLimit ) / SOL;

    // Initialize e coefficients to free space for now
    coeff_e_xe = 1.0f;
    coeff_e_xh = mCourantLimit * mFSI;
    coeff_e_ye = 1.0f;
    coeff_e_yh = mCourantLimit * mFSI;
    coeff_e_ze = 1.0f;
    coeff_e_zh = mCourantLimit * mFSI;
    // Initialize h coefficients to free space for now as well
    coeff_h_xe = mCourantLimit / mFSI;
    coeff_h_xh = 1.0f;
    coeff_h_ye = mCourantLimit / mFSI;
    coeff_h_yh = 1.0f;
    coeff_h_ze = mCourantLimit / mFSI;
    coeff_h_zh = 1.0f;

    setupSourceNode();
}

FDTD::~FDTD()
{
}

void FDTD::setMinimumWavelength( float minWavelength )
{
    mMinWavelength = minWavelength;
}

void FDTD::setMaxFreq( float maxFreq )
{
    mMinWavelength = SOL / maxFreq;
}

void FDTD::setupSourceNode( uint32_t x, uint32_t y, uint32_t z )
{
    mSourceLocX = x;
    mSourceLocY = y;
    mSourceLocZ = z;

    // TODO: change this
    mSourceAmplitude = 1.0f;
}

void FDTD::setupSourceNode( void )
{
    mSourceLocX = mDepthX / 2;
    mSourceLocY = mDepthY / 2;
    mSourceLocZ = mDepthZ / 2;

    mSourceAmplitude = 1.0f;
}

uint32_t FDTD::getXDim( void )
{
    return mDepthX;
}

uint32_t FDTD::getYDim( void )
{
    return mDepthY;
}

uint32_t FDTD::getZDim( void )
{
    return mDepthZ;
}

void FDTD::applySourceNode( void )
{
    //float ret = -2 * ( ( af::Pi * ( mCourantLimit * mCurTimeStep ) / 20 )
    //  * ( af::Pi * ( mCourantLimit * mCurTimeStep ) / 20 ) )
    //    * exp( -( af::Pi * ( mCourantLimit * mCurTimeStep ) / 20 ) );
    //e_x( mDepthX / 2, mDepthY / 2, mDepthZ / 2 ) += ret;
    e_x( mDepthX / 2, mDepthY / 2, mDepthZ / 2 ) += 500.0f * cos( 2*af::Pi*0.05f*mCurTimeStep );

}

// TODO: generalize this out to member geometry
void FDTD::applyArrayNodes( void )
{
    // x is depth, y is far field, z is array direction
    uint32_t arrayCenterX = mDepthX / 2; // keep in middle of vertical
    uint32_t arrayCenterY = mDepthY / 8; // shift it to edge
    uint32_t arrayCenterZ = mDepthZ / 2; // keep lengthwise center of array centered

    float arraySpacing = 0.03; // 30 mm
    float gridSpacing = arraySpacing / mDeltaS;

    double curFreq = 6e9; // 6 GHz
    double curPeriod = 1 / curFreq;
    double f0 = 4e9;
    double f1 = 6e9;
    double T = 15e-9;
    double k = ( f1 - f0 ) / T;
    double timeDelay = 1.5626e-10; // 1 sample time delay
    double constPhase = -1.5;// 2.3617; // 15 deg @ 4 GHz
    double curPhase;
    if( ( mCurTimeStep * mDeltaT ) < T )
    {
        curPhase = 0.03 * cos( 0.2617994 ) * 2 * af::Pi * ( mCurTimeStep * mDeltaT * k ) / 3e8; // 15 deg
        for( int i = -4; i < 4; ++i )
        {
            // TODO replace with actual freq chirp
            e_x( arrayCenterX, arrayCenterY, arrayCenterZ + ( i * gridSpacing ) ) +=
                //500.0f * cos( 2 * af::Pi * ( mCurTimeStep * mDeltaT * curFreq ) );
                //500.0f * cos( 2 * af::Pi * ( f0 * ( mCurTimeStep * mDeltaT + i * timeDelay ) + ( k / 2 )* ( mCurTimeStep * mDeltaT + i * timeDelay ) * ( mCurTimeStep * mDeltaT + i * timeDelay ) ) );
                //500.0f * cos( 2 * af::Pi * ( f0 * ( mCurTimeStep * mDeltaT ) + ( k / 2 )* ( mCurTimeStep * mDeltaT ) * ( mCurTimeStep * mDeltaT ) ) + constPhase * i );
                500.0f * cos( 2 * af::Pi * ( f0 * ( mCurTimeStep * mDeltaT ) + ( k / 2 )* ( mCurTimeStep * mDeltaT ) * ( mCurTimeStep * mDeltaT ) ) + curPhase * i );

        }
    }
    //else
    //{
    //    for( int i = -4; i < 4; ++i )
    //    {
    //        // TODO replace with actual freq chirp
    //        e_x( arrayCenterX, arrayCenterY, arrayCenterZ + ( i * gridSpacing ) ) +=
    //            500.0f * cos( 2 * af::Pi * ( mCurTimeStep * mDeltaT * f1 ) );
    //            //500.0f * cos( 2 * af::Pi * ( f0 * ( mCurTimeStep * mDeltaT + i * timeDelay ) + ( k / 2 )* ( mCurTimeStep * mDeltaT + i * timeDelay ) * ( mCurTimeStep * mDeltaT + i * timeDelay ) ) );
    //    }
    //}

}

double FDTD::getCurrentTime( void )
{
    return ( mCurTimeStep * mDeltaT );
}

void FDTD::enforcePEC( void )
{
    // Back wall
    e_x( af::seq( mDepthX / 2 - 10, mDepthX / 2 + 10 ),
         af::seq( mDepthY / 2 - 10, mDepthY / 2 + 10 ),
         mDepthZ / 2 - 10 ) = 0.0f;
    e_y( af::seq( mDepthX / 2 - 10, mDepthX / 2 + 10 ),
         af::seq( mDepthY / 2 - 10, mDepthY / 2 + 10 ),
         mDepthZ / 2 - 10 ) = 0.0f;
    e_z( af::seq( mDepthX / 2 - 10, mDepthX / 2 + 10 ),
         af::seq( mDepthY / 2 - 10, mDepthY / 2 + 10 ),
         mDepthZ / 2 - 10 ) = 0.0f;

    // Top wall
    e_x( af::seq( mDepthX / 2 - 10, mDepthX / 2 + 10 ),
         mDepthY / 2 + 10,
         af::seq( mDepthZ / 2 - 10, mDepthZ / 2 + 10 ) ) = 0.0f;
    e_y( af::seq( mDepthX / 2 - 10, mDepthX / 2 + 10 ),
         mDepthY / 2 + 10,
         af::seq( mDepthZ / 2 - 10, mDepthZ / 2 + 10 ) ) = 0.0f;
    e_z( af::seq( mDepthX / 2 - 10, mDepthX / 2 + 10 ),
         mDepthY / 2 + 10,
         af::seq( mDepthZ / 2 - 10, mDepthZ / 2 + 10 ) ) = 0.0f;

    // Bottom wall
    e_x( af::seq( mDepthX / 2 - 10, mDepthX / 2 + 10 ),
         mDepthY / 2 - 10,
         af::seq( mDepthZ / 2 - 10, mDepthZ / 2 + 10 ) ) = 0.0f;
    e_y( af::seq( mDepthX / 2 - 10, mDepthX / 2 + 10 ),
         mDepthY / 2 - 10,
         af::seq( mDepthZ / 2 - 10, mDepthZ / 2 + 10 ) ) = 0.0f;
    e_z( af::seq( mDepthX / 2 - 10, mDepthX / 2 + 10 ),
         mDepthY / 2 - 10,
         af::seq( mDepthZ / 2 - 10, mDepthZ / 2 + 10 ) ) = 0.0f;

    // Right wall
    e_x( mDepthX / 2 - 10,
         af::seq( mDepthY / 2 - 10, mDepthY / 2 + 10 ),
         af::seq( mDepthZ / 2 - 10, mDepthZ / 2 + 10 ) ) = 0.0f;
    e_y( mDepthX / 2 - 10,
         af::seq( mDepthY / 2 - 10, mDepthY / 2 + 10 ),
         af::seq( mDepthZ / 2 - 10, mDepthZ / 2 + 10 ) ) = 0.0f;
    e_z( mDepthX / 2 - 10,
         af::seq( mDepthY / 2 - 10, mDepthY / 2 + 10 ),
         af::seq( mDepthZ / 2 - 10, mDepthZ / 2 + 10 ) ) = 0.0f;

    // Left wall
    e_x( mDepthX / 2 + 10,
         af::seq( mDepthY / 2 - 10, mDepthY / 2 + 10 ),
         af::seq( mDepthZ / 2 - 10, mDepthZ / 2 + 10 ) ) = 0.0f;
    e_y( mDepthX / 2 + 10,
         af::seq( mDepthY / 2 - 10, mDepthY / 2 + 10 ),
         af::seq( mDepthZ / 2 - 10, mDepthZ / 2 + 10 ) ) = 0.0f;
    e_z( mDepthX / 2 + 10,
         af::seq( mDepthY / 2 - 10, mDepthY / 2 + 10 ),
         af::seq( mDepthZ / 2 - 10, mDepthZ / 2 + 10 ) ) = 0.0f;
}

void FDTD::timeStep( void )
{
    //enforcePEC();

    updateHField();
    updateEField();

    if( mTimeAvg )
    {
        mAvgPlane = mAvgPlane + ( af::abs( getYZPlane() ) );// / 100.0f );// ( mCurTimeStep * mAvgPlane + af::abs( getYZPlane() ) ) / ( (double) mCurTimeStep + 2 );
    }

    //applySourceNode();
    applyArrayNodes();

    enforceABC();

    mCurTimeStep++;
}

//--------------------------------------------------------------------
// updateEField - Update the E field through one time step
//--------------------------------------------------------------------
void FDTD::updateEField( void )
{
    // x
    e_x( af::seq( 0, mDepthX - 2 ), af::seq( 1, mDepthY - 2 ), af::seq( 1, mDepthZ - 2 ) ) =
        ( coeff_e_xe( af::seq( 0, mDepthX - 2 ), af::seq( 1, mDepthY - 2 ), af::seq( 1, mDepthZ - 2 ) ) 
          * e_x( af::seq( 0, mDepthX - 2 ), af::seq( 1, mDepthY - 2 ), af::seq( 1, mDepthZ - 2 ) ) )
        + coeff_e_xh( af::seq( 0, mDepthX - 2 ), af::seq( 1, mDepthY - 2 ), af::seq( 1, mDepthZ - 2 ) )
        * ( ( h_z( af::seq( 0, mDepthX - 2 ), af::seq( 1, mDepthY - 2 ), af::seq( 1, mDepthZ - 2 ) ) 
            - h_z( af::seq( 0, mDepthX - 2 ), af::seq( 0, mDepthY - 3 ), af::seq( 1, mDepthZ - 2 ) ) )
          - ( h_y( af::seq( 0, mDepthX - 2 ), af::seq( 1, mDepthY - 2 ), af::seq( 1, mDepthZ - 2 ) )
            - h_y( af::seq( 0, mDepthX - 2 ), af::seq( 1, mDepthY - 2 ), af::seq( 0, mDepthZ - 3 ) ) ) );
    // y
    e_y( af::seq( 1, mDepthX - 2 ), af::seq( 0, mDepthY - 2 ), af::seq( 1, mDepthZ - 2 ) ) =
        ( coeff_e_ye( af::seq( 1, mDepthX - 2 ), af::seq( 0, mDepthY - 2 ), af::seq( 1, mDepthZ - 2 ) ) 
          * e_y( af::seq( 1, mDepthX - 2 ), af::seq( 0, mDepthY - 2 ), af::seq( 1, mDepthZ - 2 ) ) )
        + coeff_e_yh( af::seq( 1, mDepthX - 2 ), af::seq( 0, mDepthY - 2 ), af::seq( 1, mDepthZ - 2 ) )
        * ( ( h_x( af::seq( 1, mDepthX - 2 ), af::seq( 0, mDepthY - 2 ), af::seq( 1, mDepthZ - 2 ) )
            - h_x( af::seq( 1, mDepthX - 2 ), af::seq( 0, mDepthY - 2 ), af::seq( 0, mDepthZ - 3 ) ) )
          - ( h_z( af::seq( 1, mDepthX - 2 ), af::seq( 0, mDepthY - 2 ), af::seq( 1, mDepthZ - 2 ) )
            - h_z( af::seq( 0, mDepthX - 3 ), af::seq( 0, mDepthY - 2 ), af::seq( 1, mDepthZ - 2 ) ) ) );
    // z
    e_z( af::seq( 1, mDepthX - 2 ), af::seq( 1, mDepthY - 2 ), af::seq( 0, mDepthZ - 2 ) ) =
        ( coeff_e_ze( af::seq( 1, mDepthX - 2 ), af::seq( 1, mDepthY - 2 ), af::seq( 0, mDepthZ - 2 ) )
          * e_z( af::seq( 1, mDepthX - 2 ), af::seq( 1, mDepthY - 2 ), af::seq( 0, mDepthZ - 2 ) ) )
        + coeff_e_zh( af::seq( 1, mDepthX - 2 ), af::seq( 1, mDepthY - 2 ), af::seq( 0, mDepthZ - 2 ) )
        * ( ( h_y( af::seq( 1, mDepthX - 2 ), af::seq( 1, mDepthY - 2 ), af::seq( 0, mDepthZ - 2 ) )
            - h_y( af::seq( 0, mDepthX - 3 ), af::seq( 1, mDepthY - 2 ), af::seq( 0, mDepthZ - 2 ) ) )
          - ( h_x( af::seq( 1, mDepthX - 2 ), af::seq( 1, mDepthY - 2 ), af::seq( 0, mDepthZ - 2 ) )
            - h_x( af::seq( 1, mDepthX - 2 ), af::seq( 0, mDepthY - 3 ), af::seq( 0, mDepthZ - 2 ) ) ) );
}

//--------------------------------------------------------------------
// updateHField - Update the H field through one time step
//--------------------------------------------------------------------
void FDTD::updateHField( void )
{
    h_x( af::span, af::seq( 0, mDepthY - 2 ), af::seq( 0, mDepthZ - 2 ) ) =
        ( coeff_h_xh( af::span, af::seq( 0, mDepthY - 2 ), af::seq( 0, mDepthZ - 2 ) ) *
          h_x( af::span, af::seq( 0, mDepthY - 2 ), af::seq( 0, mDepthZ - 2 ) ) )
        + coeff_h_xe( af::span, af::seq( 0, mDepthY - 2 ), af::seq( 0, mDepthZ - 2 ) )
        * ( ( e_y( af::span, af::seq( 0, mDepthY - 2 ), af::seq( 1, mDepthZ - 1 ) )
            - e_y( af::span, af::seq( 0, mDepthY - 2 ), af::seq( 0, mDepthZ - 2 ) ) )
            - ( e_z( af::span, af::seq( 1, mDepthY - 1 ), af::seq( 0, mDepthZ - 2 ) )
            - e_z( af::span, af::seq( 0, mDepthY - 2 ), af::seq( 0, mDepthZ - 2 ) ) ) );

    h_y( af::seq( 0, mDepthX - 2 ), af::span, af::seq( 0, mDepthZ - 2 ) ) =
        ( coeff_h_yh( af::seq( 0, mDepthX - 2 ), af::span, af::seq( 0, mDepthZ - 2 ) ) *
          h_y( af::seq( 0, mDepthX - 2 ), af::span, af::seq( 0, mDepthZ - 2 ) ) ) 
        + coeff_h_ye( af::seq( 0, mDepthX - 2 ), af::span, af::seq( 0, mDepthZ - 2 ) )
        * ( ( e_z( af::seq( 1, mDepthX - 1 ), af::span, af::seq( 0, mDepthZ - 2 ) )
            - e_z( af::seq( 0, mDepthX - 2 ), af::span, af::seq( 0, mDepthZ - 2 ) ) ) 
          - ( e_x( af::seq( 0, mDepthX - 2 ), af::span, af::seq( 1, mDepthZ - 1 ) )
            - e_x( af::seq( 0, mDepthX - 2 ), af::span, af::seq( 0, mDepthZ - 2 ) ) ) );

    h_z( af::seq( 0, mDepthX - 2 ), af::seq( 0, mDepthY - 2 ), af::span ) =
        ( coeff_h_zh( af::seq( 0, mDepthX - 2 ), af::seq( 0, mDepthY - 2 ), af::span ) *
          h_z( af::seq( 0, mDepthX - 2 ), af::seq( 0, mDepthY - 2 ), af::span ) ) 
        + coeff_h_ze( af::seq( 0, mDepthX - 2 ), af::seq( 0, mDepthY - 2 ), af::span )
        * ( ( e_x( af::seq( 0, mDepthX - 2 ), af::seq( 1, mDepthY - 1 ), af::span )
            - e_x( af::seq( 0, mDepthX - 2 ), af::seq( 0, mDepthY - 2 ), af::span ) ) 
          - ( e_y( af::seq( 1, mDepthX - 1 ), af::seq( 0, mDepthY - 2 ), af::span )
            - e_y( af::seq( 0, mDepthX - 2 ), af::seq( 0, mDepthY - 2 ), af::span ) ) );
}

//--------------------------------------------------------------------
// enforceABC - Apply a first order ABC to the box
//          TODO: implement second order
//--------------------------------------------------------------------
void FDTD::enforceABC( void )
{
    // X plane 0
    e_y( 0, af::span, af::span ) = abc_xp0_ty + mABCcoeff * ( e_y( 1, af::span, af::span ) - e_y( 0, af::span, af::span ) );
    abc_xp0_ty = e_y( 1, af::span, af::span );
    e_z( 0, af::span, af::span ) = abc_xp0_tz + mABCcoeff * ( e_z( 1, af::span, af::span ) - e_z( 0, af::span, af::span ) );
    abc_xp0_tz = e_z( 1, af::span, af::span );
    // X plane 1
    e_y( mDepthX - 1, af::span, af::span ) = abc_xp1_ty + mABCcoeff * ( e_y( mDepthX - 2, af::span, af::span ) - e_y( mDepthX - 1, af::span, af::span ) );
    abc_xp1_ty = e_y( mDepthX - 2, af::span, af::span );
    e_z( mDepthX - 1, af::span, af::span ) = abc_xp1_tz + mABCcoeff * ( e_z( mDepthX - 2, af::span, af::span ) - e_z( mDepthX - 1, af::span, af::span ) );
    abc_xp1_tz = e_z( mDepthX - 2, af::span, af::span );
    // Y Plane 0
    e_x( af::span, 0, af::span ) = abc_yp0_tx + mABCcoeff * ( e_x( af::span, 1, af::span ) - e_x( af::span, 0, af::span ) );
    abc_yp0_tx = e_x( af::span, 1, af::span );
    e_z( af::span, 0, af::span ) = abc_yp0_tz + mABCcoeff * ( e_z( af::span, 1, af::span ) - e_z( af::span, 0, af::span ) );
    abc_yp0_tz = e_z( af::span, 1, af::span );
    // Y Plane 1
    e_x( af::span, mDepthY - 1, af::span ) = abc_yp1_tx + mABCcoeff * ( e_x( af::span, mDepthY - 2, af::span ) - e_x( af::span, mDepthY - 1, af::span ) );
    abc_yp1_tx = e_x( af::span, mDepthY - 2, af::span );
    e_z( af::span, mDepthY - 1, af::span ) = abc_yp1_tz + mABCcoeff * ( e_z( af::span, mDepthY - 2, af::span ) - e_z( af::span, mDepthY - 1, af::span ) );
    abc_yp1_tz = e_z( af::span, mDepthY - 2, af::span );
    // Z Plane 0
    e_x( af::span, af::span, 0 ) = abc_zp0_tx + mABCcoeff * ( e_x( af::span, af::span, 1 ) - e_x( af::span, af::span, 0 ) );
    abc_zp0_tx = e_x( af::span, af::span, 1 );
    e_y( af::span, af::span, 0 ) = abc_zp0_ty + mABCcoeff * ( e_y( af::span, af::span, 1 ) - e_y( af::span, af::span, 0 ) );
    abc_zp0_ty = e_y( af::span, af::span, 1 );
    // Z Plane 1
    e_x( af::span, af::span, mDepthZ - 1 ) = abc_zp1_tx + mABCcoeff * ( e_x( af::span, af::span, mDepthZ - 2 ) - e_x( af::span, af::span, mDepthZ - 1 ) );
    abc_zp1_tx = e_x( af::span, af::span, mDepthZ - 2 );
    e_y( af::span, af::span, mDepthZ - 1 ) = abc_zp1_ty + mABCcoeff * ( e_y( af::span, af::span, mDepthZ - 2 ) - e_y( af::span, af::span, mDepthZ - 1 ) );
    abc_zp1_ty = e_y( af::span, af::span, mDepthZ - 2 );

}

af::array FDTD::getYZPlane( uint32_t x )
{
    return e_x( x, af::span, af::span );// +e_x( x - 1, af::span, af::span ) + e_x( x + 1, af::span, af::span );// +e_z( x, af::span, af::span ) + e_x( x, af::span, af::span );
}

af::array FDTD::getYZPlane( void )
{
    return getYZPlane( mDepthX / 2 );
}

af::array FDTD::getAvgPlane( void )
{
    return mAvgPlane;
}

void FDTD::setTimeAvg( void )
{
    mTimeAvg = true;
}

//// set material permittivity arbitrarily
//// by updating coefficients
//// In order to match across boundaries, sigma_m / mu needs to = sigma / epsilon
//for( i = depth / 2; i < 0.75*depth; ++i )
//{
//    coeff_Ez_e[i] = ( 1.0 - d_loss ) / ( 1.0 + d_loss );
//    coeff_Ez_h[i] = ( fsi / 9.0 ) / ( 1.0 + d_loss );

//    coeff_Hy_h[i] = ( 1.0 - d_loss ) / ( 1.0 + d_loss );
//    coeff_Hy_e[i] = ( 1.0 / fsi ) / ( 1.0 + d_loss );
//}


void FDTD::add2dAntenna( Shape2D outline, float thickness, Vec3 centerPoint, Material mat )
{

}

