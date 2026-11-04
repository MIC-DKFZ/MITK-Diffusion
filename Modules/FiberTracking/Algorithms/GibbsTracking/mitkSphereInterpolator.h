/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef _SPHEREINTERPOLATOR
#define _SPHEREINTERPOLATOR

#include <MitkFiberTrackingExports.h>
#include <vnl/vnl_vector_fixed.h>
#include <vector>

/**
 * \brief Lookuptable based trilinear interpolation of spherically arranged scalar values.
 */
class MITKFIBERTRACKING_EXPORT SphereInterpolator
{

public:

    int   size;     // size of LUT
    int   sN;       // (sizeofLUT-1)/2
    int   nverts;   // number of data vertices
    float beta;
    float inva;
    float b;
    bool    m_ValidState;

    std::vector< float > barycoords;
    std::vector< int >   indices;

    /** The three ODF directions around a direction (1-based indices) and their weights. Returned by
     *  value, so that threads of the parallel sampler can share one interpolator. */
    struct Interpolation
    {
        vnl_vector_fixed< int, 3 >     idx;
        vnl_vector_fixed< float, 3 >   interpw;
    };

    SphereInterpolator(const std::string& lutPath);

    inline bool IsInValidState() const
    {
        return m_ValidState;
    }

    ~SphereInterpolator();

    inline Interpolation getInterpolation(const vnl_vector_fixed<float, 3>& N) const
    {
        float nx = N[0];
        float ny = N[1];
        float nz = N[2];

        if (nz > 0.5)
        {
            int x = float2int(nx);
            int y = float2int(ny);
            int i = 3*6*(x+y*size);  // (:,1,x,y)
            return lookup(i);
        }
        if (nz < -0.5)
        {
            int x = float2int(nx);
            int y = float2int(ny);
            int i = 3*(1+6*(x+y*size));  // (:,2,x,y)
            return lookup(i);
        }
        if (nx > 0.5)
        {
            int z = float2int(nz);
            int y = float2int(ny);
            int i = 3*(2+6*(z+y*size));  // (:,2,x,y)
            return lookup(i);
        }
        if (nx < -0.5)
        {
            int z = float2int(nz);
            int y = float2int(ny);
            int i = 3*(3+6*(z+y*size));  // (:,2,x,y)
            return lookup(i);
        }
        if (ny > 0)
        {
            int x = float2int(nx);
            int z = float2int(nz);
            int i = 3*(4+6*(x+z*size));  // (:,1,x,y)
            return lookup(i);
        }
        else
        {
            int x = float2int(nx);
            int z = float2int(nz);
            int i = 3*(5+6*(x+z*size));  // (:,1,x,y)
            return lookup(i);
        }
    }

protected:

    inline Interpolation lookup(int i) const
    {
        Interpolation result;
        result.idx[0] = indices[i];
        result.idx[1] = indices[i+1];
        result.idx[2] = indices[i+2];
        result.interpw[0] = barycoords[i];
        result.interpw[1] = barycoords[i+1];
        result.interpw[2] = barycoords[i+2];
        return result;
    }

    bool LoadLookuptables(const std::string& lutPath);

    bool LoadLookuptables();

    bool LoadLookuptables(std::istream& BaryCoordsStream, std::istream& IndicesStream);

    inline float invrescale(float f) const
    {
        float x = (fabs(f)-b)*inva;
        if (f>0)
            return (x*x-beta);
        else
            return beta - x*x;
    }

    inline int float2int(float x) const
    {
        return int((invrescale(x)+1)*sN-0.5);
    }

};

#endif
