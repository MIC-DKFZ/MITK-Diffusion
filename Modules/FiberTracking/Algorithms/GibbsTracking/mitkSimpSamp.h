/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef _SIMPSAMP
#define _SIMPSAMP

#include <MitkFiberTrackingExports.h>
#include <mitkParticle.h>
#include <algorithm>
#include <vector>

namespace mitk
{

/**
* rief Samples new tract from surrounding fiber segments: draws one of the added end points with
* probability proportional to its weight. The buffers grow with the number of neighbours (up to
* two end points per particle of the eight neighbouring grid cells) and keep their capacity
* across clear(). */

class MITKFIBERTRACKING_EXPORT SimpSamp
{
    std::vector<float> P;   ///< P[0] = 0, P[i] = sum of the first i weights

public:
    std::vector<EndPoint> objs;

    SimpSamp()
    {
        clear();
    }

    inline void clear()
    {
        P.assign(1, 0.0f);
        objs.clear();
    }

    inline void add(float p, EndPoint obj)
    {
        P.push_back(P.back() + p);
        objs.push_back(obj);
    }

    /** Index of the end point whose interval [P[i], P[i+1]) contains prob*total; prob in [0,1]. */
    inline int draw(float prob)
    {
        const float r = prob*P.back();
        const int bin = static_cast<int>(std::upper_bound(P.begin()+1, P.end(), r) - P.begin()) - 1;
        return std::min(bin, static_cast<int>(objs.size())-1);  // r at the total: last end point
    }

    float probFor(int idx)
    {
        return (P[idx+1]-P[idx])/P.back();
    }

    float probFor(EndPoint& t)
    {
        for (std::size_t i = 0; i < objs.size(); i++)
        {
            if (t == objs[i])
                return probFor(static_cast<int>(i));
        }
        return 0;
    }
};

}

#endif
