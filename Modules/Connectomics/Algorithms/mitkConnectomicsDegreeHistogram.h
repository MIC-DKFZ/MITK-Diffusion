
/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/


#ifndef _MITK_ConnectomicsDegreeHistogram_H
#define _MITK_ConnectomicsDegreeHistogram_H

#include<mitkConnectomicsHistogramBase.h>

namespace mitk {

    /**
  * \brief A class to calculate and store the degree of each node   */
  class ConnectomicsDegreeHistogram : public mitk::ConnectomicsHistogramBase
  {

  public:

    ConnectomicsDegreeHistogram();
    ~ConnectomicsDegreeHistogram() override;



  protected:

    /** @brief Creates a new histogram from the network source. */
    void ComputeFromConnectomicsNetwork( ConnectomicsNetwork* source ) override;

  };

}

#endif /* _MITK_ConnectomicsDegreeHistogram_H */
