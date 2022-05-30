/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#include "mitkStreamlineFeatureExtractor.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <boost/progress.hpp>
#include <vnl/vnl_sparse_matrix.h>
#include <mitkIOUtil.h>

namespace mitk{

StreamlineFeatureExtractor::StreamlineFeatureExtractor()
  : m_NumPoints(20)
{

}

StreamlineFeatureExtractor::~StreamlineFeatureExtractor()
{

}



void StreamlineFeatureExtractor::SetTractogramPlus(const mitk::FiberBundle::Pointer &TractogramPlus)
{
  m_TractogramPlus = TractogramPlus;
}


void StreamlineFeatureExtractor::SetTractogramMinus(const mitk::FiberBundle::Pointer &TractogramMinus)
{
  m_TractogramMinus = TractogramMinus;
}

void StreamlineFeatureExtractor::SetTractogramTest(const mitk::FiberBundle::Pointer &TractogramTest)
{
  m_TractogramTest = TractogramTest;
}


std::vector<vnl_matrix<float> > StreamlineFeatureExtractor::ResampleFibers(mitk::FiberBundle::Pointer tractogram)
{
  mitk::FiberBundle::Pointer temp_fib = tractogram->GetDeepCopy();
  temp_fib->ResampleToNumPoints(m_NumPoints);

  std::vector< vnl_matrix<float> > out_fib;

  for (int i=0; i<temp_fib->GetFiberPolyData()->GetNumberOfCells(); i++)
  {
    vtkCell* cell = temp_fib->GetFiberPolyData()->GetCell(i);
    int numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    vnl_matrix<float> streamline;
    streamline.set_size(3, m_NumPoints);
    streamline.fill(0.0);

    for (int j=0; j<numPoints; j++)
    {
      double cand[3];
      points->GetPoint(j, cand);

      vnl_vector_fixed< float, 3 > candV;
      candV[0]=cand[0]; candV[1]=cand[1]; candV[2]=cand[2];
      streamline.set_column(j, candV);
    }

    out_fib.push_back(streamline);
  }

  return out_fib;
}

std::vector<vnl_matrix<float> > StreamlineFeatureExtractor::CalculateDmdf(std::vector<vnl_matrix<float> > tractogram, std::vector<vnl_matrix<float> > prototypes)
{
    std::vector< vnl_matrix<float> >  dist_vec;


    for (unsigned int i=0; i<tractogram.size(); i++)
    {

        vnl_matrix<float> distances;
        distances.set_size(1, prototypes.size());
        distances.fill(0.0);


        for (unsigned int j=0; j<prototypes.size(); j++)
        {
            vnl_matrix<float> single_distances;
            single_distances.set_size(1, tractogram.at(0).cols());
            single_distances.fill(0.0);
            vnl_matrix<float> single_distances_flip;
            single_distances_flip.set_size(1, tractogram.at(0).cols());
            single_distances_flip.fill(0.0);
            for (unsigned int ik=0; ik<tractogram.at(0).cols(); ik++)
            {
                double cur_dist;
                double cur_dist_flip;

                cur_dist = sqrt(pow(tractogram.at(i).get(0,ik) - prototypes.at(j).get(0,ik), 2.0) +
                                       pow(tractogram.at(i).get(1,ik) - prototypes.at(j).get(1,ik), 2.0) +
                                       pow(tractogram.at(i).get(2,ik) - prototypes.at(j).get(2,ik), 2.0));

                cur_dist_flip = sqrt(pow(tractogram.at(i).get(0,ik) - prototypes.at(j).get(0,prototypes.at(0).cols()-(ik+1)), 2.0) +
                                       pow(tractogram.at(i).get(1,ik) - prototypes.at(j).get(1,prototypes.at(0).cols()-(ik+1)), 2.0) +
                                       pow(tractogram.at(i).get(2,ik) - prototypes.at(j).get(2,prototypes.at(0).cols()-(ik+1)), 2.0));

                single_distances.put(0,ik, cur_dist);
                single_distances_flip.put(0,ik, cur_dist_flip);

            }

            if (single_distances_flip.mean()> single_distances.mean())
            {
                distances.put(0,j, single_distances.mean());
            }
            else {
                distances.put(0,j, single_distances_flip.mean());
            }



        }
        dist_vec.push_back(distances);
    }

    return dist_vec;
}

void StreamlineFeatureExtractor::GetData()
{
    MITK_INFO << "Start Function Ged Data";
    int labels_arr [m_DistancesPlus.size()+m_DistancesMinus.size()];
    float data_arr [m_DistancesPlus.size()+m_DistancesMinus.size()][m_DistancesPlus.at(0).size()];


    int size_plus = 0;
    for ( unsigned int i=0; i<m_DistancesPlus.size(); i++)
    {
        MITK_INFO << "Plus: Print i";
        MITK_INFO << i;
        MITK_INFO << "Adding 1";
        labels_arr[i]=1;
        for ( unsigned int j=0; j<m_DistancesPlus.at(0).cols(); j++)
        {
            data_arr[i][j] = m_DistancesPlus.at(i).get(0,j);
        }
        size_plus++;
    }
    MITK_INFO << "Positive Labels added";
    for ( unsigned int i=m_DistancesPlus.size(); i<m_DistancesPlus.size()+m_DistancesMinus.size(); i++)
    {
        int it = i - size_plus;
        MITK_INFO << "Minus: Print i";
        MITK_INFO << it;
        labels_arr[i]=0;
        MITK_INFO << "Adding 0";
        for ( unsigned int j=0; j<m_DistancesMinus.at(0).cols(); j++)
        {
            data_arr[i][j] = m_DistancesMinus.at(it).get(0,j);
        }
    }
    MITK_INFO << "Negavtive Labels added";



    cv::Mat labels(m_DistancesPlus.size()+m_DistancesMinus.size(), 1, CV_32S, labels_arr);
    cv::Mat data(m_DistancesPlus.size()+m_DistancesMinus.size(), m_DistancesPlus.at(0).size(), CV_32F, data_arr);

    MITK_INFO << "data";
    MITK_INFO << data;
    MITK_INFO << "labels";
    MITK_INFO << labels;
    cv::Ptr<cv::ml::TrainData> m_traindata = cv::ml::TrainData::create(data, cv::ml::ROW_SAMPLE, labels);

    cv::Ptr<cv::ml::RTrees> m_rtrees;
    m_rtrees = cv::ml::RTrees::create();
    m_rtrees->setMaxDepth(3000);
    m_rtrees->setMinSampleCount(1);
    m_rtrees->setUse1SERule(false);
    m_rtrees->setUseSurrogates(false);
    m_rtrees->setPriors(cv::Mat());
    m_rtrees->setCVFolds(1);

    MITK_INFO << "Start Training";
    m_rtrees->train(m_traindata);

    MITK_INFO << "Predicting";
//    std::vector<std::vector<float>> test_arr(m_DistancesTest.size(), std::vector<float> (m_DistancesPlus.at(0).size(), 0));
//    MITK_INFO << "test_arr";
//    for ( unsigned int i=m_DistancesTest.size(); i<m_DistancesTest.size(); i++)
//    {
//        for ( unsigned int j=0; j<m_DistancesTest.at(0).cols(); j++)
//        {
//            test_arr[i][j] = m_DistancesTest.at(i).get(0,j);
//        }
//    }
//    cv::Mat dataTest(m_DistancesTest.size(), m_DistancesPlus.at(0).size(), CV_32F);
//        for ( unsigned int i=m_DistancesTest.size(); i<m_DistancesTest.size(); i++)
//        {
//               for ( unsigned int j=0; j<m_DistancesTest.at(0).cols(); j++)
//                    {
//                        dataTest[i][j] = m_DistancesTest.at(i).get(0,j);
//                    }
//        }

//    int val;
//    for (int i = 0; i < 20; i++)
//    {
//        //val = m_rtrees->predict(testingMat.row(i), predictLabels);
//        val = std::lround(m_rtrees->predict(dataTest.row(i)));
//        MITK_INFO << val;
//    }





}

//void StreamlineFeatureExtractor::CreateClassifier()
//{
//    MITK_INFO << m_traindata->getResponses();


//}


void StreamlineFeatureExtractor::GenerateData()
{
    MITK_INFO << "Update";
  mitk::FiberBundle::Pointer inputPrototypes = mitk::IOUtil::Load<mitk::FiberBundle>("/home/r948e/E132-Projekte/Projects/2022_Peretzke_Interactive_Fiber_Dissection/mitk_diff/prototypes_599671.trk");

  T_Prototypes = ResampleFibers(inputPrototypes);
  T_TractogramMinus= ResampleFibers(m_TractogramMinus);
  T_TractogramPlus= ResampleFibers(m_TractogramPlus);
  T_TractogramTest= ResampleFibers(m_TractogramTest);

  MITK_INFO << "Calculate Features...";
  m_DistancesPlus = CalculateDmdf(T_TractogramPlus, T_Prototypes);
  m_DistancesMinus = CalculateDmdf(T_TractogramMinus, T_Prototypes);
  m_DistancesTest= CalculateDmdf(T_TractogramTest, T_Prototypes);

  MITK_INFO << "Sizes of Plus and Minus";
  MITK_INFO << m_DistancesPlus.size() + m_DistancesMinus.size();
  MITK_INFO << "Size of Test Data";
  MITK_INFO << m_DistancesTest.size();
  MITK_INFO << "Done with Datacreation";

}


void StreamlineFeatureExtractor::StartAlgorithm()
{
    MITK_INFO << "Printing";

//    MITK_INFO << dataset->getResponses();
}



}




