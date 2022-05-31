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
    cv::Mat data;

    int size_plus = 0;



    for ( unsigned int i=0; i<m_DistancesPlus.size(); i++)
    {
        float data_arr [m_DistancesPlus.at(0).size()];
        MITK_INFO << "Plus: Print i";
        MITK_INFO << i;
        MITK_INFO << "Adding 1";
        labels_arr[i]=1;

        for ( unsigned int j=0; j<m_DistancesPlus.at(0).cols(); j++)
        {
            data_arr[j] = m_DistancesPlus.at(i).get(0,j);
        }
        cv::Mat curdata(1, m_DistancesPlus.at(0).size(), CV_32F, data_arr);
        MITK_INFO << curdata;
        data.push_back(curdata);
        size_plus++;
    }

    for ( unsigned int i=m_DistancesPlus.size(); i<m_DistancesPlus.size()+m_DistancesMinus.size(); i++)
    {
        int it = i - size_plus;
        float data_arr [m_DistancesMinus.at(0).size()];
        MITK_INFO << "Plus: Print i";
        MITK_INFO << i;
        MITK_INFO << "Adding 0";
        labels_arr[i]=0;

        for ( unsigned int j=0; j<m_DistancesMinus.at(0).cols(); j++)
        {
            data_arr[j] = m_DistancesMinus.at(it).get(0,j);
        }
        cv::Mat curdata(1, m_DistancesPlus.at(0).size(), CV_32F, data_arr);
        data.push_back(curdata);
        size_plus++;
    }

    cv::Mat labels(m_DistancesPlus.size()+m_DistancesMinus.size(), 1, CV_32S, labels_arr);


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


    cv::Mat dataTest;
    for ( unsigned int i=0; i<m_DistancesTest.size(); i++)
    {
        float data_arr [m_DistancesTest.at(0).size()];

        for ( unsigned int j=0; j<m_DistancesTest.at(0).cols(); j++)
        {
            data_arr[j] = m_DistancesTest.at(i).get(0,j);
        }
        cv::Mat curdata(1, m_DistancesTest.at(0).size(), CV_32F, data_arr);
        dataTest.push_back(curdata);
    }


//      parallel_for<size_t> (size_t(0), m_DistancesTest.size(), [&](size_t i))
//      {
//         float data_arr [m_DistancesTest.at(0).size()];
//         for (int j = 0; j < size; j++)
//         {
//            data_arr[j] = m_DistancesTest.at(i).get(0,j);
//         }
//         cv::Mat curdata(1, m_DistancesTest.at(0).size(), CV_32F, data_arr);
//         MITK_INFO << curdata;
//         dataTest.push_back(curdata);

//       });



    int val;
    std::vector<int> preds;
    for (unsigned int i = 0; i < m_DistancesTest.size(); i++)
    {
        //val = m_rtrees->predict(testingMat.row(i), predictLabels);
        val = std::lround(m_rtrees->predict(dataTest.row(i)));
        preds.push_back(val);
    }

    MITK_INFO << preds.size();

    auto it = std::find(preds.begin(), preds.end(), 1);

    // If element was found
    if (it != preds.end())
    {

        // calculating the index
        // of K
        int index = it - preds.begin();
        cout << index << endl;
    }
    else {
        // If the element is not
        // present in the vector
        cout << "-1" << endl;
    }





}





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


cv::Mat StreamlineFeatureExtractor::StartAlgorithm()
{
    MITK_INFO << "Printing";
    float data_arr [10] = {1, 2.4 ,4 ,4.5 ,6 ,7, 120.5, 100, 120, 100};
    cv::Mat curdata(1, 10, CV_32F, data_arr);
    float data_arr2 [10] = {10, 20.4 ,40 ,40.5 ,60 ,70, 1200.5, 1000, 1200, 1000};
    cv::Mat curdata2(1, 10, CV_32F, data_arr2);

    cv::Mat data;
//    cv::Mat data2;
    //    data.row(1) = curdata.clone();
    data.push_back(curdata);
    data.push_back(curdata2);
//    cv::add(curdata,data2,data2);
    cout << curdata;
    cout << data;
//    cout << data2;
    return curdata.clone();
}



}




