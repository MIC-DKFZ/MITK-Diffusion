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

std::vector<unsigned int>  StreamlineFeatureExtractor::GetData()
{
    MITK_INFO << "Start Function Ged Data";
    int labels_arr [m_DistancesPlus.size()+m_DistancesMinus.size()];
    cv::Mat data;

    int size_plus = 0;



    for ( unsigned int i=0; i<m_DistancesPlus.size(); i++)
    {
        float data_arr [m_DistancesPlus.at(0).size()];
        labels_arr[i]=1;

        for ( unsigned int j=0; j<m_DistancesPlus.at(0).cols(); j++)
        {
            data_arr[j] = m_DistancesPlus.at(i).get(0,j);
        }
        cv::Mat curdata(1, m_DistancesPlus.at(0).size(), CV_32F, data_arr);
        data.push_back(curdata);
        size_plus++;
    }

    for ( unsigned int i=m_DistancesPlus.size(); i<m_DistancesPlus.size()+m_DistancesMinus.size(); i++)
    {
        int it = i - size_plus;
        float data_arr [m_DistancesMinus.at(0).size()];
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
    MITK_INFO << labels;
    std::fstream outputFile;
    outputFile.open( "/home/r948e/labels.csv", std::ios::out ) ;

    for(int i=0; i<labels.rows; i++)
    {
        for(int j=0; j<labels.cols; j++)
        {
            outputFile << labels.at<int>(i,j) << ", ";
        }
        outputFile << endl;

    }


    std::fstream outputFile2;
    outputFile2.open( "/home/r948e/samples.csv", std::ios::out ) ;

    for(int i=0; i<data.rows; i++)
    {
        for(int j=0; j<data.cols; j++)
        {
            outputFile2 << data.at<float>(i,j) << ", ";
        }
        outputFile2 << endl;

    }


//    MITK_INFO << "data";
//    MITK_INFO << data;
    MITK_INFO << "labels";
    MITK_INFO << labels;

    cv::Ptr<cv::ml::TrainData> m_traindata = cv::ml::TrainData::create(data, cv::ml::ROW_SAMPLE, labels);
    m_traindata->shuffleTrainTest();
    MITK_INFO << m_traindata->getResponses();


    MITK_INFO << "Start Training";

    auto random_forest = cv::ml::RTrees::create();

//    auto criteria = cv::TermCriteria();
//    criteria.type = cv::TermCriteria::EPS + cv::TermCriteria::MAX_ITER;
//    criteria.epsilon = 1e-8;
//    criteria.maxCount = 5000;

    random_forest->setMaxCategories(2);
    random_forest->setMaxDepth(50);
    random_forest->setMinSampleCount(1);
    random_forest->setTruncatePrunedTree(true);
    random_forest->setUse1SERule(true);
    random_forest->setUseSurrogates(false);
    random_forest->setTermCriteria(cv::TermCriteria(5000, 800, 1e-6));
    random_forest->setCVFolds(1);

    random_forest->train(m_traindata);


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


    std::fstream outputFile3;
    outputFile3.open( "/home/r948e/test.csv", std::ios::out ) ;

    for(int i=0; i<dataTest.rows; i++)
    {
        for(int j=0; j<dataTest.cols; j++)
        {
            outputFile3 << dataTest.at<float>(i,j) << ", ";
        }
        outputFile3 << endl;

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



//    cv::Mat pred;

    int val;
    std::vector<int> pred;
    std::vector<unsigned int> index;
    for (unsigned int i = 0; i < m_DistancesTest.size(); i++)
    {
        val = random_forest->predict(dataTest.row(i));
        pred.push_back(val);
        MITK_INFO << val;

        //    random_forest->getVotes(dataTest, );
        if (val==1)
        {
           index.push_back(i);
        }
    }



    return index;

}

void StreamlineFeatureExtractor::CreatePrediction(std::vector<unsigned int> &index)
{




    MITK_INFO << "Create Bundle";

    vtkSmartPointer<vtkPolyData> FibersData;
    FibersData = vtkSmartPointer<vtkPolyData>::New();
    FibersData->SetPoints(vtkSmartPointer<vtkPoints>::New());
    FibersData->SetLines(vtkSmartPointer<vtkCellArray>::New());

    vtkSmartPointer<vtkPolyData> vNewPolyData = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkCellArray> vNewLines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> vNewPoints = vtkSmartPointer<vtkPoints>::New();

    unsigned int indexSize = index.size();
    unsigned int counter = 0;
    MITK_INFO << "Start Loop";
    for (unsigned int i=0; i<indexSize; i++)
    {

        vtkCell* cell = m_TractogramTest->GetFiberPolyData()->GetCell(index[i]);
        auto numPoints = cell->GetNumberOfPoints();
        vtkPoints* points = cell->GetPoints();

        vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
        for (unsigned int j=0; j<numPoints; j++)
        {
            double p[3];
            points->GetPoint(j, p);

            vtkIdType id = vNewPoints->InsertNextPoint(p);
            container->GetPointIds()->InsertNextId(id);
        }
        //    weights->InsertValue(counter, fib->GetFiberWeight(i));
        vNewLines->InsertNextCell(container);
        counter++;

    }
    MITK_INFO << "Counter";
    MITK_INFO << counter;


     vNewPolyData->SetLines(vNewLines);
     vNewPolyData->SetPoints(vNewPoints);

     FibersData = vtkSmartPointer<vtkPolyData>::New();
     FibersData->SetPoints(vtkSmartPointer<vtkPoints>::New());
     FibersData->SetLines(vtkSmartPointer<vtkCellArray>::New());
     FibersData->SetPoints(vNewPoints);
     FibersData->SetLines(vNewLines);

     m_Prediction = mitk::FiberBundle::New(vNewPolyData);

//     Bundle->SetFiberColors(255, 255, 255);
    MITK_INFO << "Cells Prediciton";
    MITK_INFO << m_Prediction->GetFiberPolyData()->GetNumberOfCells();
    MITK_INFO << "Cells Tractorgram";
    MITK_INFO << m_TractogramTest->GetFiberPolyData()->GetNumberOfCells();
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
  std::vector<unsigned int> index =GetData();
  CreatePrediction (index);

}


//cv::Mat StreamlineFeatureExtractor::StartAlgorithm()
//{
//    MITK_INFO << "Printing";
//    float data_arr [10] = {1, 2.4 ,4 ,4.5 ,6 ,7, 120.5, 100, 120, 100};
//    cv::Mat curdata(1, 10, CV_32F, data_arr);
//    float data_arr2 [10] = {10, 20.4 ,40 ,40.5 ,60 ,70, 1200.5, 1000, 1200, 1000};
//    cv::Mat curdata2(1, 10, CV_32F, data_arr2);

//    cv::Mat data;
////    cv::Mat data2;
//    //    data.row(1) = curdata.clone();
//    data.push_back(curdata);
//    data.push_back(curdata2);
////    cv::add(curdata,data2,data2);
//    cout << curdata;
//    cout << data;
////    cout << data2;
//    return curdata.clone();
//}



}




