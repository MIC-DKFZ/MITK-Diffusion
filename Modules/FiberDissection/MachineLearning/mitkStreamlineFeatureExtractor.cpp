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


void StreamlineFeatureExtractor::SetTractogramTest(const mitk::FiberBundle::Pointer &TractogramTest, std::string TractogramTestName)
{
    std::string path = "/home/r948e/E132-Projekte/Projects/2022_Peretzke_Interactive_Fiber_Dissection/mitk_diff/storage/";
    path.append(TractogramTestName);
    m_TractogramTest= TractogramTest;
    m_DistancesTestName= path.append("_distances.csv");
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

std::vector<std::vector<unsigned int>>  StreamlineFeatureExtractor::GetData()
{
    MITK_INFO << "Start Function Get Data";
    std::vector<std::vector<unsigned int>> index_vec;
    float labels_arr [m_DistancesPlus.size()+m_DistancesMinus.size()];
    cv::Mat data;
    cv::Mat labels_arr_vec;

    int size_plus = 0;



    for ( unsigned int i=0; i<m_DistancesPlus.size(); i++)
    {
        float data_arr [m_DistancesPlus.at(0).size()];
        labels_arr[i]=1;
        labels_arr_vec.push_back(1);

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
        labels_arr_vec.push_back(0);

         for ( unsigned int j=0; j<m_DistancesMinus.at(0).cols(); j++)
        {
            data_arr[j] = m_DistancesMinus.at(it).get(0,j);
        }
        cv::Mat curdata(1, m_DistancesPlus.at(0).size(), CV_32F, data_arr);
        data.push_back(curdata);
        size_plus++;
    }

    cv::Mat weights;
    int zerosgt = labels_arr_vec.rows - cv::countNonZero(labels_arr_vec);
    int onesgt = cv::countNonZero(labels_arr_vec);
    float plusval = labels_arr_vec.rows / (2.0 * onesgt );
    float minusval = labels_arr_vec.rows / (2.0 * zerosgt );

    // Create sample weights
//    for (int i=0; i<labels_arr_vec.rows; i++ )
//    {
//        if (i<onesgt)
//        {
//            weights.push_back(plusval);
//        }
//        else {
//            weights.push_back(minusval);
//        }
//    }
    weights.push_back(minusval);
    weights.push_back(plusval);
    MITK_INFO << "Weights";

    cv::Mat newweight;

    newweight.push_back(minusval);
    newweight.push_back(plusval);
    MITK_INFO << "Weights";
    MITK_INFO << newweight;


//    cv::Mat labels(m_DistancesPlus.size()+m_DistancesMinus.size(), 1, CV_32S, labels_arr);
    cv::Mat labels(m_DistancesPlus.size()+m_DistancesMinus.size(), 1, CV_32F, labels_arr);

    //    }


//    MITK_INFO << "data";
//    MITK_INFO << data;
//    MITK_INFO << "labels";
//    MITK_INFO << labels;

    cv::Ptr<cv::ml::TrainData> m_traindata = cv::ml::TrainData::create(data, cv::ml::ROW_SAMPLE, labels);
//    m_traindata->setTrainTestSplitRatio(1, true);
    m_traindata->shuffleTrainTest();


    MITK_INFO << "Start Training";

    auto statistic_model = cv::ml::RTrees::create();

////    auto criteria = cv::TermCriteria();
////    criteria.type = cv::TermCriteria::EPS + cv::TermCriteria::MAX_ITER;
////    criteria.epsilon = 1e-8;
////    criteria.maxCount = 5000;

    statistic_model->setMaxCategories(80);
    statistic_model->setMaxDepth(50);
    statistic_model->setMinSampleCount(3);
    statistic_model->setTruncatePrunedTree(true);
    statistic_model->setUse1SERule(true);
    statistic_model->setUseSurrogates(false);
    statistic_model->setTermCriteria(cv::TermCriteria(1500, 1000, 1e-8));
    statistic_model->setCVFolds(1);
    statistic_model->setPriors(newweight);

    statistic_model->train(m_traindata);

//    auto logistic_regression = cv::ml::LogisticRegression::create();

//    statistic_model->setLearningRate(0.001);
//    statistic_model->setIterations(100);
//    statistic_model->setRegularization(cv::ml::LogisticRegression::REG_L2);
//    statistic_model->setTrainMethod(cv::ml::LogisticRegression::MINI_BATCH);
//    statistic_model->setMiniBatchSize(100);


//    statistic_model->train(m_traindata);



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
//    MITK_INFO << "DataTest";
//    MITK_INFO << dataTest.cols;
//    MITK_INFO << dataTest.rows;
//    MITK_INFO << dataTest;
    MITK_INFO << "________";







    int val;
    std::vector<int> pred(m_DistancesTest.size());
    std::vector<unsigned int> index(m_DistancesTest.size());
    float one_prob;
    float zero_prob;
    int ones;
    int zeros;



    std::vector<float> e(m_DistancesTest.size()) ;




    cv::parallel_for_(cv::Range(0, m_DistancesTest.size()), [&](const cv::Range &range)
    {
    for (int i = range.start; i < range.end; i++)
    {


        cv::Mat vote;
        val = statistic_model->predict(dataTest.row(i));
        statistic_model->getVotes(dataTest.row(i), vote, 0);

        ones  = cv::countNonZero(vote);
        zeros = vote.cols - cv::countNonZero(vote);

        one_prob = ones/ (vote.cols * 1.0);
        zero_prob = zeros / (vote.cols * 1.0);

        e.at(i) = ( -one_prob * log2(one_prob) - zero_prob * log2(zero_prob));



        pred.at(i) = val;


        if (val==1)
        {
           index.at(i) = i;
        }

    }
    });
    index.erase(
        std::remove(index.begin(), index.end(), 0),
        index.end());
    index.shrink_to_fit();


    MITK_INFO << "--------------";
    e.erase(
        std::remove(e.begin(), e.end(), 0),
        e.end());
    e.shrink_to_fit();
    auto it = std::minmax_element(e.begin(), e.end());
    int min_idx = std::distance(e.begin(), it.first);
    int max_idx = std::distance(e.begin(), it.second);
    std::cout << min_idx << ", " << max_idx << std::endl; // 1, 5

    MITK_INFO << "--------------";

    MITK_INFO << "Start the ordering";
    std::vector<unsigned int> indextolabel;
    std::priority_queue<std::pair<float, int>> q;
      for (unsigned int i = 0; i < e.size(); ++i) {
        q.push(std::pair<float, int>(e[i], i));
      }
      int k = m_DistancesTest.size(); // number of indices we need
//      int k = 500; // number of indices we need
      for (int i = 0; i < k; ++i) {
        int ki = q.top().second;
//        std::cout << "index[" << i << "] = " << ki << std::endl;x
        indextolabel.push_back(ki);
        q.pop();
      }
      MITK_INFO << "Done";

//    MITK_INFO << statistic_model->getTermCriteria().maxCount;
//    MITK_INFO << statistic_model->getTermCriteria().type;
//    MITK_INFO << statistic_model->getTermCriteria().epsilon;
//    statistic_model->getVotes(dataTest.row(0), vote, 0);
//    statistic_model->predict(dataTest, pred);
//    cv::Mat vote;
//    statistic_model->getVotes(dataTest, vote, 0);
//    MITK_INFO << "vote";
//    MITK_INFO << vote.rows;
//    MITK_INFO << vote.cols;
//    MITK_INFO << pred;
//    mitke_INFO << vote;
////    MITK_INFO << vote;s
//    MITK_INFO << vote;
    MITK_INFO << "_______";


    MITK_INFO << statistic_model->getPriors();
//    MITK_INFO << statistic_model->getNodes();

//    index.push_back(1);
    index_vec.push_back(index);
    index_vec.push_back(indextolabel);

    return index_vec;

}

mitk::FiberBundle::Pointer StreamlineFeatureExtractor::CreatePrediction(std::vector<unsigned int> &index)
{



    mitk::FiberBundle::Pointer Prediction;
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

     Prediction = mitk::FiberBundle::New(vNewPolyData);

//     Bundle->SetFiberColors(255, 255, 255);
    MITK_INFO << "Cells Prediciton";
    MITK_INFO << Prediction->GetFiberPolyData()->GetNumberOfCells();
    MITK_INFO << "Cells Tractorgram";
    MITK_INFO << m_TractogramTest->GetFiberPolyData()->GetNumberOfCells();

    return Prediction;
}



void  StreamlineFeatureExtractor::GenerateData()
{
    MITK_INFO << "Update";
    mitk::FiberBundle::Pointer inputPrototypes = mitk::IOUtil::Load<mitk::FiberBundle>("/home/r948e/E132-Projekte/Projects/2022_Peretzke_Interactive_Fiber_Dissection/mitk_diff/prototypes_599671.trk");

    T_Prototypes = ResampleFibers(inputPrototypes);
    T_TractogramMinus= ResampleFibers(m_TractogramMinus);
    T_TractogramPlus= ResampleFibers(m_TractogramPlus);



    m_DistancesMinus = CalculateDmdf(T_TractogramMinus, T_Prototypes);
    m_DistancesPlus = CalculateDmdf(T_TractogramPlus, T_Prototypes);


    std::ifstream f(m_DistancesTestName);


    if (f.good())
    {
        MITK_INFO << "File exists";
        m_DistancesTest.clear();
        std::ifstream myFile(m_DistancesTestName);

        if(!myFile.is_open()) throw std::runtime_error("Could not open file");
        std::string line;
        vnl_matrix<float> curline;
        curline.set_size(1, m_DistancesPlus.at(0).cols());
        curline.fill(0.0);

        float val;

        while(std::getline(myFile, line))
            {



                // Create a stringstream of the current line
                std::stringstream ss(line);
//                MITK_INFO << ss;

                // Keep track of the current column index
                int colIdx = 0;

                // Extract each integer
                while(ss >> val){

//                    // Add the current integer to the 'colIdx' column's values vector
                    curline.put(0,colIdx, val);

//                    // If the next token is a comma, ignore it and move on
//                    if(ss.peek() == ',') ss.ignore();

//                    // Increment the column index
                    colIdx++;
                }
                m_DistancesTest.push_back(curline);
            }

            // Close file
            myFile.close();
    }
    else
    {
        MITK_INFO << m_DistancesTestName;
        T_TractogramTest= ResampleFibers(m_TractogramTest);
        m_DistancesTest= CalculateDmdf(T_TractogramTest, T_Prototypes);

        std::ofstream myFile(m_DistancesTestName);
//        myFile << colname << "\n";
        for(long unsigned int i = 0; i < m_DistancesTest.size(); ++i)
        {
            myFile << m_DistancesTest.at(i);
        }
        myFile.close();
    }
    MITK_INFO << m_DistancesTest.size();


    MITK_INFO << "Sizes of Plus and Minus";
    MITK_INFO << m_DistancesPlus.size() + m_DistancesMinus.size();
    MITK_INFO << "Size of Test Data";
    MITK_INFO << m_DistancesTest.size();
    MITK_INFO << "Done with Datacreation";
    m_index =GetData();

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




