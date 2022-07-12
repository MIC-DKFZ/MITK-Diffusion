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
  : m_NumPoints(40)
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
    auto s = std::to_string(m_NumPoints);
    m_DistancesTestName= path.append("_distances" + s + ".csv");
}

std::vector<vnl_matrix<float> > StreamlineFeatureExtractor::ResampleFibers(mitk::FiberBundle::Pointer tractogram)
{
  mitk::FiberBundle::Pointer temp_fib = tractogram->GetDeepCopy();
  temp_fib->ResampleToNumPoints(m_NumPoints);
  MITK_INFO << "Resampling Done";

  std::vector< vnl_matrix<float> > out_fib(temp_fib->GetFiberPolyData()->GetNumberOfCells());
//  std::vector< vnl_matrix<float> > out_fib();

//  cv::parallel_for_(cv::Range(0, temp_fib->GetFiberPolyData()->GetNumberOfCells()), [&](const cv::Range &range)
//  {
//      for (int i = range.start; i < range.end; i++)
//    #pragma omp parallel for
//  #pragma omp parallel for num_threads(10) collapse(1)
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

//        out_fib.push_back(streamline);
        out_fib.at(i)=streamline;
      }
//      });



  return out_fib;
}

std::vector<vnl_matrix<float> > StreamlineFeatureExtractor::CalculateDmdf(std::vector<vnl_matrix<float> > tractogram, std::vector<vnl_matrix<float> > prototypes,
                                                                          std::vector<vnl_matrix<float> > local_prototypes)
{
    unsigned int locals;
    if (local_prototypes.size() >= 100)
    {
        locals = 100;
    }
    else {
        locals = local_prototypes.size();
    }
    MITK_INFO << "Locals:";
    MITK_INFO << locals;


    std::vector<vnl_matrix<float> > merged_prototypes;


    for (unsigned int k=0; k<prototypes.size(); k++ )
    {
        merged_prototypes.push_back(prototypes.at(k));
    }
    for (unsigned int k=0; k<locals; k++ )
    {
        merged_prototypes.push_back(local_prototypes.at(k));
    }

//    MITK_INFO << "merged_prototypes";
//    MITK_INFO << merged_prototypes.size();


    std::vector< vnl_matrix<float> >  dist_vec(tractogram.size());//
    MITK_INFO << "Start Calculating Dmdf";
    cv::parallel_for_(cv::Range(0, tractogram.size()), [&](const cv::Range &range)
    {
    for (int i = range.start; i < range.end; i++)

//    for (unsigned int i=0; i<tractogram.size(); i++)
    {

        vnl_matrix<float> distances;
        distances.set_size(1, local_prototypes.size());
        distances.fill(0.0);


        for (unsigned int j=0; j<local_prototypes.size(); j++)
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

                cur_dist = sqrt(pow(tractogram.at(i).get(0,ik) - local_prototypes.at(j).get(0,ik), 2.0) +
                                       pow(tractogram.at(i).get(1,ik) - local_prototypes.at(j).get(1,ik), 2.0) +
                                       pow(tractogram.at(i).get(2,ik) - local_prototypes.at(j).get(2,ik), 2.0));

                cur_dist_flip = sqrt(pow(tractogram.at(i).get(0,ik) - local_prototypes.at(j).get(0,local_prototypes.at(0).cols()-(ik+1)), 2.0) +
                                       pow(tractogram.at(i).get(1,ik) - local_prototypes.at(j).get(1,local_prototypes.at(0).cols()-(ik+1)), 2.0) +
                                       pow(tractogram.at(i).get(2,ik) - local_prototypes.at(j).get(2,local_prototypes.at(0).cols()-(ik+1)), 2.0));

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
//        dist_vec.push_back(distances);
        dist_vec.at(i) = distances;
    }
    });
    MITK_INFO << "Done Calculation";
    MITK_INFO << dist_vec.at(0).size();

    return dist_vec;
}

std::vector<std::vector<unsigned int>>  StreamlineFeatureExtractor::GetData()
{
    MITK_INFO << "Start Function Get Data";

    /*Vector which saves Prediction and Fibers to label based on uncertainty*/
    std::vector<std::vector<unsigned int>> index_vec;

    int labels_arr [m_DistancesPlus.size()+m_DistancesMinus.size()];
    cv::Mat data;
    cv::Mat labels_arr_vec;

    int size_plus = 0;


    /*Create Trainingdata: Go through positive and negative Bundle and save distances as cv::Mat and create vector with labels*/
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

    }

    /*Calculate weights*/
    int zerosgt = labels_arr_vec.rows - cv::countNonZero(labels_arr_vec);
    int onesgt = cv::countNonZero(labels_arr_vec);
    float plusval = labels_arr_vec.rows / (2.0 * onesgt );
    float minusval = labels_arr_vec.rows / (2.0 * zerosgt );

    cv::Mat newweight;
    newweight.push_back(minusval);
    newweight.push_back(plusval);
    MITK_INFO << "Weights";
    MITK_INFO << newweight;


    cv::Mat labels(m_DistancesPlus.size()+m_DistancesMinus.size(), 1, CV_32S, labels_arr);

//    std::ofstream myfile1;
//    myfile1.open("/home/r948e/mycsv/labels.csv");
//    myfile1<< cv::format(labels, cv::Formatter::FMT_CSV) << std::endl;
//    myfile1.close();

//    std::ofstream myfile2;
//    myfile2.open("/home/r948e/mycsv/features.csv");
//    myfile2<< cv::format(data, cv::Formatter::FMT_CSV) << std::endl;
//    myfile2.close();

    cv::Ptr<cv::ml::TrainData> m_traindata = cv::ml::TrainData::create(data, cv::ml::ROW_SAMPLE, labels);

//    m_traindata->setTrainTestSplitRatio(0.95, true);
    m_traindata->shuffleTrainTest();
    MITK_INFO << m_traindata->getClassLabels();

    MITK_INFO << "Start Training";

    auto statistic_model = cv::ml::RTrees::create();

    auto criteria = cv::TermCriteria();
    criteria.type = cv::TermCriteria::MAX_ITER;
//    criteria.epsilon = 1e-8;
    criteria.maxCount = 800;

//    statistic_model->setMaxCategories(2);
    statistic_model->setMaxDepth(50); //set to three
//    statistic_model->setMinSampleCount(m_traindata->getNTrainSamples()*0.01);
    statistic_model->setMinSampleCount(2);
    statistic_model->setTruncatePrunedTree(false);
    statistic_model->setUse1SERule(false);
    statistic_model->setUseSurrogates(false);
    statistic_model->setTermCriteria(criteria);
    statistic_model->setCVFolds(1);
    statistic_model->setPriors(newweight);


    statistic_model->train(m_traindata);




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


    std::vector<unsigned int> index;
    std::vector<float> e(m_DistancesTest.size());

//    cv::setNumThreads(16);
//    cv::parallel_for_(cv::Range(0, m_DistancesTest.size()), [&](const cv::Range &range)
//    {
//    for (int i = range.start; i < range.end; i++)
#pragma omp parallel for
    for (unsigned int i=0; i<m_DistancesTest.size(); i++)
    {


        int val = statistic_model->predict(dataTest.row(i));
        if (val==1)
        {
           index.push_back(i);
        }


        cv::Mat vote;
        statistic_model->getVotes(dataTest.row(i), vote, 0);
        e.at(i) = ( -(vote.at<int>(1,0)*1.0)/ (vote.at<int>(1,0)+vote.at<int>(1,1)) * log2((vote.at<int>(1,0)*1.0)/ (vote.at<int>(1,0)+vote.at<int>(1,1))) -
                    (vote.at<int>(1,1)*1.0)/ (vote.at<int>(1,0)+vote.at<int>(1,1))* log2((vote.at<int>(1,1)*1.0)/ (vote.at<int>(1,0)+vote.at<int>(1,1))));

        if (isnan(e.at(i)))
        {
            e.at(i)=0;
//            MITK_INFO << e.at(i);
        }
        if (i==1)
        {
            MITK_INFO<< val;
            MITK_INFO << vote;

        }
    }
//    });
//    std::ofstream myfile3;
//    myfile3.open("/home/r948e/mycsv/entropydata.csv");

//  for (unsigned int i = 0; i < e.size(); i++) {
//          myfile3 << e.at(i) << ' ';
//      }
//  myfile3.close();




    MITK_INFO << "--------------";
    MITK_INFO << "Prediction vector size:";
    MITK_INFO << index.size();
    MITK_INFO << "Entropy vector size:";
    MITK_INFO << e.size();

    MITK_INFO << "--------------";

    auto it = std::minmax_element(e.begin(), e.end());
    int min_idx = std::distance(e.begin(), it.first);
    int max_idx = std::distance(e.begin(), it.second);
    std::cout << min_idx << ", " << max_idx << std::endl; // 1, 5
    MITK_INFO << e.at(max_idx);


    MITK_INFO << "Start the ordering";
    std::vector<unsigned int> indextolabel;
    std::priority_queue<std::pair<float, int>> q;
    for (unsigned int i = 0; i < e.size(); ++i)
    {
        q.push(std::pair<float, int>(e[i], i));
    }
//    int k = m_DistancesTest.size(); // number of indices we need
    int lengths=500;
//    int k = lengths; // number of indices we need
    for (int i = 0; i < lengths; ++i)
    {
        int ki = q.top().second;
        indextolabel.push_back(ki);
        q.pop();
    }

//    std::ofstream myfile4;
//    myfile4.open("/home/r948e/mycsv/indextolabel.csv");
//    for (unsigned int i = 0; i < indextolabel.size(); i++)
//    {
//        myfile4 << indextolabel.at(i) << ' ';
//    }
//    myfile4.close();
    // Sorted depent on entropy now sort the 1/5
    vnl_matrix<float> distances_matrix;

    distances_matrix.set_size(lengths, lengths);
    distances_matrix.fill(0.0);
    MITK_INFO << "Dist_start";

    std::vector<float> distances_matrix_mean;


    for (int i=0; i<lengths; i++)
    {
        for (int k=0; k<lengths; k++)
        {
            vnl_matrix<float> diff =  m_DistancesTest.at(indextolabel.at(i)) - m_DistancesTest.at(indextolabel.at(k));

            distances_matrix.put(i,k,diff.absolute_value_sum()/m_DistancesTest.at(0).size());

        }
        distances_matrix_mean.push_back(distances_matrix.get_row(i).mean());
//        MITK_INFO << meanval.at(i);

    }

    vnl_vector<float> sum_matrix;
    sum_matrix.set_size(lengths);
    sum_matrix.fill(0.0);

    /*Index to find values is distancematrix*/
    std::vector<unsigned int> myidx;
    /*Index to find actual streamlines using indextolabel*/
    std::vector<unsigned int> distindextolabel;
    myidx.push_back(0);

//    MITK_INFO << distances_matrix.get_row(myidx.at(i)+ sum_matrix.get_row(0)


    for (int i=0; i<lengths; i++)
    {
//        unsigned int cur_i = indextolabel.at(myidx.at(i));
        sum_matrix = (sum_matrix + distances_matrix.get_row(myidx.at(i)))/=(i+1);

//        myidx.push_back(distances_matrix.get_row(myidx.at(i)).arg_max());
        myidx.push_back(sum_matrix.arg_max());
        distindextolabel.push_back(indextolabel.at(myidx.at(i)));
        distances_matrix.set_column(myidx.at(i), 0.0001);

    }

    std::vector<unsigned int> indextolabeldist;
    std::priority_queue<std::pair<float, int>> qq;
    for (unsigned int i = 0; i < distances_matrix_mean.size(); ++i)
    {
        qq.push(std::pair<float, int>(distances_matrix_mean[i], i));
    }
//    int k = m_DistancesTest.size(); // number of indices we need
//    int k = lengths; // number of indices we need
    for (int i = 0; i < lengths; ++i)
    {
        int kki = qq.top().second;

        indextolabeldist.push_back(indextolabel.at(kki));
        qq.pop();
    }

    MITK_INFO << "Dist_stop";

//    std::ofstream myfile6;
//    myfile6.open("/home/r948e/mycsv/distances_matrix_mean.csv");
//    for (unsigned int i = 0; i < distances_matrix_mean.size(); i++)
//    {
//        myfile6 << distances_matrix_mean.at(i) << ' ';
//    }
//    myfile6.close();

//    std::ofstream myfile5;
//    myfile5.open("/home/r948e/mycsv/indextolabeldist.csv");
//    for (unsigned int i = 0; i < indextolabeldist.size(); i++)
//    {
//        myfile5 << indextolabeldist.at(i) << ' ';
//    }
//    myfile5.close();

//    MITK_INFO << distances_matrix;
//    MITK_INFO << distances_matrix.max_value();
//    MITK_INFO << distances_matrix.arg_max();


//    vnl_matrix<float> myx =  m_DistancesTest.at(indextolabel.at(0)) - m_DistancesTest.at(indextolabel.at(1));
////    myx = (m_DistancesTest.at(indextolabel.at(0)) - m_DistancesTest.at(indextolabel.at(1)));
//    MITK_INFO << m_DistancesTest.at(indextolabel.at(0));
//    MITK_INFO << m_DistancesTest.at(indextolabel.at(1));
//    MITK_INFO << myx.get(0,0);
//    MITK_INFO << sqrt(pow(myx.get(0,0),2));
//    MITK_INFO << myx.get(0,1);
//    MITK_INFO << myx.absolute_value_sum()/m_DistancesTest.at(0).size();
//    MITK_INFO << "Done";


    index_vec.push_back(index);
    index_vec.push_back(indextolabel);
    index_vec.push_back(distindextolabel);

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


    MITK_INFO << "Calculate Features";
    m_DistancesMinus = CalculateDmdf(T_TractogramMinus, T_Prototypes, T_TractogramPlus);
    m_DistancesPlus = CalculateDmdf(T_TractogramPlus, T_Prototypes, T_TractogramPlus);




    std::ifstream f(m_DistancesTestName);

    MITK_INFO << m_DistancesTestName;
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
        MITK_INFO << "Resample Test Data";
        T_TractogramTest= ResampleFibers(m_TractogramTest);
        MITK_INFO << "Calculate Features of Test Data";
        m_DistancesTest= CalculateDmdf(T_TractogramTest, T_Prototypes, T_TractogramPlus);

//        std::ofstream myFile(m_DistancesTestName);
////        myFile << colname << "\n";
//        for(long unsigned int i = 0; i < m_DistancesTest.size(); ++i)
//        {
//            myFile << m_DistancesTest.at(i);
//        }
//        myFile.close();
    }
    MITK_INFO << m_DistancesTest.size();


    MITK_INFO << "Sizes of Plus and Minus";
    MITK_INFO << m_DistancesPlus.size() + m_DistancesMinus.size();
    MITK_INFO << "Size of Test Data";
    MITK_INFO << m_DistancesTest.size();
    MITK_INFO << "Done with Datacreation";
    m_index =GetData();

}



}




