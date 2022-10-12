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

void StreamlineFeatureExtractor::SetTractogramPrediction(const mitk::FiberBundle::Pointer &TractogramPrediction)
{
  m_TractogramPrediction= TractogramPrediction;
}

void StreamlineFeatureExtractor::SetTractogramGroundtruth(const mitk::FiberBundle::Pointer &TractogramGroundtruth)
{
  m_TractogramGroundtruth= TractogramGroundtruth;
}

void StreamlineFeatureExtractor::SetTractogramPlus(const mitk::FiberBundle::Pointer &TractogramPlus)
{
  m_TractogramPlus = TractogramPlus;
}

void StreamlineFeatureExtractor::SetTractogramMinus(const mitk::FiberBundle::Pointer &TractogramMinus)
{
  m_TractogramMinus = TractogramMinus;
}

void StreamlineFeatureExtractor::SetTractogramPrototypes(const mitk::FiberBundle::Pointer &TractogramPrototypes, bool standard)
{
  if (standard)
  {
      MITK_INFO << "Use Standard Prototypes...";
    m_inputPrototypes = mitk::IOUtil::Load<mitk::FiberBundle>("/home/r948e/E132-Projekte/Projects/2022_Peretzke_Interactive_Fiber_Dissection/mitk_diff/prototypes_599671_40.trk");
  }
  else {
      MITK_INFO << "Use individual Prototypes...";
      m_inputPrototypes = TractogramPrototypes;
  }


}

void StreamlineFeatureExtractor::SetActiveCycle(int &activeCycle)
{
  m_activeCycle= activeCycle;
}

void StreamlineFeatureExtractor::SetInitRandom(int &initRandom)
{
  m_initRandom= initRandom;
}

void StreamlineFeatureExtractor::SetTractogramTest(const mitk::FiberBundle::Pointer &TractogramTest, std::string TractogramTestName)
{
//    std::string path = "/home/r948e/E132-Projekte/Projects/2022_Peretzke_Interactive_Fiber_Dissection/mitk_diff/storage/";
//    path.append(TractogramTestName);
    MITK_INFO << TractogramTestName;
    m_TractogramTest= TractogramTest;
//    auto s = std::to_string(m_NumPoints);
//    m_DistancesTestName= path.append("_distances" + std::to_string(m_NumPoints) + "_" + std::to_string(m_activeCycle) + ".csv");
}

std::vector<vnl_matrix<float> > StreamlineFeatureExtractor::ResampleFibers(mitk::FiberBundle::Pointer tractogram)
{
  MITK_INFO << "Infunction";
  mitk::FiberBundle::Pointer temp_fib = tractogram->GetDeepCopy();
//  temp_fib->ResampleToNumPoints(m_NumPoints);
  MITK_INFO << "Resampling Done";

  std::vector< vnl_matrix<float> > out_fib(temp_fib->GetFiberPolyData()->GetNumberOfCells());

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

std::vector<vnl_matrix<float> > StreamlineFeatureExtractor::CalculateDmdf(std::vector<vnl_matrix<float> > tractogram,
                                                                          std::vector<vnl_matrix<float> > prototypes)
{

    std::vector< vnl_matrix<float> >  dist_vec(tractogram.size());//
    MITK_INFO << "Start Calculating Dmdf";
    cv::parallel_for_(cv::Range(0, tractogram.size()), [&](const cv::Range &range)
    {
    for (int i = range.start; i < range.end; i++)

//    for (unsigned int i=0; i<tractogram.size(); i++)
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
//        dist_vec.push_back(distances);
        dist_vec.at(i) = distances;
    }
    });
    MITK_INFO << "Done Calculation";
    MITK_INFO << dist_vec.at(0).size();

    return dist_vec;
}

std::vector<vnl_matrix<float> > StreamlineFeatureExtractor::MergeTractogram(std::vector<vnl_matrix<float> > prototypes,
                                                                           std::vector<vnl_matrix<float> > positive_local_prototypes,
                                                                           std::vector<vnl_matrix<float> > negative_local_prototypes)
{
    unsigned int pos_locals;
    unsigned int neg_locals;

//    if (positive_local_prototypes.size() >= 50)
//    {
//        pos_locals= 50;
//    }
//    else {
//        pos_locals= positive_local_prototypes.size();
//    }

//    if (pos_locals <= negative_local_prototypes.size())
//    {
//        neg_locals = pos_locals;
//    }
//    else {
//        neg_locals= negative_local_prototypes.size();
//    }



    std::vector< vnl_matrix<float> > merged_prototypes;

    for (unsigned int k=0; k<prototypes.size(); k++)
    {
        merged_prototypes.push_back(prototypes.at(k));
    }

    for (unsigned int k=0; k<positive_local_prototypes.size(); k++)
    {
        merged_prototypes.push_back(positive_local_prototypes.at(k));
    }

    for (unsigned int k=0; k<negative_local_prototypes.size(); k++)
    {
        merged_prototypes.push_back(negative_local_prototypes.at(k));
    }

    MITK_INFO << "Number of prototypes:";
    MITK_INFO << prototypes.size();
    MITK_INFO << "Number of positive local prototypes:";
    MITK_INFO << pos_locals;
    MITK_INFO << "Number of negative local prototypes:";
    MITK_INFO << neg_locals;

    return merged_prototypes;

}


std::vector<unsigned int> StreamlineFeatureExtractor::Sort(std::vector<float> sortingVector, int lengths, int start)
{
    std::vector<unsigned int> index;
    std::priority_queue<std::pair<float, int>> q;

    for (unsigned int i = 0; i < sortingVector.size(); ++i)
    {
        q.push(std::pair<float, int>(sortingVector[i], i));
    }


    for (int i = 0; i < lengths; ++i)
    {
        int ki = q.top().second;
        if (i>=start)
        {
        index.push_back(ki);
        }
        q.pop();
    }
    return index;
}


std::vector<std::vector<unsigned int>>  StreamlineFeatureExtractor::GetData()
{
    MITK_INFO << "Start Function Get Data";

    /*Vector which saves Prediction and Fibers to label based on uncertainty*/
    std::vector<std::vector<unsigned int>> index_vec;

//    int labels_arr [m_DistancesPlus.size()+m_DistancesMinus.size()];
    cv::Mat data;
    cv::Mat labels_arr_vec;

    int size_plus = 0;


    /*Create Trainingdata: Go through positive and negative Bundle and save distances as cv::Mat and create vector with labels*/
    for ( unsigned int i=0; i<m_DistancesPlus.size(); i++)
    {
        float data_arr [m_DistancesPlus.at(0).size()];
//        labels_arr[i]=1;
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
//        labels_arr[i]=0;
        labels_arr_vec.push_back(0);

         for ( unsigned int j=0; j<m_DistancesMinus.at(0).cols(); j++)
        {
            data_arr[j] = m_DistancesMinus.at(it).get(0,j);
        }
        cv::Mat curdata(1, m_DistancesPlus.at(0).size(), CV_32F, data_arr);
        data.push_back(curdata);

    }
    MITK_INFO << data.cols;
    MITK_INFO << data.rows;

    /*Calculate weights*/
    int zerosgt = labels_arr_vec.rows - cv::countNonZero(labels_arr_vec);
    int onesgt = cv::countNonZero(labels_arr_vec);
    float plusval = labels_arr_vec.rows / (2.0 * onesgt ) ;
    float minusval = labels_arr_vec.rows / (2.0 * zerosgt );
    float w[2] = {minusval, plusval};
    cv::Mat newweight = cv::Mat(1,2, CV_32F, w);
    MITK_INFO << "Weights";
    MITK_INFO << newweight;


    /*Shuffle Data*/
    std::vector <int> seeds;
    for (int cont = 0; cont < labels_arr_vec.rows; cont++)
    {
        seeds.push_back(cont);
    }

    cv::randShuffle(seeds);
    cv::Mat labels_shuffled;
    cv::Mat samples_shuffled;




    for (int cont = 0; cont < labels_arr_vec.rows; cont++)
    {
        labels_shuffled.push_back(labels_arr_vec.row(seeds[cont]));
    }

    for (int cont = 0; cont < labels_arr_vec.rows; cont++)
    {
        samples_shuffled.push_back(data.row(seeds[cont]));
    }


//    std::ofstream labelsfile;
//    labelsfile.open("/home/r948e/E132-Projekte/Projects/2022_Peretzke_Interactive_Fiber_Dissection/mitk_diff/storage/Labels_" + std::to_string(m_activeCycle) + ".csv");
//    labelsfile<< cv::format(labels_shuffled, cv::Formatter::FMT_CSV) << std::endl;
//    labelsfile.close();

//    std::ofstream featuresfile;
//    featuresfile.open("/home/r948e/E132-Projekte/Projects/2022_Peretzke_Interactive_Fiber_Dissection/mitk_diff/storage/Features_" + std::to_string(m_activeCycle) + ".csv");
//    featuresfile<< cv::format(samples_shuffled, cv::Formatter::FMT_CSV) << std::endl;
//    featuresfile.close();

    /*Create Dataset and initialize Classifier*/
    cv::Ptr<cv::ml::TrainData> m_traindata = cv::ml::TrainData::create(samples_shuffled, cv::ml::ROW_SAMPLE, labels_shuffled);




    statistic_model= cv::ml::RTrees::create();
    auto criteria = cv::TermCriteria();
    criteria.type = cv::TermCriteria::MAX_ITER;
//    criteria.epsilon = 1e-8;
    criteria.maxCount = 400;

    statistic_model->setMaxDepth(10); //set to three
//    statistic_model->setMinSampleCount(m_traindata->getNTrainSamples()*0.01);
    statistic_model->setMinSampleCount(2);
    statistic_model->setTruncatePrunedTree(false);
    statistic_model->setUse1SERule(false);
    statistic_model->setUseSurrogates(false);
    statistic_model->setTermCriteria(criteria);
    statistic_model->setCVFolds(1);
    statistic_model->setPriors(newweight);


    /*Train Classifier*/
    MITK_INFO << "Start Training";
    statistic_model->train(m_traindata);


    /*Predict on Test Data*/
    MITK_INFO << "Predicting";

    /*Create Dataset as cv::Mat*/
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


    std::vector<unsigned int> indexPrediction;
    std::vector<float> e(m_DistancesTest.size());
    std::vector<int> pred(m_DistancesTest.size());


    /*For every Sample/Streamline get Prediction and entropy (=based on counts of Random Forest)*/
    MITK_INFO << "Predicting on all cores";
#pragma omp parallel for
    for (unsigned int i=0; i<m_DistancesTest.size(); i++)
    {


        int val = statistic_model->predict(dataTest.row(i));
        pred.at(i)=val;

        #pragma omp critical
        if (val==1)
        {
           indexPrediction.push_back(i);
        }


        cv::Mat vote;
        statistic_model->getVotes(dataTest.row(i), vote, 0);
        e.at(i) = ( -(vote.at<int>(1,0)*1.0)/ (vote.at<int>(1,0)+vote.at<int>(1,1)) * log2((vote.at<int>(1,0)*1.0)/ (vote.at<int>(1,0)+vote.at<int>(1,1))) -
                    (vote.at<int>(1,1)*1.0)/ (vote.at<int>(1,0)+vote.at<int>(1,1))* log2((vote.at<int>(1,1)*1.0)/ (vote.at<int>(1,0)+vote.at<int>(1,1))));

        if (isnan(e.at(i)))
        {
            e.at(i)=0;
        }

    }
    MITK_INFO << "Done";

    /*Save entropy values for analysis*/
    std::ofstream entropyfile;
    entropyfile.open("/home/r948e/mycsv/entropydata_" + std::to_string(m_activeCycle) + ".csv");
    for (unsigned int i = 0; i < e.size(); i++)
    {
        entropyfile << e.at(i) << ' ';
    }
    entropyfile.close();

    std::ofstream predictionfile;
    predictionfile.open("/home/r948e/mycsv/predictiondata_" + std::to_string(m_activeCycle) + ".csv");
    for (unsigned int i = 0; i < pred.size(); i++)
    {
        predictionfile << pred.at(i) << ' ';
    }
    predictionfile.close();




    MITK_INFO << "--------------";
    MITK_INFO << "Prediction vector size:";
    MITK_INFO << indexPrediction.size();
    MITK_INFO << "Entropy vector size:";
    entropy_vector = e;
    MITK_INFO << e.size();

    MITK_INFO << "--------------";

    /*Get index of most unertain data (lengths defines how many data is saved)*/
//    int lengths=500;
  int lengths = std::count_if(e.begin(), e.end(),[&](auto const& val){ return val >= 0.95; });
  if (lengths>1500)
  {
      lengths=1500;
  }
//  else if (lengths>3000)
//  {
//      lengths=3000;
//  }

    int lengthsCertain = std::count_if(e.begin(), e.end(),[&](auto const& val){ return val < 0.1; });

    std::vector<unsigned int> indexUnc = Sort(e, lengths, 0);

    std::vector<unsigned int> indexCertain = Sort(e, e.size() , e.size()-lengthsCertain );

//    std::vector<unsigned int> indexCertainBetween = Sort(e, e.size()-lengthsCertain , lengths);

    MITK_INFO << "Index Certainty Vector size";
    MITK_INFO << indexCertain.size();

//    for (unsigned int i=indexCertain.size(); i>=0; --i)
    std::vector<unsigned int> indexCertainNeg;
    std::vector<unsigned int> indexCertainPos;

    for (unsigned int i=0; i<indexCertain.size(); i++)
    {
        if(pred.at(indexCertain.at(i))==0 )
        {
            indexCertainNeg.push_back(indexCertain.at(i));
        }
        else {
            indexCertainPos.push_back(indexCertain.at(i));
        }

    }


    //    for (unsigned int i=indexCertain.size(); i>=0; --i)
//    std::vector<unsigned int> indexCertainBetweenNeg;
//    std::vector<unsigned int> indexCertainBetweenPos;

//    for (unsigned int i=0; i<indexCertainBetween.size(); i++)
//    {
//        if(pred.at(indexCertainBetween.at(i))==0 )
//        {
//            indexCertainBetweenNeg.push_back(indexCertainBetween.at(i));
//        }
//        else {
//            indexCertainBetweenPos.push_back(indexCertainBetween.at(i));
//        }

//    }

    MITK_INFO << "Index Uncertainty Vector size";
    MITK_INFO << indexUnc.size();
    MITK_INFO << "Index Certainty Vector size";
    MITK_INFO << indexCertainPos.size();
    MITK_INFO << indexCertainNeg.size();
    MITK_INFO << indexCertainNeg.size() +indexCertainPos.size();
    MITK_INFO << "Index Certainty between Vector size";
//    MITK_INFO << indexCertainBetweenPos.size();
//    MITK_INFO << indexCertainBetweenNeg.size();
//    MITK_INFO << indexCertainBetweenNeg.size() +indexCertainBetweenPos.size();


//    vnl_matrix<float> distances_matrix;

//    distances_matrix.set_size(lengths, lengths);
//    distances_matrix.fill(0.0);

//    std::vector<float> distances_matrix_mean;


//    for (int i=0; i<lengths; i++)
//    {
//        for (int k=0; k<lengths; k++)
//        {
//            /*From the length.size() Samples with the highest Entropey calculate the differen between the Features*/
//            vnl_matrix<float> diff =  m_DistancesTest.at(indexUnc.at(i)) - m_DistancesTest.at(indexUnc.at(k));

//            /*Into the eucledean difference matrix, put the distance in Feature Space between every sample pair*/
//            distances_matrix.put(i,k,diff.absolute_value_sum()/m_DistancesTest.at(0).size());

//        }
//        /*For every Sample/Streamline get the mean eucledean distance to all other Samples => one value for every Sample*/
////        distances_matrix_mean.push_back(distances_matrix.get_row(i).mean());
////        MITK_INFO << meanval.at(i);

//    }



//    /*Index to find values in distancematrix*/
//    std::vector<unsigned int> myidx;
//    /*Index to find actual streamlines using indexUnc*/
//    std::vector<unsigned int> indexUncDist;
//    /*Start with the Streamline of the highest entropy, which is in distance_matrix at idx 0*/
//    myidx.push_back(0);
//    indexUncDist.push_back(indexUnc.at(myidx.at(0)));

//    /*Vecotr that stores minvalues of current iteration*/
//    vnl_matrix<float> cur_vec;
//    cur_vec.set_size(1,lengths);
//    cur_vec.fill(0.0);
//    for (int i=0; i<lengths; i++)
//    {

////        unsigned int cur_i = indexUnc.at(myidx.at(i));

//        /*Save mean distance of all used Samples*/
//        vnl_matrix<float> sum_matrix;
//        sum_matrix.set_size(myidx.size(), lengths);
//        sum_matrix.fill(0);
//        for (unsigned int ii=0; ii<myidx.size(); ii++)
//        {

//            sum_matrix.set_row(ii, distances_matrix.get_column(myidx.at(ii)));
//        }

//        for (unsigned int k=0; k<sum_matrix.columns(); k++)
//        {
//            cur_vec.put(0,k, sum_matrix.get_column(k).min_value());
//        }
//        myidx.push_back(cur_vec.arg_max());


//        indexUncDist.push_back(indexUnc.at(myidx.at(i+1)));
//        sum_matrix.clear();

//    }

//    MITK_INFO << "Dist_stop";


    /*Save Prediction*/
    index_vec.push_back(indexPrediction);
    /*Save index of uncertainty measures*/
    index_vec.push_back(indexUnc);
    /*Save index of uncertainty measures influenced by distance*/
//    index_vec.push_back(indexUncDist);
    /*Save index of certain measures*/
    index_vec.push_back(indexCertainNeg);
    /*Save index of certain measures*/
//    index_vec.push_back(indexCertainPos);
    /*Save index of certain measures*/
//    index_vec.push_back(indexCertainBetweenNeg);
    /*Save index of certain measures*/
//    index_vec.push_back(indexCertainBetweenPos);

    return index_vec;

}

std::vector<std::vector<unsigned int>>  StreamlineFeatureExtractor::GetDistanceData(float &value)
{
    /*Vector which saves Fibers to be labeled based on fft subset uncertainty*/
    std::vector<std::vector<unsigned int>> index_vec;

    /*Get index of most unertain data (lengths defines how many data is saved)*/
//    int lengths=500;
    MITK_INFO << entropy_vector.size();
  int lengths = std::count_if(entropy_vector.begin(), entropy_vector.end(),[&](auto const& val){ return val >= value; });
  if (lengths>1500)
  {
      lengths=1500;
  }
  MITK_INFO << lengths;
    /*Maybe shuffling of length so not the most uncertain values are chosen*/
    std::vector<unsigned int> indexUnc = Sort(entropy_vector, lengths, 0);

    vnl_matrix<float> distances_matrix;

    distances_matrix.set_size(lengths, lengths);
    distances_matrix.fill(0.0);

    std::vector<float> distances_matrix_mean;


    for (int i=0; i<lengths; i++)
    {
        for (int k=0; k<lengths; k++)
        {
            /*From the length.size() Samples with the highest Entropey calculate the differen between the Features*/
            vnl_matrix<float> diff =  m_DistancesTest.at(indexUnc.at(i)) - m_DistancesTest.at(indexUnc.at(k));

            /*Into the eucledean difference matrix, put the distance in Feature Space between every sample pair*/
            distances_matrix.put(i,k,diff.absolute_value_sum()/m_DistancesTest.at(0).size());

        }
        /*For every Sample/Streamline get the mean eucledean distance to all other Samples => one value for every Sample*/
//        distances_matrix_mean.push_back(distances_matrix.get_row(i).mean());
//        MITK_INFO << meanval.at(i);

    }

    MITK_INFO << "Distance Matrix Calculated";

    /*Index to find values in distancematrix*/
    std::vector<unsigned int> myidx;
    /*Index to find actual streamlines using indexUnc*/
    std::vector<unsigned int> indexUncDist;
    /*Start with the Streamline of the highest entropy, which is in distance_matrix at idx 0*/
    myidx.push_back(0);
    indexUncDist.push_back(indexUnc.at(myidx.at(0)));

    /*Vecotr that stores minvalues of current iteration*/
    vnl_matrix<float> cur_vec;
    cur_vec.set_size(1,lengths);
    cur_vec.fill(0.0);
    for (int i=0; i<lengths; i++)
    {

//        unsigned int cur_i = indexUnc.at(myidx.at(i));

        /*Save mean distance of all used Samples*/
        vnl_matrix<float> sum_matrix;
        sum_matrix.set_size(myidx.size(), lengths);
        sum_matrix.fill(0);
        for (unsigned int ii=0; ii<myidx.size(); ii++)
        {

            sum_matrix.set_row(ii, distances_matrix.get_column(myidx.at(ii)));
        }

        for (unsigned int k=0; k<sum_matrix.columns(); k++)
        {
            cur_vec.put(0,k, sum_matrix.get_column(k).min_value());
        }
        myidx.push_back(cur_vec.arg_max());


        indexUncDist.push_back(indexUnc.at(myidx.at(i+1)));
        sum_matrix.clear();

    }

    MITK_INFO << "Dist_stop";

    index_vec.push_back(indexUncDist);

    return index_vec;
}

mitk::FiberBundle::Pointer StreamlineFeatureExtractor::CreatePrediction(std::vector<unsigned int> &index, bool removefrompool)
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

//    vtkSmartPointer<vtkFloatArray> weights = vtkSmartPointer<vtkFloatArray>::New();
//    weights->SetNumberOfValues(this->GetNumFibers()+fib->GetNumFibers());

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
//            weights->InsertValue(counter, fib->GetFiberWeight(i));
        vNewLines->InsertNextCell(container);
        counter++;

    }
    if (removefrompool)
    {
        for (unsigned int i=0; i<indexSize; i++)
        {
            m_TractogramTest->GetFiberPolyData()->DeleteCell(index[i]);
        }
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

std::vector<int> StreamlineFeatureExtractor::CreateLabels(std::vector<vnl_matrix<float> > Testdata,
                                                              std::vector<vnl_matrix<float> > Prediction)
{
//    vnl_vector<int> labels;
//    vnl_vector.set_size(Testdata.size());
//    vnl_vector.fill(0);
    std::vector<int>  labels(Testdata.size(), 0);


#pragma omp parallel for
    for (unsigned int i=0; i<Prediction.size(); i++)
    {
        for (unsigned int k=0; k<Testdata.size(); k++)
        {
            if (Testdata.at(k)==Prediction.at(i))
            {
                labels.at(k)=1;
            }
        }
    }


    return labels;




}

void StreamlineFeatureExtractor::GenerateData()
{
    MITK_INFO << "Update";

//    mitk::FiberBundle::Pointer inputPrototypes = mitk::IOUtil::Load<mitk::FiberBundle>("/home/r948e/E132-Projekte/Projects/2022_Peretzke_Interactive_Fiber_Dissection/data/Synt_tract_40_prototypes.trk");

    std::vector<vnl_matrix<float> >             T_Prototypes;
    std::vector<vnl_matrix<float> >             T_TractogramPlus;
    std::vector<vnl_matrix<float> >             T_TractogramMinus;
    std::vector<vnl_matrix<float> >             T_TractogramTest;
    std::vector<vnl_matrix<float> >             T_mergedPrototypes;


    MITK_INFO << "Resample Input Prototypes";
    T_Prototypes = ResampleFibers(m_inputPrototypes);
    MITK_INFO << "Resample Input Tractogram Minus";
    T_TractogramMinus= ResampleFibers(m_TractogramMinus);
    MITK_INFO << "Resample Input Tractogram Plus";
    T_TractogramPlus= ResampleFibers(m_TractogramPlus);

    /* Merge T_Prototypes, T_TractogramMinus and T_TractogramPlus for extra Features*/
    MITK_INFO << "Merging Prototypes";
    T_mergedPrototypes = MergeTractogram(T_Prototypes, T_TractogramPlus, T_TractogramMinus);




    MITK_INFO << "Calculate Features";
    MITK_INFO << "Calculate Minus Features";
    m_DistancesMinus = CalculateDmdf(T_TractogramMinus, T_mergedPrototypes);
//    std::ofstream myFileminus("/home/r948e/E132-Projekte/Projects/2022_Peretzke_Interactive_Fiber_Dissection/mitk_diff/storage/dist_minus.csv");
//    for(long unsigned int i = 0; i < m_DistancesMinus.size(); ++i)
//    {
//        myFileminus << m_DistancesMinus.at(i);
//    }
//    myFileminus.close();

    MITK_INFO << "Calculate Plus Features";
    m_DistancesPlus = CalculateDmdf(T_TractogramPlus, T_mergedPrototypes);

//    std::ofstream myFileplus("/home/r948e/E132-Projekte/Projects/2022_Peretzke_Interactive_Fiber_Dissection/mitk_diff/storage/dist_plus.csv");
//    for(long unsigned int i = 0; i < m_DistancesPlus.size(); ++i)
//    {
//        myFileplus << m_DistancesPlus.at(i);
//    }
//    myFileplus.close();



//    MITK_INFO << m_DistancesTestName;
    MITK_INFO << "Resample Test Data";


    T_TractogramTest= ResampleFibers(m_TractogramTest);

//    std::ofstream Prototypesfibers;
//    Prototypesfibers.open("/home/r948e/mycsv/Prototypefibers_" + std::to_string(m_activeCycle) + ".csv");
//    for (unsigned int i = 0; i < T_mergedPrototypes.size(); i++)
//    {
//        Prototypesfibers<< T_mergedPrototypes.at(i) << ' ';
//    }
//    Prototypesfibers.close();

//    std::ofstream TestfbiersFile;
//    TestfbiersFile.open("/home/r948e/mycsv/Testfibers_" + std::to_string(m_activeCycle) + ".csv");
//    for (unsigned int i = 0; i < T_TractogramTest.size(); i++)
//    {
//        TestfbiersFile << T_TractogramTest.at(i) << ' ';
//    }
//    TestfbiersFile.close();



//    MITK_INFO << "Calculate Features of Test Data";
    m_DistancesTest= CalculateDmdf(T_TractogramTest, T_mergedPrototypes);

//            std::ofstream myFile(m_DistancesTestName);
//            for(long unsigned int i = 0; i < m_DistancesTest.size(); ++i)
//            {
//                myFile << m_DistancesTest.at(i);
//            }
//            myFile.close();

//    std::ifstream f(m_DistancesTestName);


//    if (f.good() && m_activeCycle!=0)
//    {
//        MITK_INFO << "Loading Features of Tractogram";
//        m_DistancesTest.clear();
//        std::ifstream myFile(m_DistancesTestName);

//        if(!myFile.is_open()) throw std::runtime_error("Could not open file");
//        std::string line;
//        vnl_matrix<float> curline;
//        curline.set_size(1, m_DistancesPlus.at(0).cols());
//        curline.fill(0.0);

//        float val;

//        while(std::getline(myFile, line))
//            {
//                // Create a stringstream of the current line
//                std::stringstream ss(line);
////                MITK_INFO << ss;

//                // Keep track of the current column index
//                int colIdx = 0;

//                // Extract each integer
//                while(ss >> val){

////                    // Add the current integer to the 'colIdx' column's values vector
//                    curline.put(0,colIdx, val);

////                    // If the next token is a comma, ignore it and move on
////                    if(ss.peek() == ',') ss.ignore();
////                    // Increment the column index
//                    colIdx++;
//                }
//                m_DistancesTest.push_back(curline);
//            }

//            // Close file
//            myFile.close();
//    }
//    else
//    {
//        MITK_INFO << m_DistancesTestName;
//        MITK_INFO << "Resample Test Data";
//        T_TractogramTest= ResampleFibers(m_TractogramTest);
//        MITK_INFO << "Calculate Features of Test Data";
//        m_DistancesTest= CalculateDmdf(T_TractogramTest, T_mergedPrototypes);

//                std::ofstream myFile(m_DistancesTestName);
//        //        myFile << colname << "\n";
//                for(long unsigned int i = 0; i < m_DistancesTest.size(); ++i)
//                {
//                    myFile << m_DistancesTest.at(i);
//                }
//                myFile.close();

//    }
//    MITK_INFO << m_DistancesTest.size();


    MITK_INFO << "Sizes of Plus and Minus";
    MITK_INFO << m_DistancesPlus.size() + m_DistancesMinus.size();
    MITK_INFO << "Sizes of Prototypes";
    MITK_INFO << T_mergedPrototypes.size();
    MITK_INFO << T_mergedPrototypes.at(0).rows();
    MITK_INFO << "Size of Test Data";
    MITK_INFO << m_DistancesTest.size();
    MITK_INFO << "Done with Datacreation";
    m_index =GetData();

}

vnl_vector<float> StreamlineFeatureExtractor::ValidationPipe()
{
    std::vector<vnl_matrix<float> >             T_Prototypes;
    std::vector<vnl_matrix<float> >             T_TractogramPrediction;
    std::vector<vnl_matrix<float> >             T_TractogramGroundtruth;
    std::vector<vnl_matrix<float> >             T_TractogramTest;
    std::vector<vnl_matrix<float> >             DistancesPrediction;
    std::vector<vnl_matrix<float> >             DistancesGroundtruth;
    std::vector<vnl_matrix<float> >             DistancesTest;
    std::vector<int>                            LabelsPrediction;
    std::vector<int>                            LabelsGroundtruth;

    MITK_INFO << "Start Resampling";
    T_Prototypes = ResampleFibers(m_inputPrototypes);
    T_TractogramPrediction= ResampleFibers(m_TractogramPrediction);
    T_TractogramGroundtruth= ResampleFibers(m_TractogramGroundtruth);
    T_TractogramTest= ResampleFibers(m_TractogramTest);



    MITK_INFO << "Calculate Features";
    DistancesPrediction = CalculateDmdf(T_TractogramPrediction, T_Prototypes);
    DistancesGroundtruth = CalculateDmdf(T_TractogramGroundtruth, T_Prototypes);
    DistancesTest = CalculateDmdf(T_TractogramTest, T_Prototypes);

    LabelsGroundtruth = CreateLabels(DistancesTest, DistancesGroundtruth);
    LabelsPrediction = CreateLabels(DistancesTest, DistancesPrediction);

    std::ofstream LabelsPredictionFile;
    LabelsPredictionFile.open("/home/r948e/mycsv/predictionlabels_" + std::to_string(m_activeCycle) + ".csv");
    for (unsigned int i = 0; i < LabelsPrediction.size(); i++)
    {
        LabelsPredictionFile << LabelsPrediction.at(i) << ' ';
    }
    LabelsPredictionFile.close();

    std::ofstream LabelsGroundtruthFile;
    LabelsGroundtruthFile.open("/home/r948e/mycsv/groundtruthlabels_" + std::to_string(m_activeCycle) + ".csv");
    for (unsigned int i = 0; i < LabelsGroundtruth.size(); i++)
    {
        LabelsGroundtruthFile << LabelsGroundtruth.at(i) << ' ';
    }
    LabelsGroundtruthFile.close();

    float FP = 0.0;
    float FN = 0.0;
    float TP = 0.0;
    float TN = 0.0;
    std::vector<int> indexfp;
//#pragma omp parallel for
    for (unsigned int i=0; i<LabelsGroundtruth.size(); i++)
    {
        if (LabelsGroundtruth.at(i)==1 && LabelsPrediction.at(i)==1)
        {
            TP +=1;
        }
        else if (LabelsGroundtruth.at(i)==1 && LabelsPrediction.at(i)==0)
        {
            FN +=1;
        }
        else if (LabelsGroundtruth.at(i)==0 && LabelsPrediction.at(i)==1)
        {
            FP +=1;
            indexfp.push_back(i);
        }
        else if (LabelsGroundtruth.at(i)==0 && LabelsPrediction.at(i)==0)
        {
            TN +=1;
        }
    }

    float Precision;
    float Recall;
    float F1_score;

    MITK_INFO << "TP";
    MITK_INFO << TP;
    MITK_INFO << "FP";
    MITK_INFO << FP;
    MITK_INFO << "TN";
    MITK_INFO << TN;
    MITK_INFO << "FN";
    MITK_INFO << FN;
    Precision = TP/(TP + FP);
    MITK_INFO << "Precision";
    MITK_INFO << Precision;
    Recall = TP/(TP + FN);
    MITK_INFO << "Recall";
    MITK_INFO << Recall;
    F1_score = (2*Precision*Recall)/(Precision+Recall);
    MITK_INFO << "F1_score";
    MITK_INFO << F1_score;



    vnl_vector<float> metrics(7);

    metrics.put(0, TP);
    metrics.put(1, FP);
    metrics.put(2, TN);
    metrics.put(3, FN);
    metrics.put(4, Precision);
    metrics.put(5, Recall);
    metrics.put(6, F1_score);




    return metrics;

}


}



