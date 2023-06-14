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

#include <itkTractDensityImageFilter.h>
#include <mitkImagePixelReadAccessor.h>
#include <boost/progress.hpp>
#include <vnl/vnl_sparse_matrix.h>
#include <mitkIOUtil.h>
#include <random>

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
  m_activeCycle = activeCycle;
}


void StreamlineFeatureExtractor::SetTractogramTest(const mitk::FiberBundle::Pointer &TractogramTest, std::string TractogramTestName)
{
    MITK_INFO << TractogramTestName;
    m_TractogramTest= TractogramTest;
}

std::vector<vnl_matrix<float> > StreamlineFeatureExtractor::ResampleFibers(mitk::FiberBundle::Pointer tractogram)
{
  MITK_INFO << "Infunction";
   mitk::FiberBundle::Pointer temp_fib = tractogram->GetDeepCopy();
   vtkPolyData* polyData = temp_fib->GetFiberPolyData();

   for (vtkIdType i = 0; i < polyData->GetNumberOfCells(); i++)
   {
     vtkCell* cell = polyData->GetCell(i);
     if (cell->GetNumberOfPoints() != polyData->GetCell(0)->GetNumberOfPoints())
     {
       throw std::runtime_error("Not all cells have an equal number of points!");
 //        temp_fib->ResampleToNumPoints(m_NumPoints);
     }
   }

  MITK_INFO << "Resampling Done";

  std::vector< vnl_matrix<float> > out_fib(tractogram->GetFiberPolyData()->GetNumberOfCells());//temp_fib
//#pragma omp parallel for
  for (int i=0; i<tractogram->GetFiberPolyData()->GetNumberOfCells(); i++)//temp_fib
  {
    vtkCell* cell = tractogram->GetFiberPolyData()->GetCell(i);//temp_fib
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

    out_fib.at(i)=streamline;
  }



  return out_fib;
}

std::vector<vnl_matrix<float> > StreamlineFeatureExtractor::CalculateDmdf(std::vector<vnl_matrix<float> > tractogram,
                                                                          std::vector<vnl_matrix<float> > prototypes)
{
// Initialize a vector to store distance matrices for each tractogram item
std::vector<vnl_matrix<float>> dist_vec(tractogram.size());

// Loop over each tractogram item
#pragma omp parallel for
for (std::size_t i = 0; i < tractogram.size(); ++i)
{
    // Initialize a distance matrix for the current tractogram item
    vnl_matrix<float> distances(1, 2 * prototypes.size(), 0.0);

    // Loop over each prototype
    for (std::size_t j = 0; j < prototypes.size(); ++j)
    {
        // Initialize matrices to store distances between the current tractogram item and the current prototype
        vnl_matrix<float> single_distances(1, tractogram[0].cols(), 0.0);
        vnl_matrix<float> single_distances_flip(1, tractogram[0].cols(), 0.0);
        vnl_matrix<float> single_end_distances(1, 2, 0.0);
        vnl_matrix<float> single_end_distances_flip(1, 2, 0.0);

        // Loop over each point in the current tractogram item
        for (std::size_t ik = 0; ik < tractogram[0].cols(); ++ik)
        {
            // Calculate the Euclidean distance between the current point in the current tractogram item and the current point in the current prototype
            const double cur_dist = std::sqrt(std::pow(tractogram[i](0, ik) - prototypes[j](0, ik), 2.0)
                                              + std::pow(tractogram[i](1, ik) - prototypes[j](1, ik), 2.0)
                                              + std::pow(tractogram[i](2, ik) - prototypes[j](2, ik), 2.0));
            single_distances(0, ik) = cur_dist;

            // Calculate the Euclidean distance between the current point in the current tractogram item and the mirrored current point in the current prototype
            const double cur_dist_flip = std::sqrt(std::pow(tractogram[i](0, ik) - prototypes[j](0, prototypes[0].cols() - (ik + 1)), 2.0)
                                                   + std::pow(tractogram[i](1, ik) - prototypes[j](1, prototypes[0].cols() - (ik + 1)), 2.0)
                                                   + std::pow(tractogram[i](2, ik) - prototypes[j](2, prototypes[0].cols() - (ik + 1)), 2.0));
            single_distances_flip(0, ik) = cur_dist_flip;

            // Store the distances between the first and last points of the current tractogram item and the current prototype
            if (ik == 0)
            {
                single_end_distances(0, 0) = cur_dist;
                single_end_distances_flip(0, 0) = cur_dist_flip;
            }
            else if (ik == tractogram[0].cols() - 1)
            {
                single_end_distances(0, 1) = cur_dist;
                single_end_distances_flip(0, 1) = cur_dist_flip;
            }
        }

        // Calculate the mean distance between the current tractogram item and the current prototype
        const double mean_single_distances = single_distances.mean();
        const double mean_single_distances_flip = single_distances_flip.mean();
        // distances(0, j) = mean_single_distances_flip > mean_single_distances ? mean_single_distances : mean_single_distances_flip;

        // // Calculate the mean END distance between the current tractogram item and the current prototype
        const double mean_single_end_distances = single_end_distances.mean();
        const double mean_single_end_distances_flip = single_end_distances_flip.mean();
        const std::size_t end_distances_index = prototypes.size() + j;
        distances(0, j) = mean_single_distances_flip > mean_single_distances ? mean_single_distances : mean_single_distances_flip;

        distances(0, end_distances_index) = mean_single_end_distances_flip > mean_single_end_distances ? mean_single_end_distances : mean_single_end_distances_flip;
    }

        dist_vec[i] = distances;
    }

    MITK_INFO << "The size of the distances is " <<dist_vec[0].cols();
    MITK_INFO << "The size of the distances is " <<dist_vec[0].rows();
    MITK_INFO << "Done calculating Dmdf";

    return dist_vec;
}

std::vector<vnl_matrix<float> > StreamlineFeatureExtractor::MergeTractogram(std::vector<vnl_matrix<float> > prototypes,
                                                                           std::vector<vnl_matrix<float> > positive_local_prototypes,
                                                                           std::vector<vnl_matrix<float> > negative_local_prototypes)
{
    unsigned int pos_locals;
    unsigned int neg_locals;

    if (positive_local_prototypes.size() >= 50)
    {
        pos_locals= 50;
    }
    else {
        pos_locals= positive_local_prototypes.size();
    }

    if (negative_local_prototypes.size() >= 50)
    {
        neg_locals = 50;
    }
    else {
        neg_locals= negative_local_prototypes.size();
    }



    std::vector< vnl_matrix<float> > merged_prototypes;

    for (unsigned int k=0; k<prototypes.size(); k++)
    {
        merged_prototypes.push_back(prototypes.at(k));
    }

    for (unsigned int k=0; k<negative_local_prototypes.size(); k++)
    {
        merged_prototypes.push_back(negative_local_prototypes.at(k));
    }


    for (unsigned int k=0; k<positive_local_prototypes.size(); k++)
    {
        merged_prototypes.push_back(positive_local_prototypes.at(k));
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
    MITK_INFO << "Calculate Plus Features";
    m_DistancesPlus = CalculateDmdf(T_TractogramPlus, T_mergedPrototypes);

    MITK_INFO << "Resample Test Data";
    T_TractogramTest= ResampleFibers(m_TractogramTest);
    MITK_INFO << "Calculate Features of Test Data";
    m_DistancesTest= CalculateDmdf(T_TractogramTest, T_mergedPrototypes);
    MITK_INFO << m_DistancesTest.size();

    // if (m_activeCycle == 0){
    myindex = GetIndex(m_DistancesTest);    
    // }
    
    




    TrainModel();
    // m_index = PredictLabels();

    m_index = Predict3();

}
std::vector<unsigned int> StreamlineFeatureExtractor::GetIndex(std::vector< vnl_matrix<float>> distances)
{
    std::vector<unsigned int> indices;
    unsigned int num_matrices = distances.size();

    // Calculate the threshold based on the mean of the last row in the last matrix
    float threshold = 0.0;
    if (num_matrices > 0) {
        const vnl_matrix<float>& lastMatrix = distances.back();
        unsigned int lastRow = lastMatrix.rows() - 1;

        std::vector<float> values;
        for (unsigned int col = 0; col < lastMatrix.cols(); ++col) {
            values.push_back(lastMatrix(lastRow, col));
        }

        // Only take 20% of data for higher speed
        std::sort(values.begin(), values.end());
        unsigned int firstQuartileIndex = std::floor(values.size() / 5.0);
        threshold = values[firstQuartileIndex];
    }

    // Check each matrix to find indices that meet the condition
    for (unsigned int i = 0; i < num_matrices; ++i) {
        const vnl_matrix<float>& matrix = distances[i];
        unsigned int num_cols = matrix.cols();
        unsigned int last_col_index = num_cols - 1;
        bool last_col_below_threshold = true;

        for (unsigned int row = 0; row < matrix.rows(); ++row) {
            if (matrix(row, last_col_index) >= threshold) {
                last_col_below_threshold = false;
                break;
            }
        }

        if (last_col_below_threshold) {
            indices.push_back(i);
        }
    }

    // Print debug information
    std::cout << "Number of matrices: " << num_matrices << std::endl;
    std::cout << "Number of indices: " << indices.size() << std::endl;

    return indices;
}
//std::vector<unsigned int> StreamlineFeatureExtractor::GetIndex(std::vector< vnl_matrix<float>> distances)
//{
//    std::vector<unsigned int> indices;
//    unsigned int num_matrices = distances.size();


//    // Calculate the mean of the last row in the distances matrix
//    float threshold = 0.0;
//    if (num_matrices > 0) {
//        const vnl_matrix<float>& lastMatrix = distances[num_matrices - 1];
//        unsigned int lastRow = lastMatrix.rows() - 1;

//        std::vector<float> values;
//        for (unsigned int col = 0; col < lastMatrix.cols(); ++col) {
//            values.push_back(lastMatrix(lastRow, col));
//        }

//        std::sort(values.begin(), values.end());
//        unsigned int firstQuartileIndex = std::floor(values.size() / 4.0);
//        threshold = values[firstQuartileIndex];
//    }

//    // for (unsigned int i = 0; i < num_matrices; ++i) {

//    //     const vnl_matrix<float>& matrix = distances[i];
//    //     unsigned int num_cols = matrix.cols();
//    //     unsigned int last_col_index = num_cols - 1;
//    //     bool last_col_below_10 = true;
//    //     if (matrix(0, last_col_index) >= 25) {
//    //         last_col_below_10 = false;
//    //         // break;
//    //     }

//    //     if (last_col_below_10) {
//    //         indices.push_back(i);
//    //     }
//    // }


//    for (unsigned int i = 0; i < num_matrices; ++i) {
//        const vnl_matrix<float>& matrix = distances[i];
//        unsigned int num_cols = matrix.cols();
//        unsigned int last_col_index = num_cols - 1;
//        bool last_col_below_threshold = true;

//        for (unsigned int row = 0; row < matrix.rows(); ++row) {
//            if (matrix(row, last_col_index) >= threshold) {
//                last_col_below_threshold = false;
//                break;
//            }
//        }

//        if (last_col_below_threshold) {
//            indices.push_back(i);
//        }
//    }
//    std::cout << "Number of matrices: " << num_matrices << std::endl;
//    std::cout << "Number of indices: " << indices.size() << std::endl;
//    return indices;
//}

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

void StreamlineFeatureExtractor::TrainModel()
{
    MITK_INFO << "Start Function Get Data";


    cv::Mat data;
    cv::Mat labels_arr_vec;

    int size_plus = 0;


    /*Create Trainingdata: Go through positive and negative Bundle and save distances as cv::Mat and create vector with labels*/
    for ( unsigned int i=0; i<m_DistancesPlus.size(); i++)
    {
        float data_arr [m_DistancesPlus.at(0).size()];

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

    /*Create Dataset and initialize Classifier*/
    cv::Ptr<cv::ml::TrainData> m_traindata = cv::ml::TrainData::create(samples_shuffled, cv::ml::ROW_SAMPLE, labels_shuffled);




    statistic_model= cv::ml::RTrees::create();
    auto criteria = cv::TermCriteria();
    criteria.type = cv::TermCriteria::MAX_ITER;
//    criteria.epsilon = 1e-8;
    criteria.maxCount = 300;

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

    // // Create a new MITK data node and set its data to the RTrees instance
    // mitk::DataNode::Pointer newNode = mitk::DataNode::New();
    // newNode->SetName("My RTrees model");
    // newNode->SetData(rtrees);

    // // Add the new data node to the MITK data storage
    // this->GetDataStorage()->Add(newNode);

}

std::vector<std::vector<unsigned int>>  StreamlineFeatureExtractor::Predict()
{
     std::vector<std::vector<unsigned int>> index_vec;


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
  if (lengths>500)
  {
      lengths=500;
  }
  else
  {
    lengths=100;
  }



    // int lengthsCertain = std::count_if(e.begin(), e.end(),[&](auto const& val){ return val < 0.1; });

    std::vector<unsigned int> indexUnc = Sort(e, lengths, 0);

    // std::vector<unsigned int> indexCertain = Sort(e, e.size() , e.size()-lengthsCertain );

//    std::vector<unsigned int> indexCertainBetween = Sort(e, e.size()-lengthsCertain , lengths);

    // MITK_INFO << "Index Certainty Vector size";
    // MITK_INFO << indexCertain.size();

    // std::vector<unsigned int> indexCertainNeg;
    // std::vector<unsigned int> indexCertainPos;

    // for (unsigned int i=0; i<indexCertain.size(); i++)
    // {
    //     if(pred.at(indexCertain.at(i))==0 )
    //     {
    //         indexCertainNeg.push_back(indexCertain.at(i));
    //     }
    //     else {
    //         indexCertainPos.push_back(indexCertain.at(i));
    //     }

    // }



    MITK_INFO << "Index Uncertainty Vector size";
    MITK_INFO << indexUnc.size();
    // MITK_INFO << "Index Certainty Vector size";
    // MITK_INFO << indexCertainPos.size();
    // MITK_INFO << indexCertainNeg.size();
    // MITK_INFO << indexCertainNeg.size() +indexCertainPos.size();
    MITK_INFO << "Index Certainty between Vector size";
    index_vec.push_back(indexPrediction);
    index_vec.push_back(indexUnc);
    // index_vec.push_back(indexCertainNeg);

    return index_vec;

}

std::vector<std::vector<unsigned int>>  StreamlineFeatureExtractor::Predict2()
{

    std::vector<std::vector<unsigned int>> index_vec;
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

    std::vector<float> raw_scores;
    cv::Mat votes;

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

    MITK_INFO << "--------------";
    MITK_INFO << "Prediction vector size:";
    MITK_INFO << indexPrediction.size();
    MITK_INFO << "Entropy vector size:";
    entropy_vector = e;
    MITK_INFO << e.size();

    MITK_INFO << "--------------";

    /*Get index of most unertain data (lengths defines how many data is saved)*/

      int lengths = std::count_if(e.begin(), e.end(),[&](auto const& val){ return val >= 0.75; });
    if (lengths>1500)
    {
        lengths=1500;
    }


    std::vector<unsigned int> indexUnc = Sort(e, lengths, 0);

    mitk::FiberBundle::Pointer fib =  CreatePrediction(indexUnc, false);

    typedef unsigned char OutPixType;
    typedef itk::Image<OutPixType, 3> OutImageType;

    itk::TractDensityImageFilter< OutImageType >::Pointer generator = itk::TractDensityImageFilter< OutImageType >::New();
    generator->SetFiberBundle(fib);
    generator->SetMode(TDI_MODE::VISITATION_COUNT);
    generator->SetOutputAbsoluteValues(true);
    generator->SetUpsamplingFactor(0.5);

    generator->Update();


    // get output image
    typedef itk::Image<OutPixType,3> OutType;
    OutType::Pointer outImg = generator->GetOutput();
    mitk::Image::Pointer img = mitk::Image::New();
    img->InitializeByItk(outImg.GetPointer());
    img->SetVolume(outImg->GetBufferPointer());


    vtkSmartPointer<vtkPolyData> fiberPolyData = fib->GetFiberPolyData();
    mitk::BaseGeometry::Pointer geometry = fib->GetGeometry();


    vtkSmartPointer<vtkFloatArray> weights = vtkSmartPointer<vtkFloatArray>::New();

    std::vector<float> intensities_vec;
    float intensitiesMax = 0;
    for (unsigned int i=0; i<fiberPolyData->GetNumberOfCells(); i++)
    {

        vtkCell* cell = m_TractogramTest->GetFiberPolyData()->GetCell(indexUnc.at(i));
        auto numPoints = cell->GetNumberOfPoints();
        vtkPoints* points = cell->GetPoints();

        int intensities = 0;
        for (unsigned int j=0; j<numPoints; j++)
        {
            double p[3];
            points->GetPoint(j, p);
            mitk::Point3D point; // define the point where you want to get the pixel value
            point[0] = p[0];       // set the x coordinate
            point[1] = p[1];       // set the y coordinate
            point[2] = p[2];       // set the z coordinate
            mitk::ImagePixelReadAccessor<OutPixType, 3> readAccess(img);
            float value(readAccess.GetPixelByWorldCoordinates(point));
            intensities += value;

        }
        intensities_vec.push_back(intensities);
            if (intensitiesMax<intensities){
                intensitiesMax = intensities;
            }
        

    }
    unsigned int counter = 0;
    
    std::vector<float> e_mix;
    for (unsigned int i=0; i<indexUnc.size(); i++){
        float val;
        val = 1-intensities_vec.at(i)/ intensitiesMax;
        e_mix.push_back(e.at(indexUnc.at(i)) * val);
        MITK_INFO << "E_mix: " << e_mix.at(i) << "E: " << e.at(indexUnc.at(i)) << "Densitiy " << val;
        MITK_INFO << val; 
        weights->InsertValue(counter, e_mix.at(i));
        counter++;
    }

    
    std::vector<unsigned int> indexmix = Sort(e_mix, lengths, 0);


    std::vector<unsigned int> newidx;

    for (unsigned int i = 0; i < indexmix.size(); i++) {
        newidx.push_back(indexUnc.at(indexmix.at(i)));
    }

    fib->SetFiberWeights(weights);
    fib->ColorFibersByFiberWeights(false, mitk::LookupTable::JET_TRANSPARENT);

    
    
    m_fiberNode = mitk::DataNode::New();
    m_fiberNode->SetData(fib);
    m_fiberNode->SetName("Colorbydensity");


     // Get the dimensions of the density image
    // typename OutputImageType::SizeType imageSize = outImg->GetLargestPossibleRegion().GetSize();
    // int numVoxels = imageSize[0] * imageSize[1] * imageSize[2];

    // // init data node
    m_imgNode = mitk::DataNode::New();
    m_imgNode->SetData(img);
    m_imgNode->SetName("Density");
    // ... initialize and set the image ...

    // mitk::Point3D point; // define the point where you want to get the pixel value
    // point[0] = 0;       // set the x coordinate
    // point[1] = 0;       // set the y coordinate
    // point[2] = 0;       // set the z coordinate
    // mitk::Point3D index;
    // img->GetGeometry()->WorldToIndex(point, index);
    // MITK_INFO << index;

    // get the pixel value at the given point
    // mitk::ImagePixelReadAccessor<OutPixType, 3> readAccess(img);
    // float value(readAccess.GetPixelByWorldCoordinates(point));



    index_vec.push_back(indexPrediction);
    index_vec.push_back(newidx);
    

    return index_vec;

}

std::vector<std::vector<unsigned int>>  StreamlineFeatureExtractor::Predict3()
{
     std::vector<std::vector<unsigned int>> index_vec;


    /*Predict on Test Data*/
    MITK_INFO << "Predicting";

    /*Create Dataset as cv::Mat*/
    cv::Mat dataTest;
    for ( unsigned int i=0; i<myindex.size(); i++)
    {
        float data_arr [m_DistancesTest.at(0).size()];

        for ( unsigned int j=0; j<m_DistancesTest.at(myindex.at(0)).cols(); j++)
        {
            data_arr[j] = m_DistancesTest.at(myindex.at(i)).get(0,j);
        }
        cv::Mat curdata(1, m_DistancesTest.at(myindex.at(0)).size(), CV_32F, data_arr);
        dataTest.push_back(curdata);
    }


    std::vector<unsigned int> indexPrediction;
    std::vector<float> e(myindex.size());
    // std::fill(e.begin(), e.end(), 0.0f);
    std::vector<int> pred(myindex.size());
    std::fill(pred.begin(), pred.end(), 0.0f);


    /*For every Sample/Streamline get Prediction and entropy (=based on counts of Random Forest)*/
    MITK_INFO << "Predicting on all cores";
#pragma omp parallel for
    for (unsigned int i=0; i<myindex.size(); i++)
    {


        int val = statistic_model->predict(dataTest.row(i));
        pred.at(i)=val;

        #pragma omp critical
        if (val==1)
        {
           indexPrediction.push_back(myindex.at(i));
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

    MITK_INFO << "--------------";
    MITK_INFO << "Prediction vector size:";
    MITK_INFO << indexPrediction.size();
    MITK_INFO << "Entropy vector size:";
    entropy_vector = e;
    MITK_INFO << e.size();

    MITK_INFO << "--------------";


    int lengths = std::count_if(e.begin(), e.end(),[&](auto const& val){ return val >= 0.95; });
    if (lengths>500)
    {
        lengths=500;
    }
    else
    {
        lengths=100;
    }


    std::vector<unsigned int> indexUnc = Sort(e, lengths, 0);

    std::vector<unsigned int> newidx;
    
    for (unsigned int i = 0; i < indexUnc.size(); i++) {
        newidx.push_back(myindex.at(indexUnc.at(i)));
    }

    MITK_INFO << indexUnc.size();
    index_vec.push_back(indexPrediction);
    index_vec.push_back(newidx);

    return index_vec;

}

std::vector<std::vector<unsigned int>>  StreamlineFeatureExtractor::GetDistanceData(float &value)
{

    // /*Vector which saves Fibers to be labeled based on fft subset uncertainty*/
    // std::vector<std::vector<unsigned int>> index_vec;

    // /*Get index of most unertain data (lengths defines how many data is saved)*/
    // MITK_INFO << entropy_vector.size();
    // int lengths = std::count_if(entropy_vector.begin(), entropy_vector.end(), [&](auto const& val){ return val >= value; });
    // if (lengths>1500)
    // {
    //     lengths=1500;
    // }
    // MITK_INFO << lengths;
    // /*Maybe shuffling of length so not the most uncertain values are chosen*/
    // std::vector<unsigned int> indexUnc = Sort(entropy_vector, lengths, 0);

    // vnl_matrix<float> distances_matrix;

    // distances_matrix.set_size(lengths, lengths);
    // distances_matrix.fill(0.0);

    // std::vector<float> distances_matrix_mean;

    // #pragma omp parallel for schedule(dynamic)
    // for (int i=0; i<lengths; i++)
    // {
    //     for (int k=0; k<lengths; k++)
    //     {
    //         /*From the length.size() Samples with the highest Entropey calculate the differen between the Features*/
    //         vnl_matrix<float> diff =  m_DistancesTest.at(indexUnc.at(i)) - m_DistancesTest.at(indexUnc.at(k));

    //         /*Into the eucledean difference matrix, put the distance in Feature Space between every sample pair*/
    //         distances_matrix.put(i,k,diff.absolute_value_sum()/m_DistancesTest.at(0).size());

    //     }
    //     /*For every Sample/Streamline get the mean eucledean distance to all other Samples => one value for every Sample*/
    //     //distances_matrix_mean.push_back(distances_matrix.get_row(i).mean());
    //     //MITK_INFO << meanval.at(i);

    // }

    // MITK_INFO << "Distance Matrix Calculated";

    // /*Index to find values in distancematrix*/
    // std::vector<unsigned int> myidx;
    // /*Index to find actual streamlines using indexUnc*/
    // std::vector<unsigned int> indexUncDist;
    // /*Start with the Streamline of the highest entropy, which is in distance_matrix at idx 0*/
    // myidx.push_back(0);
    // indexUncDist.push_back(indexUnc.at(myidx.at(0)));

    // /*Vecotr that stores minvalues of current iteration*/
    // vnl_matrix<float> cur_vec;
    // cur_vec.set_size(1,lengths);
    // cur_vec.fill(0.0);

    // // #pragma omp parallel for schedule(dynamic)
    // for (int i=0; i<lengths; i++)
    // {

    //     /*Save mean distance of all used Samples*/
    //     vnl_matrix<float> sum_matrix;
    //     sum_matrix.set_size(myidx.size(), lengths);
    //     sum_matrix.fill(0);

    //     // #pragma omp parallel for schedule(dynamic)
    //     for (unsigned int ii=0; ii<myidx.size(); ii++)
    //     {
    //         sum_matrix.set_row(ii, distances_matrix.get_column(myidx.at(ii)));
    //     }

    //     for (unsigned int k=0; k<sum_matrix.columns(); k++)
    //     {
    //         cur_vec.put(0,k, sum_matrix.get_column(k).min_value());
    //     }
    //     myidx.push_back(cur_vec.arg_max());

    //     indexUncDist.push_back(indexUnc.at(myidx.at(i+1)));
    //     sum_matrix.clear();

    // }

    // MITK_INFO << "Dist_stop";

    // index_vec.push_back(indexUncDist);


    /*Vector which saves Fibers to be labeled based on fft subset uncertainty*/
    std::vector<std::vector<unsigned int>> index_vec;

    /*Get index of most unertain data (lengths defines how many data is saved)*/
//    int lengths=500;
    MITK_INFO << entropy_vector.size();
  int lengths = std::count_if(entropy_vector.begin(), entropy_vector.end(),[&](auto const& val){ return val >= value; });
  if (lengths>500)
  {
      lengths=500;
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

std::vector<std::vector<unsigned int>>  StreamlineFeatureExtractor::GetDistanceData2(float &value)
{

    /*Vector which saves Fibers to be labeled based on fft subset uncertainty*/
    std::vector<std::vector<unsigned int>> index_vec;

    /*Get index of most unertain data (lengths defines how many data is saved)*/
//    int lengths=500;
    MITK_INFO << entropy_vector.size();
  int lengths = std::count_if(entropy_vector.begin(), entropy_vector.end(),[&](auto const& val){ return val >= value; });
  if (lengths>250)
  {
      lengths=250;
  }
  MITK_INFO << lengths;

    /*Maybe shuffling of length so not the most uncertain values are chosen*/
    std::vector<unsigned int> indexUnc = Sort(entropy_vector, lengths, 0);

    std::vector<unsigned int> newidx;
    for (unsigned int i = 0; i < indexUnc.size(); i++) {
        newidx.push_back(myindex.at(indexUnc.at(i)));
    }

    vnl_matrix<float> distances_matrix;

    distances_matrix.set_size(lengths, lengths);
    distances_matrix.fill(0.0);

    std::vector<float> distances_matrix_mean;


    for (int i=0; i<lengths; i++)
    {
        for (int k=0; k<lengths; k++)
        {
            /*From the length.size() Samples with the highest Entropey calculate the differen between the Features*/
            vnl_matrix<float> diff =  m_DistancesTest.at(newidx.at(i)) - m_DistancesTest.at(newidx.at(k));

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
    indexUncDist.push_back(newidx.at(myidx.at(0)));

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


        indexUncDist.push_back(newidx.at(myidx.at(i+1)));
        sum_matrix.clear();

    }

    MITK_INFO << "Dist_stop";

    index_vec.push_back(indexUncDist);

    return index_vec;
}


}





 
