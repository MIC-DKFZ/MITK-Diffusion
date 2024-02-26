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
#include <boost/timer/progress_display.hpp>
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

void StreamlineFeatureExtractor::SetTractogramPrototypes(const mitk::FiberBundle::Pointer &TractogramPrototypes)
{
      MITK_INFO << "Use individual Prototypes...";
      m_inputPrototypes = TractogramPrototypes;

}

void StreamlineFeatureExtractor::SetActiveCycle(int &activeCycle)
{
  m_activeCycle = activeCycle;
}

void StreamlineFeatureExtractor::SetTractogramTest(const mitk::FiberBundle::Pointer &TractogramTest)
{
    m_TractogramTest= TractogramTest;
}

std::vector<vnl_matrix<float> > StreamlineFeatureExtractor::TractToMatrixFibers(mitk::FiberBundle::Pointer tractogram)
{

  // Create a deep copy of the input tractogram
  mitk::FiberBundle::Pointer temp_fib = tractogram->GetDeepCopy();

  // Get the fiber polydata associated with the copied tractogram
  vtkPolyData* polyData = temp_fib->GetFiberPolyData();

  // Check if all cells in the fiber polydata have the same number of points
  for (vtkIdType i = 0; i < polyData->GetNumberOfCells(); i++)
  {
    vtkCell* cell = polyData->GetCell(i);

    // If the number of points in the current cell is different from the first cell,
    // throw a runtime error indicating the need for preprocessing
    if (cell->GetNumberOfPoints() != polyData->GetCell(0)->GetNumberOfPoints())
    {

      throw std::runtime_error("Not all cells have an equal number of points! Resample the tractogram first in the preprocessing view and reset the Classifier before restart");
    }
  }

  // Log a message to indicate that the resampling is done
  MITK_INFO << "All fibers have equal number of points";

  // Create a vector to store the output matrices, one for each fiber
  std::vector< vnl_matrix<float> > out_fib(tractogram->GetFiberPolyData()->GetNumberOfCells());

  // Process each cell/fiber in the fiber polydata
  for (int i = 0; i < tractogram->GetFiberPolyData()->GetNumberOfCells(); i++)
  {
    // Get the current cell/fiber
    vtkCell* cell = tractogram->GetFiberPolyData()->GetCell(i);

    // Get the number of points in the current fiber
    int numPoints = cell->GetNumberOfPoints();

    // Get the points of the current fiber
    vtkPoints* points = cell->GetPoints();

    // Create a matrix to store the coordinates of the fiber points
    vnl_matrix<float> streamline;
    streamline.set_size(3, m_NumPoints); // Set the size of the matrix to 3 rows and m_NumPoints columns
    streamline.fill(0.0); // Initialize all matrix elements to zero

    // Process each point in the fiber
    for (int j = 0; j < numPoints; j++)
    {
      double cand[3];
      points->GetPoint(j, cand);

      // Create a vnl_vector_fixed object to store the coordinates of the current point
      vnl_vector< float > candV(3);
      candV[0] = cand[0];
      candV[1] = cand[1];
      candV[2] = cand[2];

      // Set the j-th column of the streamline matrix to the current point's coordinates
      streamline.set_column(j, candV);
    }

    // Store the completed streamline matrix in the output vector
    out_fib.at(i) = streamline;
  }

  // Return the vector of streamline matrices
  return out_fib;
}

std::vector<vnl_matrix<float> > StreamlineFeatureExtractor::CalculateDmdf(std::vector<vnl_matrix<float> > tractogram,
                                                                          std::vector<vnl_matrix<float> > prototypes)
{
// Initialize a vector to store distance matrices for each tractogram item
std::vector<vnl_matrix<float>> dist_vec(tractogram.size());

// Loop over each tractogram item
#pragma omp parallel for
for (int i = 0; i < static_cast<int>(tractogram.size()); ++i)
{
    // Initialize a distance matrix for the current tractogram item
    vnl_matrix<float> distances(1, 2 * prototypes.size(), 0.0);

    // Loop over each prototype
    for (int j = 0; j < static_cast<int>(prototypes.size()); ++j)
    {
        // Initialize matrices to store distances between the current tractogram item and the current prototype
        vnl_matrix<float> single_distances(1, tractogram[0].cols(), 0.0);
        vnl_matrix<float> single_distances_flip(1, tractogram[0].cols(), 0.0);
        vnl_matrix<float> single_end_distances(1, 2, 0.0);
        vnl_matrix<float> single_end_distances_flip(1, 2, 0.0);

        // Loop over each point in the current tractogram item
        for (int ik = 0; ik < static_cast<int>(tractogram[0].cols()); ++ik)
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
            else if (ik == static_cast<int>(tractogram[0].cols()) - 1)
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


    return dist_vec;
}

std::vector<vnl_matrix<float> > StreamlineFeatureExtractor::MergeTractogram(std::vector<vnl_matrix<float> > prototypes,
                                                                           std::vector<vnl_matrix<float> > positive_local_prototypes,
                                                                           std::vector<vnl_matrix<float> > negative_local_prototypes)
{
    // Variables to store the number of positive and negative local prototypes
    unsigned int pos_locals;
    unsigned int neg_locals;

    // Determine the number of positive local prototypes to consider
    if (positive_local_prototypes.size() >= 50)
    {
        pos_locals = 50;
    }
    else
    {
        pos_locals = positive_local_prototypes.size();
    }

    // Determine the number of negative local prototypes to consider
    if (negative_local_prototypes.size() >= 50)
    {
        neg_locals = 50;
    }
    else
    {
        neg_locals = negative_local_prototypes.size();
    }

    // Create a vector to store the merged prototypes
    std::vector< vnl_matrix<float> > merged_prototypes;

    // Copy the prototypes into the merged prototypes vector
    for (unsigned int k = 0; k < prototypes.size(); k++)
    {
        merged_prototypes.push_back(prototypes.at(k));
    }

    // Append the negative local prototypes to the merged prototypes vector
    for (unsigned int k = 0; k < negative_local_prototypes.size(); k++)
    {
        merged_prototypes.push_back(negative_local_prototypes.at(k));
    }

    // Append the positive local prototypes to the merged prototypes vector
    for (unsigned int k = 0; k < positive_local_prototypes.size(); k++)
    {
        merged_prototypes.push_back(positive_local_prototypes.at(k));
    }

    // Log the number of prototypes, positive local prototypes, and negative local prototypes
    MITK_INFO << "Number of prototypes:";
    MITK_INFO << prototypes.size();
    MITK_INFO << "Number of positive local prototypes:";
    MITK_INFO << pos_locals;
    MITK_INFO << "Number of negative local prototypes:";
    MITK_INFO << neg_locals;

    // Return the merged prototypes vector
    return merged_prototypes;
}

std::vector<unsigned int> StreamlineFeatureExtractor::Sort(std::vector<float> sortingVector, int lengths, int start)
{
    // Create a vector to store the sorted indices
    std::vector<unsigned int> index;

    // Create a priority queue to sort the elements in descending order
    std::priority_queue<std::pair<float, int>> q;

    // Push each element of the sortingVector along with its index to the priority queue
    for (unsigned int i = 0; i < sortingVector.size(); ++i)
    {
        q.push(std::pair<float, int>(sortingVector[i], i));
    }

    // Retrieve the top 'lengths' elements from the priority queue
    // and store their indices in the 'index' vector
    for (int i = 0; i < lengths; ++i)
    {
        int ki = q.top().second;

        // Check if the current index is greater than or equal to the 'start' index
        // and only then add it to the 'index' vector
        if (i >= start)
        {
            index.push_back(ki);
        }

        // Remove the top element from the priority queue
        q.pop();
    }

    // Return the sorted indices
    return index;
}

mitk::FiberBundle::Pointer StreamlineFeatureExtractor::CreatePrediction(std::vector<unsigned int>& index, bool removefrompool)
{

    // Create a new vtkPolyData for storing the fibers
    vtkSmartPointer<vtkPolyData> fibersData = vtkSmartPointer<vtkPolyData>::New();
    fibersData->SetPoints(vtkSmartPointer<vtkPoints>::New());
    fibersData->SetLines(vtkSmartPointer<vtkCellArray>::New());

    // Create temporary variables for fiber construction
    vtkSmartPointer<vtkPolyData> vNewPolyData = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkCellArray> vNewLines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> vNewPoints = vtkSmartPointer<vtkPoints>::New();

    unsigned int indexSize = index.size();
    unsigned int counter = 0;

    for (unsigned int i = 0; i < indexSize; i++)
    {
        // Get the cell (fiber) from the tractogram
        vtkCell* cell = m_TractogramTest->GetFiberPolyData()->GetCell(index[i]);
        auto numPoints = cell->GetNumberOfPoints();
        vtkPoints* points = cell->GetPoints();

        // Create a new polyline for the current fiber
        vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
        for (unsigned int j = 0; j < numPoints; j++)
        {
            double p[3];
            points->GetPoint(j, p);

            // Add fiber points to the new points array
            vtkIdType id = vNewPoints->InsertNextPoint(p);
            container->GetPointIds()->InsertNextId(id);
        }

        // Add the polyline (fiber) to the new lines array
        vNewLines->InsertNextCell(container);
        counter++;
    }

    if (removefrompool)
    {
        // Remove the selected cells from the tractogram
        for (unsigned int i = 0; i < indexSize; i++)
        {
            m_TractogramTest->GetFiberPolyData()->DeleteCell(index[i]);
        }
    }


    // Set the lines and points to the new vtkPolyData
    vNewPolyData->SetLines(vNewLines);
    vNewPolyData->SetPoints(vNewPoints);

    // Update the fibersData with the new vtkPolyData
    fibersData->SetPoints(vNewPoints);
    fibersData->SetLines(vNewLines);

    // Create a new FiberBundle using the fibersData
    mitk::FiberBundle::Pointer prediction = mitk::FiberBundle::New(fibersData);

    return prediction;
}

void StreamlineFeatureExtractor::GenerateData()
{

    // Initialize variables
    std::vector<vnl_matrix<float>> T_Prototypes;
    std::vector<vnl_matrix<float>> T_TractogramPlus;
    std::vector<vnl_matrix<float>> T_TractogramMinus;
    std::vector<vnl_matrix<float>> T_TractogramTest;
    std::vector<vnl_matrix<float>> T_mergedPrototypes;

    // Convert input prototypes to matrices
    T_Prototypes = TractToMatrixFibers(m_inputPrototypes);

    // Convert input tractogram minus to matrices
    T_TractogramMinus = TractToMatrixFibers(m_TractogramMinus);

    // Convert input tractogram plus to matrices
    T_TractogramPlus = TractToMatrixFibers(m_TractogramPlus);

    // Merge prototypes with tractogram plus and minus
    T_mergedPrototypes = MergeTractogram(T_Prototypes, T_TractogramPlus, T_TractogramMinus);

    // Calculate features for tractogram minus
    MITK_INFO << "Calculate Features of Training Data";
    m_DistancesMinus = CalculateDmdf(T_TractogramMinus, T_mergedPrototypes);

    // Calculate features for tractogram plus
    m_DistancesPlus = CalculateDmdf(T_TractogramPlus, T_mergedPrototypes);

    // Convert test tractogram to matrices
    T_TractogramTest = TractToMatrixFibers(m_TractogramTest);

    // Calculate features for test data
    MITK_INFO << "Calculate Features of Test Data";
    m_DistancesTest = CalculateDmdf(T_TractogramTest, T_mergedPrototypes);

    // Get indices for prediction
    myindex = GetIndex(m_DistancesTest);

    // Train the model
    TrainModel();

    // Perform prediction
    m_index = Predict();
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
        if (distances.size()>80000 && m_activeCycle == 0){
            // Only take 20% of data for higher speed
            std::sort(values.begin(), values.end());
            unsigned int firstQuartileIndex = std::floor(values.size() / 5.0);
            threshold = values[firstQuartileIndex];
        }
        else {
            std::sort(values.begin(), values.end());
            threshold = values.back();
        }
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


    return indices;
}

//void StreamlineFeatureExtractor::TrainModel()
//{


//    cv::Mat data;
//    cv::Mat labels_arr_vec;

//    int size_plus = 0;


//    /*Create Trainingdata: Go through positive and negative Bundle and save distances as cv::Mat and create vector with labels*/
//    for ( unsigned int i=0; i<m_DistancesPlus.size(); i++)
//    {
//        float data_arr [m_DistancesPlus.at(0).size()];

//        labels_arr_vec.push_back(1);

//        for ( unsigned int j=0; j<m_DistancesPlus.at(0).cols(); j++)
//        {
//            data_arr[j] = m_DistancesPlus.at(i).get(0,j);
//        }
//        cv::Mat curdata(1, m_DistancesPlus.at(0).size(), CV_32F, data_arr);
//        data.push_back(curdata);
//        size_plus++;

//    }

//    for ( unsigned int i=m_DistancesPlus.size(); i<m_DistancesPlus.size()+m_DistancesMinus.size(); i++)
//    {
//        int it = i - size_plus;
//        float data_arr [m_DistancesMinus.at(0).size()];

//        labels_arr_vec.push_back(0);

//         for ( unsigned int j=0; j<m_DistancesMinus.at(0).cols(); j++)
//        {
//            data_arr[j] = m_DistancesMinus.at(it).get(0,j);
//        }
//        cv::Mat curdata(1, m_DistancesPlus.at(0).size(), CV_32F, data_arr);
//        data.push_back(curdata);

//    }

//    /*Calculate weights*/
//    int zerosgt = labels_arr_vec.rows - cv::countNonZero(labels_arr_vec);
//    int onesgt = cv::countNonZero(labels_arr_vec);
//    float plusval = labels_arr_vec.rows / (2.0 * onesgt ) ;
//    float minusval = labels_arr_vec.rows / (2.0 * zerosgt );
//    float w[2] = {minusval, plusval};
//    cv::Mat newweight = cv::Mat(1,2, CV_32F, w);


//    /*Shuffle Data*/
//    std::vector <int> seeds;
//    for (int cont = 0; cont < labels_arr_vec.rows; cont++)
//    {
//        seeds.push_back(cont);
//    }

//    cv::randShuffle(seeds);
//    cv::Mat labels_shuffled;
//    cv::Mat samples_shuffled;




//    for (int cont = 0; cont < labels_arr_vec.rows; cont++)
//    {
//        labels_shuffled.push_back(labels_arr_vec.row(seeds[cont]));
//    }

//    for (int cont = 0; cont < labels_arr_vec.rows; cont++)
//    {
//        samples_shuffled.push_back(data.row(seeds[cont]));
//    }

//    /*Create Dataset and initialize Classifier*/
//    cv::Ptr<cv::ml::TrainData> m_traindata = cv::ml::TrainData::create(samples_shuffled, cv::ml::ROW_SAMPLE, labels_shuffled);




//    statistic_model= cv::ml::RTrees::create();
//    auto criteria = cv::TermCriteria();
//    criteria.type = cv::TermCriteria::MAX_ITER;
////    criteria.epsilon = 1e-8;
//    criteria.maxCount = 300;

//    statistic_model->setMaxDepth(10); //set to three
////    statistic_model->setMinSampleCount(m_traindata->getNTrainSamples()*0.01);
//    statistic_model->setMinSampleCount(2);
//    statistic_model->setTruncatePrunedTree(false);
//    statistic_model->setUse1SERule(false);
//    statistic_model->setUseSurrogates(false);
//    statistic_model->setTermCriteria(criteria);
//    statistic_model->setCVFolds(1);
//    statistic_model->setPriors(newweight);


//    /*Train Classifier*/
//    MITK_INFO << "Start Training";
//    statistic_model->train(m_traindata);


//}

void StreamlineFeatureExtractor::TrainModel()
{
     cv::Mat data;
     cv::Mat labels_arr_vec;

     int size_plus = 0;

     /* Determine minority and majority classes */
     int minority_class = 0;  // Initialize with one of your class labels
     int majority_class = 1;  // Initialize with the other class label

    /* Create Trainingdata: Go through positive and negative Bundle and save distances as cv::Mat and create vector with labels */
    for (int i = 0; i < static_cast<int>(m_DistancesPlus.size()); i++)
    {
        float data_arr[m_DistancesPlus.at(0).size()];

        labels_arr_vec.push_back(1);

        for (int j = 0; j < static_cast<int>(m_DistancesPlus.at(0).cols()); j++)
        {
            data_arr[j] = m_DistancesPlus.at(i).get(0, j);
        }
        cv::Mat curdata(1, m_DistancesPlus.at(0).size(), CV_32F, data_arr);
        data.push_back(curdata);
        size_plus++;
    }

    for (int i = static_cast<int>(m_DistancesPlus.size()); i < static_cast<int>(m_DistancesPlus.size() + m_DistancesMinus.size()); i++)
    {
        int it = i - size_plus;
        float data_arr[m_DistancesMinus.at(0).size()];

        labels_arr_vec.push_back(0);

        for (int j = 0; j < static_cast<int>(m_DistancesMinus.at(0).cols()); j++)
        {
            data_arr[j] = m_DistancesMinus.at(it).get(0, j);
        }
        cv::Mat curdata(1, m_DistancesPlus.at(0).size(), CV_32F, data_arr);
        data.push_back(curdata);
    }

    /* Shuffle Data */
    std::vector<int> seeds;
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

    /* Oversample the minority class to match the majority class size */
    int minority_size = cv::countNonZero(labels_shuffled == minority_class);
    int majority_size = cv::countNonZero(labels_shuffled == majority_class);

    if (minority_size < majority_size)
    {
        int oversample_factor = majority_size / minority_size;
        cv::Mat oversampled_minority_samples;
        cv::Mat oversampled_minority_labels;

        for (int i = 0; i < oversample_factor; i++)
        {
            oversampled_minority_samples.push_back(samples_shuffled);
            oversampled_minority_labels.push_back(labels_shuffled);
        }

        samples_shuffled.push_back(oversampled_minority_samples);
        labels_shuffled.push_back(oversampled_minority_labels);
    }

    /* Create Dataset and initialize Classifier */
    cv::Ptr<cv::ml::TrainData> m_traindata = cv::ml::TrainData::create(samples_shuffled, cv::ml::ROW_SAMPLE, labels_shuffled);

    statistic_model = cv::ml::RTrees::create();
    auto criteria = cv::TermCriteria();
    criteria.type = cv::TermCriteria::MAX_ITER;
    criteria.maxCount = 300;

    statistic_model->setMaxDepth(10);
    statistic_model->setMinSampleCount(2);
    statistic_model->setTruncatePrunedTree(false);
    statistic_model->setUse1SERule(false);
    statistic_model->setUseSurrogates(false);
    statistic_model->setTermCriteria(criteria);
    statistic_model->setCVFolds(1);

    /* Calculate weights - You can remove this part as it's no longer needed with oversampling */
    // int zerosgt = labels_arr_vec.rows - cv::countNonZero(labels_arr_vec);
    // int onesgt = cv::countNonZero(labels_arr_vec);
    // float plusval = labels_arr_vec.rows / (2.0 * onesgt );
    // float minusval = labels_arr_vec.rows / (2.0 * zerosgt );
    // float w[2] = {minusval, plusval};
    // cv::Mat newweight = cv::Mat(1, 2, CV_32F, w);
    // statistic_model->setPriors(newweight);

    /* Train Classifier */
    MITK_INFO << "Start Training";
    statistic_model->train(m_traindata);
}

std::vector<std::vector<unsigned int>>  StreamlineFeatureExtractor::Predict()
{
    std::vector<std::vector<unsigned int>> index_vec;
    /*Create Dataset as cv::Mat*/
    cv::Mat dataTest;
    for (unsigned int i = 0; i < myindex.size(); i++)
    {
        float data_arr[m_DistancesTest.at(0).size()];

        for (unsigned int j = 0; j < m_DistancesTest.at(myindex.at(0)).cols(); j++)
        {
            data_arr[j] = m_DistancesTest.at(myindex.at(i)).get(0, j);
        }
        cv::Mat curdata(1, m_DistancesTest.at(myindex.at(0)).size(), CV_32F, data_arr);
        dataTest.push_back(curdata);
    }

    std::vector<unsigned int> indexPrediction;
    std::vector<float> e(myindex.size());
    std::vector<int> pred(myindex.size());
    std::fill(pred.begin(), pred.end(), 0.0f);

    /*For every Sample/Streamline get Prediction and entropy (=based on counts of Random Forest)*/
    MITK_INFO << "Predicting on all cores";
    #pragma omp parallel for
    for (unsigned int i = 0; i < myindex.size(); i++)
    {
        int val = statistic_model->predict(dataTest.row(i));
        pred.at(i) = val;

        #pragma omp critical
        if (val == 1)
        {
            indexPrediction.push_back(myindex.at(i));
        }

        cv::Mat vote;
        statistic_model->getVotes(dataTest.row(i), vote, 0);
        e.at(i) = ( -(vote.at<int>(1,0)*1.0)/ (vote.at<int>(1,0)+vote.at<int>(1,1)) * log2((vote.at<int>(1,0)*1.0)/ (vote.at<int>(1,0)+vote.at<int>(1,1))) -
                    (vote.at<int>(1,1)*1.0)/ (vote.at<int>(1,0)+vote.at<int>(1,1))* log2((vote.at<int>(1,1)*1.0)/ (vote.at<int>(1,0)+vote.at<int>(1,1))));

        if (std::isnan(e.at(i)))
        {
            e.at(i) = 0;
        }
    }
    entropy_vector = e;
    MITK_INFO << "Done";
    MITK_INFO << "--------------";
    MITK_INFO << "Number of fibers in prediction: " << indexPrediction.size();
    MITK_INFO << "Entropy vector size:" << e.size();
    MITK_INFO << "--------------";

    int lengths = std::count_if(e.begin(), e.end(), [&](auto const& val) { return val >= 0.95; });
    if (lengths > 500)
    {
        lengths = 500;
    }
    else
    {
        lengths = 100;
    }

    std::vector<unsigned int> indexUnc = Sort(e, lengths, 0);
    std::vector<unsigned int> newidx;

    for (unsigned int i = 0; i < indexUnc.size(); i++)
    {
        newidx.push_back(myindex.at(indexUnc.at(i)));
    }

    index_vec.push_back(indexPrediction);
    index_vec.push_back(newidx);

    return index_vec;

}

std::vector<std::vector<unsigned int>>  StreamlineFeatureExtractor::GetDistanceData(float &value)
{
    /*Vector which saves Fibers to be labeled based on fft subset uncertainty*/
    std::vector<std::vector<unsigned int>> index_vec;

    /*Get index of most unertain data (lengths defines how many data is saved)*/

  int lengths = std::count_if(entropy_vector.begin(), entropy_vector.end(),[&](auto const& val){ return val >= value; });
  if (lengths>250)
  {
      lengths=250;
  }
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


    index_vec.push_back(indexUncDist);

    return index_vec;
}

}





 
