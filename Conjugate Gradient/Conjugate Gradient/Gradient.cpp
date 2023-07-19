#include "Gradient.h"
#include <iostream>
#include <fstream>
#include <math.h>

//Constructor
Gradient::Gradient()
{
    //Set defaults
    m_nDims = 0;
    m_stepSize = 0.0;
    m_maxIter = 1;
    m_h = 0.0001;
    m_gradientThresh = 1e-03;
}

//Destructor
Gradient::~Gradient()
{
    //dtor
}

//Function to set the object function

void Gradient::SetObjectFcn(std::function<double(std::vector<double>*)> objectFcn)
{
    m_objectFcn = objectFcn;
}

//Function to set the start point
//Note that this also sets the number of degrees of freedom
void Gradient::SetStartPoint(const std::vector<double> startPoint)
{
    //Copy the start point
    m_startPoint = startPoint;

    //Determine the number of degrees of freedom

    m_nDims = m_startPoint.size();
}

//Function to set the step size
void Gradient::SetStepSize(double stepSize)
{
    m_stepSize = stepSize;

}


//Function to Calculate the step size
double Gradient::CalculateStepSize(double gradientMagnitude, std::vector<double> gradientVector, std::vector<double> m_currentPoint)
{
    //Computing Step length
    gradientMagnitude = gradientMagnitude;
    gradientVector = gradientVector;

    std::vector<double> currentPoint = m_currentPoint;
    //  currentPoint = m_currentPoint;
    double funcValPrev = m_objectFcn(&currentPoint);

    std::vector<double> nextPoint = currentPoint;

    double funcValNext;

    double alpha{ 1.0 };
    double beta{ 0.1 };
    double gamma{ 0.6 };
    double armijo{};
    double goldstein{};

    /*   for (int i = 0; i < m_nDims; ++i)
       {
           nextPoint[i] = nextPoint[i] + (-(gradientVector[i] * alpha));
       }

       funcValNext = m_objectFcn(&nextPoint);
       armijo = (funcValPrev - (alpha * beta * gradientMagnitude * gradientMagnitude));

       while (funcValNext <= armijo)
       {
           alpha = 0.3 * alpha;
           for (int i = 0; i < m_nDims; ++i)
           {
               nextPoint[i] = nextPoint[i] + (-(gradientVector[i] * alpha));
           }

           funcValNext = m_objectFcn(&nextPoint);
           armijo = (funcValPrev - (alpha * beta * gradientMagnitude * gradientMagnitude));

       }
   */
    do
    {
        for (int i = 0; i < m_nDims; ++i)
        {
            nextPoint[i] = nextPoint[i] + (-(gradientVector[i] * alpha));
        }

        funcValNext = m_objectFcn(&nextPoint);
        //std::cout<<funcValNext<<"\n";

        armijo = (funcValPrev - (alpha * beta * gradientMagnitude * gradientMagnitude));
        //  armijo = (funcValPrev + (alpha * beta * (-1.0)));
        goldstein = (funcValPrev - (alpha * gamma * gradientMagnitude * gradientMagnitude));

        alpha = 0.3 * alpha;

    } while (funcValNext <= armijo);
    // while (funcValNext <= armijo && funcValNext >= goldstein);

    return alpha;

}


//Function to set the maximum number of iterations
void Gradient::SetMaxIterations(int maxIteration)
{
    m_maxIter = maxIteration;
}

//Function to set the gradient magnitude threshold (stopping condition)
//Optimization stops when gradient magnitude falls below this value
void Gradient::SetGradientThresh(double gradientThresh)
{
    m_gradientThresh = gradientThresh;
}

//Function to perform the actual optimization

bool Gradient::Optimize(std::vector<double>* funcLoc, double* funcVal)
{

    //Set the current point to start point
    m_currentPoint = m_startPoint;

    //Loop up to max iterations or until threshold reached
    int iterCount = 0; double stepSize = 1.0;
    double gradientMagnitude = 1.0;
    while ((iterCount < m_maxIter) && (gradientMagnitude > m_gradientThresh))
    {
        //Compute the gradient vector
        std::vector<double> gradientVector = ComputeGradientVector();
        gradientMagnitude = ComputeGradientMagnitude(gradientVector);

       

        stepSize = CalculateStepSize(gradientMagnitude, gradientVector, m_currentPoint);

        //Compute the new point
        std::vector<double> newPoint = m_currentPoint;
        for (int i = 0; i < m_nDims; ++i)
        {
            //  newPoint[i]  = newPoint[i] + (-(gradientVector[i]*m_stepSize));
            newPoint[i] = newPoint[i] + (-(gradientVector[i] * stepSize));

        }


        //Update the current point
        m_currentPoint = newPoint;

        //Increment the iteration counter
        iterCount++;
        std::cout << "Step Length: " << stepSize << "   " << iterCount << "    x=  " << m_currentPoint[0] << "    y=  " << m_currentPoint[1] << "      Gradient Magnitude   " << gradientMagnitude << "      Function Value   " << m_objectFcn(&m_currentPoint) << "\n";
    }
    std::cout << "\n\nNo of iterations: " << iterCount << "\n\n";
    //std::cout<<"\n\nNo of iterations: "<<iterCount<<"    Step Length: "<<stepSize<<"\n\n";
    //return the results
    *funcLoc = m_currentPoint;
    *funcVal = m_objectFcn(&m_currentPoint);

    return 0;
}


//Function to compute the gradient of the object function in the specified dimension
double Gradient::ComputeGradient(int dim)
{
    //make a copy of the current location
    std::vector<double> newPoint = m_currentPoint;

    //Modify the copy, according to h and dim
    newPoint[dim] += m_h;

    //Compute the two function values for these points
    double funcVal1 = m_objectFcn(&m_currentPoint);
    double funcVal2 = m_objectFcn(&newPoint);

    //Compute the approx numerical gradient
    return (funcVal2 - funcVal1) / m_h;
}

//Function to compute the gradient vector
std::vector<double> Gradient::ComputeGradientVector()
{
    std::vector<double> gradientVector = m_currentPoint;
    for (int i = 0; i < m_nDims; i++)
        gradientVector[i] = ComputeGradient(i);

    return gradientVector;
}

//Function to compute the gradient magnitude
double Gradient::ComputeGradientMagnitude(std::vector<double> gradientVector)
{
    double vectorMagnitude = 0.0;
    for (int i = 0; i < m_nDims; ++i)
        vectorMagnitude += gradientVector[i] * gradientVector[i];
    return sqrt(vectorMagnitude);
}


//Function to perform optimization using Fletcher and Reeves conjugate gradient method
bool Gradient::Fletcher_Reeves(std::vector<double>* funcLoc, double* funcVal)
{
   // std::cerr << "Start Fletcher_Reeves\n";
    //Set the current point to start point
    m_currentPoint = m_startPoint;

  /*  double funcValInit{};    // Variable for calculating initial function value f(0) given x(0) or start point
    funcValInit = m_objectFcn(&m_currentPoint);
  */

    double gradientMagnitude{};
    //Compute the initial gradient vector
    std::vector<double> gradientVector; 

    gradientVector = ComputeGradientVector();
    gradientMagnitude = ComputeGradientMagnitude(gradientVector);

    double stepSize{};

    std::vector<double> descentDir;

    descentDir.resize(m_nDims);
    for (int i = 0; i < m_nDims; ++i)
    {
        descentDir[i] = -1.0 * (gradientVector[i]);
    }

    //Loop until threshold reached
    int iterCount{0};
    
    while (gradientMagnitude > m_gradientThresh)

    {
   //     std::cerr << "Start While loop Fletcher_Reeves\n";
        std::vector<double> newPoint = m_currentPoint;

        stepSize = CalculateStepSizeCG(gradientMagnitude, gradientVector, m_currentPoint);

        //Compute the new point
        
        for (int i = 0; i < m_nDims; ++i)
        {
            newPoint[i] = newPoint[i] + (descentDir[i] * stepSize);

        }

        double beeta = 1 / (gradientMagnitude * gradientMagnitude);

        //Update the current point
        m_currentPoint = newPoint;
        gradientVector = ComputeGradientVector();
        gradientMagnitude = ComputeGradientMagnitude(gradientVector);

        beeta = beeta * gradientMagnitude * gradientMagnitude;

        for (int i = 0; i < m_nDims; ++i)
        {
            descentDir[i] = -1* gradientVector[i] + (beeta * descentDir[i]);

        }

        iterCount++;

    }

    std::cout << "\n\nNo of iterations: " << iterCount << "\n\n";
   
    //return the results
    *funcLoc = m_currentPoint;
    *funcVal = m_objectFcn(&m_currentPoint);

    return 0;
}

double Gradient::CalculateStepSizeCG(double gradientMagnitude, std::vector<double> gradientVector, std::vector<double> m_currentPoint)
{
 //   std::cerr << "inside CalculateStepSizeCG\n";
    //Computing Step length
    gradientMagnitude = gradientMagnitude;
    gradientVector = gradientVector;

    std::vector<double> currentPoint = m_currentPoint;
   
    double funcValPrev = m_objectFcn(&currentPoint);

    std::vector<double> nextPoint = currentPoint;

    double funcValNext;

    double alpha{ 0.5 };
    double C1{ 0.1 };   // 0 < C1 < C2 < 0.5
    double C2{ 0.4 };
    double armijo{0};
    double wolfeLHS{0};
    double wolfeRHS{0};

    std::vector<double> descentDir;
    descentDir.resize(m_nDims);
    for (int i = 0; i < m_nDims; ++i)
    {
        descentDir[i] = -1.0*(gradientVector[i]);
    }

    double prodArmijo{}; // Variable to calculate 2nd part of RHS product i.e. product of gradient transpose and descent dir and RHS of Strong Wolfe
    prodArmijo = vectorMultiplication(descentDir, gradientVector);



    do
    {
        for (int i = 0; i < m_nDims; ++i)
        {
            nextPoint[i] = nextPoint[i] + (descentDir[i] * alpha);
        }

        funcValNext = m_objectFcn(&nextPoint);

        //Update the current point
        m_currentPoint = nextPoint;

        //Compute the new gradient vector for Wolfe
        std::vector<double> gradientVectorNew = ComputeGradientVector();

        wolfeLHS = vectorMultiplication(gradientVectorNew, descentDir);
        wolfeRHS = -1.0*(C2 * prodArmijo);
     //   wolfeRHS =  (C2 * prodArmijo);

        armijo = (funcValPrev + (alpha * C1 * prodArmijo));

        alpha = 0.3 * alpha;

    } while (funcValNext <= armijo && wolfeLHS<= wolfeRHS);
  

    return alpha;
}

//Matrix multiplication
double Gradient::vectorMultiplication(std::vector<double> Mat1, std::vector<double> Mat2)
{
//    std::cerr << "inside vectorMultiplication\n";
    int product = 0;

    // Loop for calculate dot product
    for (int i = 0; i < m_nDims; ++i)
        product = product + Mat1[i] * Mat2[i];

    return product;
}