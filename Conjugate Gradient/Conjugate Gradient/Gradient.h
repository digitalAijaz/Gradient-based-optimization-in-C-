#ifndef GRADIENT_H
#define GRADIENT_H

#include <vector>
#include<functional>

class Gradient
{
public:
    Gradient();
    virtual ~Gradient();

    //Function to add a pointer to the object function
    void SetObjectFcn(std::function<double(std::vector<double>*)> objectFcn);

    //Function to set the start point
    //Note that this also sets the number of degrees of Freedom
    void SetStartPoint(const std::vector<double> startPoint);

    //Function to set the step size
    void SetStepSize(double stepSize);

    //Function to set the maximum number of iterations
    void SetMaxIterations(int maxIterations);

    //Function to set the gradient magnitude threshold (stopping condition)
    void SetGradientThresh(double gradientThresh);

    //Function to perform the optimization
    bool Optimize(std::vector<double>* funcLoc, double* funcVal);

    //Function to perform optimization using Fletcher and Reeves conjugate gradient method

    bool Fletcher_Reeves(std::vector<double>* funcLoc, double* funcVal);


private:

    //Function to calculate gradient in specified dimension (partial differentiation)
    double ComputeGradient(int dim);

    //Function to compute the gradient vector
    std::vector<double> ComputeGradientVector();

    //Function to compute the gradient magnitude
    double ComputeGradientMagnitude(std::vector<double> gradientVector);

    //Function to Calculate the step size
    double CalculateStepSize(double gradientMagnitude, std::vector<double> gradientVector, std::vector<double> m_currentPoint);

    //Function to Calculate the step size for Conjugate Gradient Using strong Wolfe conditions
    double CalculateStepSizeCG(double gradientMagnitude, std::vector<double> gradientVector, std::vector<double> m_currentPoint);

    //Matrix multiplication
    double vectorMultiplication(std::vector<double> Mat1, std::vector<double> Mat2);


private:

    //The number of degrees of freedom of the system
    int m_nDims;

    //The step size
    double m_stepSize;

    //The maximum number of iterations
    int m_maxIter;

    //The gradient step size (h)
    double m_h;

    //The gradient magnitude threshold (stopping condition)
    double m_gradientThresh;

    //The start point
    std::vector<double> m_startPoint;

    //The current point
    std::vector<double> m_currentPoint;

    //Function pointer to the object function
    std::function<double(std::vector<double>*)> m_objectFcn;
};

#endif // GRADIENT_H

